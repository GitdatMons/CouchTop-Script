import clr
import os
import re
import shutil
import sys

from connect import *  # Interact w/ RS
from pydicom import dcmread  # Read DICOM from file
from scripy.signal import find_peaks  #Find 'peak' values in an array

clr.AddReference('System.Windows.Forms')
from System.Windows.Forms import MessageBox   #For display errors


INNER_COUCH_NAME = 'Elakta Couch Foam Core'
OUTER_COUCH_NAME = 'Elekta Carbon Fiber Shell'
ADD_COUCH_TOP_EXPORT_PATH = os.path.join('T:', os.sep, 'Physics', 'Temp', 'Add Couch')


def add_couchtop_to_couch() -> None:
    '''Adds template couch top at correct location on the current exam.
    
    Correct location is R-L center, I-S center and on top of sim couch P-A
    Alert user of any couch-top geometries in the incorrect position that can't be moved because the geometries are approved
    Set default dose grid for any non-initial sim plan on the current exam. An initial sim plan name contains 'initial sim' or 'trial_1' (case insensitive). 
    Assumes patient position of current exam is HFS, HFP, or FFS
    
    Couch top is 51.4 cm wide, 196 cm long, and 2.8 cm tall 
    with a water equivalence of 2.0 mm at 6 MV.

    Most coding copied from scripts by Kalew White
    '''

    # Ensure a case is open
    try:
        case = get_current('Case')
    except:
        MessageBox.Show('No case is open. Click OK to abort the script.', 'No Open Case')
        sys.exit()
    
    # Ensure the case has an exam
    try:
        exam = get_current('Examination')
    except:
        MessageBox.Show('There are no exams in the open case. Clcik OK to abort the script', 'No Exams')
        sys.exit()
    
    patient_db = get_current('PatientDB')
    struct_set = case.PatientModel.StructuresSets[exam.Name]

    warnings = '' # We will display warnings at the end of the script

    # Ensure exam is not used in an approved beam set
    approved_beam_set_names = []
    for plan in case.TreatmentPlans:
        for beam_set in plan.BeamSets:
            if beam_set.GetPlanningExamination().Equals(exam) and beam_set.Review is not None and beam_set.Review.ApprovalStatus == 'Approved':
                approved_beam_set_names.append(beam_set.DicomPlanLabel)
    if approved_beam_set_names:
        MessageBox.Show('The couch geometries have materials defined, and the current exam is used in approved beam set(s) (' + ', '.join('"' + name + '"' for name in approved_beam_set_names) + '), so couch geometries cannot be added. Click OK to abort the script.', 'Exam in Approved Beam Set(s)')
        sys.exit()
    
    # Patient position
    pt_pos = exam.PatientPosition
    if pt_pos not in ['HFS', 'HFP', 'FFS']:
        MessageBox.Show('This script does not support ' + pt_pos + ' exams. Click OK to abort the script.', 'Unsupported Patient Position')
        sys.exit()
    is_supine = exam.PatientPosition in ['HFS', 'FFS']

    # Get external ROI
    try:
        ext_geom = next(geom for geom in struct_set.RoiGeometries if geom.OfRoi.Typ == 'External')
        couch_top = ext_geom.OfRoi
    except StopIteration:
        MessageBox.Show('There is no external geometry on the current exam. Click OK to abort the script', 'No External Geometry')
        sys.exit()

    # Will need to check whether various ROIs are approved
    approved_roi_names = [geom.OfRoi.Name for ss in struct_set.ApprovedStructureSets for geom in ss.ApprovedRoiStructures]

    # Density values are cleared and dose values are invalidated, so recompute dose on unapproved beam sets that have dose and that are on the current examination
    to_recompute = [bs for p in case.TreatmentPlans for bs in p.BeamSets if bs.GetPlanningExamination().Equals(exam) and (p.Review is None or p.Review.ApprovalStatus != 'Approved') and bs.FractionDose.DoseValues is not None and bs.FractionDose.DoseValues.IsClinical]

    # Ensure Couch top is approved
    if couch_top.Name in approved_roi_names:
        MessageBox.Show('Couch Top is approved on the current exam, so the effects of it cannot be changed. Click OK to abort the script', 'Couch Top Is Approved')
        sys.exit()
    
    # Get couch ROIs
    try:
        inner_couch = struct_set.RoiGeometries[INNER_COUCH_NAME]
        outer_couch = struct_set.RoiGeometries[OUTER_COUCH_NAME]
    except:
        MessageBox.Show('One or both of ' + INNER_COUCH_NAME + ' and ' + OUTER_COUCH_NAME + ' ROIs do not exist in the current case. Clcik OK to abort the script.', 'Missing Couch ROI(s)')
        sys.exit()

    # Couch template and structure names
    template_name = 'Elekta Couch' if is_supine else 'Elekta Prone Couch'
    template = patient_db.LoadTemplatePatientModel(templateName=template_name)
    couch_structs = template.PatientModel.RegionsOfInterest
    outer_couch, inner_couch = couch_structs

    # Get couch geometries
    if not inner_couch.HasContours() or not outer_couch.HasContours():
        MessageBox.Show( 'One or both of ' + INNER_COUCH_NAME + ' and ' + OUTER_COUCH_NAME + ' geometries do not exist on the current exam. Click OK to abort the script.', 'Missing Couch Geometry(ies)')
        sys.exit()

    # Apply template
    # Doesn't matter which source exam we use, so just use the first
    case.PatientModel.CreateStructuresFromTemplate(SourceTemplate=template, SourceExaminationName=template.StructureSetExaminations[0].Name, SourceRoiNames=[outer_couch.Name, inner_couch.Name], SourcePoiNames=[], TargetExamination=exam, AssociateStructuresByName=True, InitializationOption='AlignImageCenters')

    ## Export exam so we can use pydicom to extract pixel data
    # Create export folder
    if os.path.isdir(ADD_COUCH_EXPORT_PATH):
        shutil.rmtree(ADD_COUCH_EXPORT_PATH)
    os.makedirs(ADD_COUCH_EXPORT_PATH)

    # Export exam
    # No RS functionality to export a single CT slice
    # Save Patient
    patient.Save()  # Error if you don't save before export
    try:
        case.ScriptableDicomExport(ExportFolderPath=ADD_COUCH_EXPORT_PATH, Examinations=[exam.Name], IgnorePreConditionWarnings=False)
    except:
        case.ScriptableDicomExport(ExportFolderPath=ADD_COUCH_EXPORT_PATH, Examinations=[exam.Name], IgnorePreConditionWarnings=True)


  
     ## Get y-coordinate of top of sim couch
    # If entire sim couch is in scan, the couch top is the fourth 'bright' peak in the R-L middle

    # Get pixel intensities array
    first_filename = os.listdir(ADD_COUCH_EXPORT_PATH)[0]  # Doesn't matter which CT DICOM file we use, so just use the first one
    dcm = dcmread(os.path.join(ADD_COUCH_EXPORT_PATH, first_filename))
    intensities = dcm.pixel_array

    # Select intensities in R-L center
    y_axis_intensities = [row[intensities.shape[0] // 2] for row in intensities]

    # Select 'bright peaks' on y-axis
    # These are white areas on image
    # The `prominence` argument was obtained through trial and error
    # `find_peaks` returns a tuple whose first element is the indices of the peaks
    y_axis_peaks = find_peaks(y_axis_intensities, height=(0, 700), prominence=180)[0]

    if len(y_axis_peaks) > 0:
        # 2nd peak from the posterior end
        second_peak = y_axis_peaks[-1]

        # Convert pixel y-coordinate to exam coordinate
        if is_supine:
            correct_y = img_stack.Corner.y + second_peak * img_stack.PixelSize.y - 2.25
        else:
            correct_y = img_stack.Corner.y - second_peak * img_stack.PixelSize.y + 2.25
        
        # The TOP of the couch goes at the 2nd peak, so compute where the P-A CENTER of the couch should go
        # Center is half a couch height down
        outer_couch_bounds = struct_set.RoiGeometries[outer_couch.Name].GetBoundingBox()
        outer_couch_ht = outer_couch_bounds[1].y - outer_couch_bounds[0].y
        if is_supine:
            correct_y += outer_couch_ht / 2
        else:
            correct_y -= outer_couch_ht / 2
    else:
        correct_y = None
        warnings += 'Could not determine correct AP coordinate. Manually position the couch geometries in the transverse view and then run the Center Couch script to correct any RL error you may have introduced.'

    # Correct x- and z-coordinates
    # x is R-L center (always zero)
    # z is I-S center + 2.8 cm
    correct_x = 0
    img_bounds = img_stack.GetBoundingBox()
    correct_z = ((img_bounds[0].z + img_bounds[1].z) / 2) + 2.8

    # Center each couch geometry
    for roi in (inner_couch, outer_couch):  # Iterate over couch ROIs in the applied couch template
        geom = struct_set.RoiGeometries[roi.Name]  # That couch ROI's geometry on the exam

        # Current position
        geom_ctr = geom.GetCenterOfRoi()

        # If AP coordinate (top of couch) was not found above, don't move the couch in the AP direction
        y_chg = 0 if correct_y is None else correct_y - geom_ctr.y

        # Transformation matrix
        # Each x, y, and z transform is the difference between the correct and current coordinates
        mat = {'M11': 1, 'M12': 0, 'M13': 0, 'M14': correct_x - geom_ctr.x,
               'M21': 0, 'M22': 1, 'M23': 0, 'M24': y_chg,
               'M31': 0, 'M32': 0, 'M33': 1, 'M34': correct_z - geom_ctr.z, 
               'M41': 0, 'M42': 0, 'M43': 0, 'M44': 1}

        # Reposition couch structure
        geom.OfRoi.TransformROI3D(Examination=exam, TransformationMatrix=mat)

        # Show structure
        patient.SetRoiVisibility(RoiName=roi.Name, IsVisible=True)

    # For any plan on this exam that has no dose, set default dose grid
    for p in case.TreatmentPlans:
        if not re.search(r'(initial )?sim|trial.?1', p.Name, re.IGNORECASE) and p.GetTotalDoseStructureSet().OnExamination.Equals(exam) and p.TreatmentCourse.TotalDose.DoseValues is None:
            voxel_sz = p.GetTotalDoseGrid().VoxelSize
            if any(val == 0 for val in voxel_sz.values()):
                voxel_sz = {coord: 0.3 for coord in 'xyz'}
            for bs in p.BeamSets:
                bs.SetDefaultDoseGrid(VoxelSize=voxel_sz)
                bs.FractionDose.UpdateDoseGridStructures()

    # Recompute invalidated beam sets
    for bs in to_recompute:
        bs.ComputeDose(ComputeBeamDoses=True, DoseAlgorithm=bs.AccurateDoseAlgorithm.DoseAlgorithm)

    # Remove unnecessary exported files
    shutil.rmtree(ADD_COUCH_EXPORT_PATH)

    # Display warnings if they exist
    if warnings:
        MessageBox.Show(warnings, 'Warnings')
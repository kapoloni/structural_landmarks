#!/usr/bin/python3
import sys


def showUsage():
    print('Usage: ./atlas_sides.py')
    sys.exit(2)
# ********************************************************************
#
# ********************************************************************


if __name__ == "__main__":
    if len(sys.argv) < 1:
        showUsage()
        exit()

    bof = "../../filter/logGab_#_256_256_256_04_06_03_033_209_055_120_0.bof"
    template_image = "../../../template/NAC_T1_RAI.nii.gz"
    nr, na, rmax = ['5', '10', '32']
    desc = '-getl'
    auxdata_path = "../../atlas_test/sorted"
    pc_folder = "../../../images/pc/ixi"
    input_image_path = "/databases/data2/IXI/brain_extraction/"
    input_mask = "/databases/data2/IXI/mesh/r8/"
    output_path = "output"
    output_age_path = output_path+"/IXI_20-80"
    hospital_names = ["Guys", "HH", "IOP"]
    atlas_samples_prefix = "landmarks_atlas_sample"

    sys.path.append("../../")
    from atlas_utils import Utils

    utils = Utils(na, nr, rmax, desc, output_path, bof)
    utils.create_folders(pc_folder)
    utils.create_folders(output_path)

    path_in, path_mask, path_out = utils.get_directories_mask(input_image_path, input_mask, pc_folder)
    path_in = [pt+".nii.gz" for pt in path_in]
    path_mask = [pt+"_RH.nii.gz" for pt in path_mask]
    path_out = [pt for pt in path_out]

    # Check PC
    utils.multip(path_in, path_mask, path_out, utils.phase_congruency, 60)
    # Landmarks
    utils.multip(path_in, path_mask, path_out, utils.landmarks, 20)
    # Atlas construction
    utils.probabilistic_atlas_landmarks(output_age_path,
                                        utils.get_input_paths_atlas(output_age_path, ["Guys", "HH", "IOP"],
                                                                    "../../atlas_test/sorted", atlas_samples_prefix),
                                        atlas_samples_prefix, template_image)
    print("Done")

For IXI images
# Phase congruency
./pc_ixi.py -r /databases/data1_study1/AD/MRI/Katia/IXI/brain_extraction/
            -n ../images/pc/IXI

# IXI is used only for the atlas


# For ADNI images:
    - PC
    - Fast
    - landmark detections
    - landmark descritors extractions

# Phase congruency
./pc_adni.py -r /databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/brain_extraction
             -n ../images/pc/ADNI

# Matching descriptor
    # Spider-web
    ./land_adni.py -r /databases/data2/ADNI/images/intensity_bs
                -p ../../images/hippocampus/adni/pc
                -l ../../images/hippocampus/adni/land/r8 -t s

# Fast ADNI

# SVM attr descriptor
# Log-spherical
./land_adni_tissues.py -r /databases/data2/ADNI/images/intensity_bs
                       -p ../../images/hippocampus/adni/pc
                       -l ../../images/hippocampus/landcsf/32 -t csf

./land_adni_tissues.py -r /databases/data2/ADNI/images/intensity_bs
                       -p ../../images/hippocampus/adni/pc
                       -l ../../images/hippocampus/landgm/32 -t gm

./land_adni_tissues.py -r /databases/data2/ADNI/images/intensity_bs
                       -p ../../images/hippocampus/adni/pc
                       -l ../../images/hippocampus/landwm/32 -t wm

# Join Log-spherical descriptors
./join_descriptors.py -n ../../images/hippocampus/descs_tissues
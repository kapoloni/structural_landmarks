# For IXI images:
### IXI is used only for the atlas
#### Phase congruency

    ./pc_ixi.py --ref /databases/data1_study1/AD/MRI/Katia/IXI/brain_extraction/
                --new_folder ../images/pc/IXI
# For ADNI images:

    - PC
    - Fast
    - landmark detections
    - landmark descritors extractions

#### Phase congruency

     ./pc_adni.py --ref /databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/brain_extraction
                  --new_folder ../images/pc/ADNI

#### Matching descriptor

    # Spider-web
    ./landmark_detection.py --ref /databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/brain_extraction
                            --pc ../images/pc/ADNI
                            --land ../images/hippocampus/descriptor/spiderweb -type sw

#### SVM attr descriptor
#### Log-spherical space

    ./landmark_detection_prob_maps.py --ref /databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/brain_extraction
                                      --pc ../images/pc/ADNI
                                      --land ../images/hippocampus/descriptor/csf -tissue csf

    ./landmark_detection_prob_maps.py --ref /databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/brain_extraction
                                      --pc ../images/pc/ADNI
                                      --land ../images/hippocampus/descriptor/gm -tissue gm

    ./landmark_detection_prob_maps.py --ref /databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/brain_extraction
                                      --pc ../images/pc/ADNI
                                      --land ../images/hippocampus/descriptor/wm -tissue wm

#### Join Log-spherical descriptors
     
     ./join_descriptors.py --new_folder ../images/hippocampus/descriptor/tissues

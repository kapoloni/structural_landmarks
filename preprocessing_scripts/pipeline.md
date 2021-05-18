# For IXI images:
### IXI is used only for the atlas
#### Phase congruency

    ./pc_ixi.py -r /databases/data1_study1/AD/MRI/Katia/IXI/brain_extraction/
                -n ../images/pc/IXI
# For ADNI images:

    - PC
    - Fast
    - landmark detections
    - landmark descritors extractions

#### Phase congruency

     ./pc_adni.py -r /databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/brain_extraction
                  -n ../images/pc/ADNI

#### Matching descriptor

    # Spider-web
    ./land_adni.py -r /databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/brain_extraction
                   -p ../images/pc/ADNI
                   -l ../images/hippocampus/descriptor/spiderweb -t s

#### SVM attr descriptor
#### Log-spherical space

    ./land_adni_tissues.py -r /databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/brain_extraction
                           -p ../images/pc/ADNI
                           -l ../images/hippocampus/descriptor/csf -t csf

    ./land_adni_tissues.py -r /databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/brain_extraction
                           -p ../images/pc/ADNI
                           -l ../images/hippocampus/descriptor/gm -t gm

    ./land_adni_tissues.py -r /databases/data1_study1/AD/MRI/Katia/ADNI/experiment/ADNI/brain_extraction
                           -p ../images/pc/ADNI
                           -l ../images/hippocampus/descriptor/wm -t wm

#### Join Log-spherical descriptors
     
     ./join_descriptors.py -n ../images/hippocampus/descriptor/tissues

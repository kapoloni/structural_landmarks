# Implementation of the paper
## Automated detection, selection and classification of hippocampal landmark points for the diagnosis of Alzheimer's disease.

Our approach, based on a two-level classification, first detects and
selects discriminative landmark points from two diagnosis populations based on
their matching distance compared to a probabilistic atlas of 3-D labeled landmark points.
The points are classified using attributes computed in a spherical
support region around each point using information from brain probability image
tissues of gray matter, white matter, and cerebrospinal fluid as sources of
information. Next, at the second level, the images are classified based on a
quantitative evaluation obtained from the first-level classifier outputs.

<img src="https://github.com/kapoloni/structural_landmarks/blob/main/figures/Diagram_landmarks_level.png" width=70%/>

Images from the public ADNI and IXI datasets:<br>
http://adni.loni.usc.edu/ and https://brain-development.org/ixi-dataset/<br>

# This repository was tested using:

    Python                 3.6.9
    C++                    7.5.0
    ITK                    5.1

# This repository contains:
    Images filename used in each fold in experiment/cmpb/splits
 <img src="https://github.com/kapoloni/structural_landmarks/blob/main/figures/cross_validation.png" width=50%/>

    preprocessing scripts
    landmark level scripts
    image level scripts    
    c++ src code

# Results
## Best image-level results

 <table style="width:60%">
  <tr>
    <th></th>  
    <th>CNxMCI</th>
    <th>MCIxAD</th>
    <th>CNxAD</th>
  </tr>
  <tr>
    <td>AUC</td>
    <td>0.83 ± 0.05</td>
    <td>0.73 ± 0.07</td>
    <td>0.95 ± 0.03</td>
  </tr>
  <tr>
    <td>Accuracy</td>
    <td>75.58 ± 3.6</td>
    <td>69.8 ± 6.73</td>
    <td>89.24 ± 4.04</td>
  </tr>
</table> 

## Scatter plots of four (out of the ten) image attributes for one fold of the CN×MCI, MCI×AD, and CN×AD experiments
<table style="width:60%">
  <tr>    
    <td><img src="https://github.com/kapoloni/structural_landmarks/blob/main/figures/cn_mci_points.png" width=80%/>CN×MCI</td> 
    <td><img src="https://github.com/kapoloni/structural_landmarks/blob/main/figures/mci_ad_points.png" width=80%/>MCI×AD</td> 
  </tr>
  <tr>
    <td> <img src="https://github.com/kapoloni/structural_landmarks/blob/main/figures/cn_ad_points.png" width=80%/>CN×AD</td> 
  </tr>
</table>
      


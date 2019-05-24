# Registrar

Copyright 2018 Northeastern University
</br>
This project has been supported by the National Institude of Health (NIH)
</br></br>
Registrar is open source software for accurate spatial registration of multiple overlapping stacks of images, registration of stacks of images acquired in a time-lapse manner, and registration of image plains within individual stacks. Registrar provides a utility for registration based on translation, rigid, affine, and B-spline transformations. This software was developed by Seyed M.M. Kahaki and Armen Stepanyants with input from other members of the Neurogeometry group.
</br></br>
<img src="https://web.northeastern.edu/kahaki/Registrar_.PNG" alt="Registrar" align="middle"> 
</br></br></br>
<img src="https://web.northeastern.edu/kahaki/reg_before_after.PNG" alt="Registrar" align="middle">

# Functionalities
The ability to map neural circuits on the scale of an entire brain is critical for advancing our understanding of brain functions. Circuit mapping can be based on whole-brain imaging of sparsely labeled populations of neurons with 3D confocal or two-photon microscopy. Registrar is a software for stitching and registration designed for different proposes including 1. Stitching of 3D image stacks of entire brain, 2. registration of time-lapse images, and 3. registration of images within stacks. 
1.	Stitching of 3D image stacks of entire brain: whole brain imaging experiment, if applied to the mouse brain, would result in tens of thousands of image stacks, totaling several terabytes of data. Because, imaging is generally done with small overlaps between neighboring stacks, the information contained in the stack overlap regions can be used for stitching. The Registrar software can provide the ability of stitching 3D stacks in both sequential and parallel and output the registered stack positions along with the transformation for each stack for further analysis. 

2.	Registration of time-lapse images: the second goal of this application is to register 3D stack images captured in different times. These images are usually having high overlap because they are taken from a same part of the brain in different time. The current Registrar software can also register this type of images with a high accuracy. 

3.	Registration of images within stacks: the Registrar software can register stack slices within a single stack which is an important task in Neuroscience and many other fields. Large scale image stacks captured in live animal can be transformed within the slices of the stacks. This is a common issue specially in EM data. Using this ability of the software, you can register the slices based on different transformation provided in the GUI.

# Installation

Registrar can be downloaded at https://github.com/neurogeometry/registrar. Start MATLAB, navigate to the software folder, and run Registrar in MATLAB command window. Registrar is designed to run sequentially or in parallel, on a PC or a cluster.

# Sample Data

To test the software, download <a href="http://www.northeastern.edu/neurogeometry/wp-content/uploads/RegistrarSampleData.zip">RegistrarSampleData.zip</a> (~160 MB) and extract its components. This folder includes the Neocortical Layer 1 Axons dataset of the DIADEM Competition, which can be used along with spatial registration.

# A sample input CSV file content

% Stack location on disk, stack position x, stack position y, stack position z

E:\SubTiling\00735\00735-ngc.0.tif,203299.775250213,196643.100746094,19118.6238160000

E:\SubTiling\00736\00736-ngc.0.tif,204164.694690563,196643.100746094,19118.6238160000

E:\SubTiling\00759\00759-ngc.0.tif,203299.775250213,198088.934353151,19118.6238160000

E:\SubTiling\00760\00760-ngc.0.tif,204164.694690563,198088.934353151,19118.6238160000

# Code Structure

The Registrar.m is the main file which call the GUI and will other functions. The main function in this file called registration():

```
try
    registeration(StackList_csv_pth,TransformationValue,Seq_Par,Par_workers,blendingSID,handles,LogHandle)
catch ME
    LogHandle.Children(2).String = ME.getReport;
end
```
The registration() function includes <a href='https://github.com/neurogeometry/'>Stack Overlaps</a>, <a href='https://github.com/neurogeometry/'>Feature Extraction</a>, <a href='https://github.com/neurogeometry/'>Feature Matching</a>, <a href='https://github.com/neurogeometry/'>Blending</a>, <a href='https://github.com/neurogeometry/'>Retiling</a>, and <a href='https://github.com/neurogeometry/'>Create Zoom-levels</a> functions.

The detailed information about these functions can be found <a href='https://github.com/neurogeometry/'>here</a>.

# License

Copyright 2018 Northeastern University.


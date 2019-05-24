# Registrar Functions

The main functions in Registrar software are:
-<a href='https://github.com/neurogeometry/'>Registration()</a> (generates the stack overlaps information)</br>
-<a href='https://github.com/neurogeometry/'>Stack Overlaps()</a> (generates the stack overlaps information)</br>
-<a href='https://github.com/neurogeometry/'>Feature Extraction()</a> (find features in overlaps regions of all stacks)</br>
-<a href='https://github.com/neurogeometry/'>Feature Matching()</a> (pairwise feature matching) </br>
-<a href='https://github.com/neurogeometry/'>Global Registration()</a> (minimize the registration error globally and find transformation for each stack) </br>
-<a href='https://github.com/neurogeometry/'>Blending()</a> (generate seamless registered image for visualization)</br>
-<a href='https://github.com/neurogeometry/'>Retiling()</a> (generate final tiles)</br>
-<a href='https://github.com/neurogeometry/'>CreateZoomLevels()</a> (create different zoom levels)</br>

## Registration()
This is the main function in the GUI file (Registrar.m) which runs when user set parameters and click on Run on the GUI. 
```
try
    registeration(StackList_csv_pth,TransformationValue,Seq_Par,Par_workers,blendingSID,handles,LogHandle)
catch ME
    LogHandle.Children(2).String = ME.getReport;
end
```
The inputs of this function includes:

- StackList_csv_pth : the path to the input <a href='https://github.com/neurogeometry/Registrar#a-sample-input-csv-file-content'>CSV file</a>. 
- TransformationValue: this parameter defines the transformation type, 1 for Translation, 2 for Rigid, 3 for Affine, and 4 for B-Spline.
- Seq_Par: 1 for run in sequential, and 2 for parallel 
- Par_workers: sets the number of workers to be used from the local or remote cluster if Seq_Par = 1.
- blendingSID: the stack ID to be shown along with its neighboring stacks as a sample of the registration result.
- handles: the handles object of the main GUI.
- LogHandle: the handles object of the Log GUI.

## Stack Overlaps()

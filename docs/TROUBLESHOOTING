#KiT troubleshooting and diagnostics

These will be demonstrated using a test dataset of a HeLa cell imaged at a resolution of 2s per stack on a spinning disk microscope. To download this dataset (~100MB), run:
```
kitDataTest()
```

Track this movie via 
```
kitRun
```
using the default options selected in the GUI. 

Load the relevant jobset and job via:
```
jobset = kitLoadJobset();
job = kitLoadJob(jobset,1);
```

##Diagnostic plots

```
kitShowDragonTails
```
Use this to visualise tracks on an image. Helps to verify that a reasonable number of spots have been detected and tracked. Can identify tracks with label switching between kinetochores. 

Example: `kitShowDragonTails(job,’usePairs’,0,’crop',0,'timePoint',10,'tailLength',5);` Note that in the test dataset centrioles have also been labelled and have been tracked. We can choose to visualise only paired sister kinetochores: `kitShowDragonTails(job,'usePairs',1,'crop',0,‘timePoint',50,'tailLength',5,'tracks',1:60);`


```
kitShowSpots
```
Use this to visualise detected spots in a frame. A cropped region around each spot is shown. Helpful to check that spots are kinetochores as opposed to background. Biorientated kinetochores are expected to have a Gaussian shaped spot-like signal. 
Example: `kitShowSpots(job,’timePoint’,10)`

```
kitShowSisters
```
Use this to visualise paired sister kinetochores in a frame.
Example: `kitShowSisters(job,'timePoint',10)`


```
kitPlotHowManyTracksOfGivenLength
```
Use this to visualise a barcode plot of when tracks start and stop throughout a movie. This can highlight whether high quality tracking has been achieved, whether there is a certain part of the movie whether tracking is worse than elsewhere, or whether there is a constant gain and loss of tracks throughout the movie.
Example: `kitPlotHowManyTracksOfGivenLength(job)`

```
kitDiagnoseAndTest
```
Evaluate each step of the detection and tracking algorithm to show how many kinetochore sisters are detected, tracked, and grouped at each stage. 
This is helpful to identify drops in performance at a certain step which can indicate how to therefore improve tracking performance. 
Example: `kitDiagnoseAndTest(job,1,0)`
Note that since the original spot detections before the refinement step (eg. refinement of locations wiith a Gaussian mixture model) are not stored in the jobset object. Therefore for this evaluation the initial detections need to be performed again which can be time consuming.
It may be useful to obtain and store these via `[spots, movie, ~] = kitExampleDetectionTest(job);` if repeating this evaluation. The stored spot locations can then be used via `kitDiagnoseAndTest(job,1,0,spots,movie,1)`  

![performance](images/Performance.png)
![num sisters plot](images/NumSisters.png)
![barcode plot](images/Barcode.png)





KiT (Kinetochore Tracking) software
===================================

KiT is a MATLAB program for tracking kinetochores in microscopy movies in 2D or
3D, and in multiple channels.

This software is licensed for use under the GNU GPL v3 license. See LICENSE file
for details.

Installation
------------

No installation is necessary, simply copy the source code to the directory of
your choice. Enter that directory within MATLAB (use directory toolbar at top,
or 'cd' command) and then enter commands.

Usage
-----

An tracking experiment is organized in KiT as a 'jobset' stored in a
mat-file. To create a jobset:

   jobset = kitGUI

This will display a GUI to allow modification of tracking parameters. Click
'Select directory' to choose a directory to search for movie files. Select the
movies to analyse by highlighting them in the list and click Add ROI. For each
movie, select one or more ROIs around the cell(s) to be tracked. Double click
after drawing each ROI to save it. After ROI selection, configure tracking
parameters. Typically the default tracking parameters will be fine, but may need
to be changed on a case-by-case basis. Choose a name for the jobset. Finally,
click either 'Execute' to run tracking immediately, or 'Save' to run later.

To run tracking on the jobset manually:

   kitRunJobs(jobset);

or, to run as parallel batch jobs:

   kitRunJobs(jobset,'parallel',1);

or, to run on just a subset of movies, e.g. 2, 4 and 6:

   kitRunJobs(jobset,'subset',[2 4 6]);

After tracking is complete you will find a file named something like
'kittracking001_expname_moviefile.mat', where expname is replaced with the name of
jobset mat-file and moviefile is replaced with the name of the movie this result
is associated with.

To load a previously created jobset mat-file use:

   jobset = kitLoadJobset;

To save a jobset to its mat-file after making changes:

   kitSaveJobset(jobset);

To load a job's results from the 'kittracking*.mat' file, e.g. movie #2 of a
jobset:

   job = kitLoadJob(jobset,2);


Many functions within KiT have optional parameters. The 'help' command is your
friend here...


Quick analysis and diagnostics
------------------------------

Assuming you've loaded a jobset, as described above, you can quickly examine how
well tracking worked and generate some basic plots using the analysis GUI:

   kitAnalysis(jobset)

Advanced usage
--------------

Multiple jobsets may be analyzed with one command:
   kitMultiJobsetRun(filenames);
where filenames is a cell array for jobset filenames.

Quick movies showing tracked spots can be rendered using:
   kitMakeSpotMovie(job);

There are many options for this function (see 'help kitMakeSpotMovie'). In
particular try setting rotate, slow and plotPlane, e.g.
   kitMakeSpotMovie(job,'rotate',1,'slow',0.1,'plotPlane',1);

There are more options to control tracking than those available in the GUI. Look
at jobset.options, but you will have to inspect the code comments for
documentation.


Bugs
----

Send bug or crash reports and feature requests to jonathan.armond@warwick.ac.uk.

Include:
  - Brief description of steps to reproduce error.
  - Copy of the error text.
  - Jobset mat-file (if available).

Movies will be useful but should be sent via e.g. Dropbox or Files@Warwick.


Credits
-------

- KiT evolved from both Maki (Jaqaman, et al., JCB 188:665-679 (2010)) and
  u-Track (Jaqaman, et al., Nat. Meth. 5:695-702 (2008)). Large portions of the
  core tracking code is lifted directly from Maki and u-Track.

- E. Harry modified the u-Track mixture-model fitting to operate on 3D data,
  plus various other tweaks throughout.

- E. Vladimirou made modifications to the coordinate system fitting for
  non-congressed kinetochores, plus various other tweaks throughout.

- C. A. Smith contributed some bugfixes and additions to some of the plotting code.

- J. W. Armond organised the code into its current form as KiT. He added the
  adaptive wavelet spot finder, image-moment based coordinate system, spot intensity
  measurement, multiple channel tracking, and movie rendering. He also made a
  number of improvements toward code cleanliness and efficiency.


The KiT package contains code from various other sources. See file headers to
identity sources and licenses.


MATLAB requirements
--------------------

MATLAB R2014b
MATLAB Image Processing Toolbox 9.1
MATLAB Statistics Toolbox 9.1

Optional toolboxes
------------------

For parallel execution of tracking
    MATLAB Parallel Computing Toolbox 6.5
For adaptive threshold particle detection
    MATLAB Global Optimization Toolbox 3.3
For correction of photobleaching in intensity data and enhancing
adaptice particle finding algorithm
    MATLAB Curve Fitting Toolbox 3.5

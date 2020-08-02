### Background:
Parallel simulation using HPC cluster is very common in R&D activities. Often is it observed that many desired functionalities (e.g. for post-processing)
are missing in the used simulation software especially if it is an open sourced one. Similarly, an example of various output files are given here that come from
different cores of a parallel simulation. In practical case, the real problems will be run probably with some thousands of cores in supercomputer.
Here a simple serial code snippet is presented that read all output file, extract desired values at each time step and then combined the values and plotted.

### Problem description:
A very simple, parallel thermo-mechanical simulation results files are attached. These files start with Big O. For example, O10grains_0001 is the result file coming
from core-1, O10grains_0002 is from core-2 and so on. The task is to extract a particular component of stress and strain at each time step
and will make average over all 'O' files located in the same folder (e.g. 8 Ofiles here). Save the results in a text file and plot stress strain curve.

#### averageStressStrain.py:
this is the main python code (serial version) that reads all Ofiles, extracts stress and strain components along 33 direction at each time step and save the result into the text file 'avgStrainStress.txt'.
The result is plotted with matplotlib and shown in the .PNG file.

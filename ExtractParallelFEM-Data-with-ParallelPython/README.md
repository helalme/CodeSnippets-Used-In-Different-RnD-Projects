### Background:
Parallel FEM simulation using HPC cluster is very common in R&D activities. Often is it observed that many desired functionalities (e.g. for post-processing)
are missing in the used simulation software especially if it is an open sourced one. Here an example of 4 output files are given that come from
different cores of a parallel FEM simulation. In practical case, the real problems will be run probably with some thousands of cores in a supercomputer.
If a parallel simulation takes 1 hour, and its postprocessing with serial python also takes nearly 5 hour, then it would be very time comsuming. An alternate option here wouldbe
would be to parallelize the postprocessing Python-script. If it is parallelized properly, postprocessing time of 5 hours is reduced to few minutes.

### Problem description:
A very simple, parallel, thermo-mechanical simulation-results are attached (4 files). These files start with Big O. For example, O10grains_0001 is the result file coming
from core-1, O10grains_0002 is from core-2 and so on. The task is to extract a particular component of stress or strain at each time step. At each time step there might be
thousands of elements, and several Gauss-points inside each elements. To extract properly, it is necessary to extract a particular components from each Gauss-point. 
The task is to extract all necessay data and make average over all 'O' files located in the same folder (e.g. 4 Ofiles here), save the results in a text file.

#### AvgS33FromGPs-Parallel.py:
this is the parallel python code using multiprocessing package, that reads all Ofiles, extracts stress components along 33 direction from each time step, each element and from
each integration point and save the result into the text file 'StressAvgFromGPs-Parallel.txt'. **This code has been tested with several cores using 2 TB of simulation data.**

##### How to Use:
Simply run the file. Filename is to be changed if necessary.

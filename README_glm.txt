README_glm.txt
J Pillow
Created: Feb 1, 2010 

DESCRIPTION OF THE CODE AND INSTRUCTIONS FOR USE

Performs simulation and maximum-likelihood fitting for several
versions of a generalized linear point process model (GLM) for neural
spike train data.


INSTALLATION
============
(1) Unnpack the archive using a compression utility such as WinZip
    (or, from the command line in linux: 'tar -xvzf FILENAME.tgz').
(2) Launch matlab and cd into the main directory containing the code
    (e.g. '/code_GLM/').
(3) CD to the sub-directory 'tools_mexcode/' and run the initialization
    script 'initialize_mexcode' from the command line.  This will
    compile the C files in the sub-directory 'mexcode' using the
    default compiler settings

USE 
===
(1) From the main code directory, run "setpaths" to add relevant
    sub-directories to the matlab path and initialize the global
    variable ("RefreshRate").
(2) Examine scripts in sub-directory "testscripts/" to see simple
    examples of simulation and fitting to spike data using the code:

Test Scripts:
-------------
1. testscript_GLM.m - fits plain GLM (temporal stimulus kernel only).
2. testscript_GLM_coupled.m - fits GLM with two coupled neurons.m
3. testscript_GLM_spatialStim.m - fits GLM using two different parametrizations
     of stimulus kernel (linear vs. bilinear) 
4. testscript_GLM_splineNlin.m - fits GLM incorporating an arbitrary
     nonlinearity parametrized by cubic splines.
     - uses coordinate ascent of filter and nonlinearity params (not concave!)
     - illustrates code for plotting spike rasters


NOTES
=====
- time is represented in units of "stimulus frames", which is
controlled by the global variable RefreshRate (Hz).  Thus, for
example, if RefreshRate=100, each unit of time is 10ms, and a
post-spike kernel discretized in time bins of width dt=.02 has time
bins of length 0.2 ms.

- fitting code relies on the matlab optimization toolbox ("fminunc",
"fmincon").

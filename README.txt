# GLMspiketools

Fitting and simulation of Poisson generalized linear model for single and multi-neuron spike trains (Pillow et al 2008).


**Description:**  Simulates and computes maximum likelihood estimates for
the parameters of a Poisson GLM spike train model. Parameters
consist of a bank of stimulus filters ("receptive fields"),
spike-history filters, and coupling filters that capture dependencies
between neurons. The stimulus filter can be parametrized linearly or
bi-linearly, and the nonlinearity can be selected from a class
ensuring convexity of the negative log-likelihood, or parametrized using
using cubic splines. This model is a generalization of the
"Linear-Nonlinear-Poisson" model that incorporates spike-history
effects and correlations between neurons.

**Relevant publication:**
[Pillow et al, *Nature* 2008](http://pillowlab.princeton.edu/pubs/abs_ParkI_NN14.html)



Installation
==========

1. Download: clone the repository (```git clone git@github.com:pillowlab/GLMspiketools.git```) or
   [download as zip](https://github.com/pillowlab/GLMspiketools/archive/master.zip)
   and then unzip the archive.
2. Launch matlab and cd into the main directory containing the code
    (e.g. `cd  code/GLMspiketools/`).
3. CD to the sub-directory `tools_mexcode` and run the initialization
    script `initialize_mexcode` from the command line.  This will
    compile the C files in the sub-directory `mexcode` using the
    default compiler settings

Use
===

1. From the main code directory, run the `setpaths` script to add relevant
    sub-directories to the matlab path and initialize the global variable "RefreshRate".
2. Examine scripts in sub-directory `testscripts/` to see simple
    examples of simulation and fitting to spike data using the code:

**Test Scripts:**

1. `testscript_GLM.m` - fits plain GLM (temporal stimulus kernel only).
2. `testscript_GLM_coupled.m` - fits GLM with two coupled neurons.m
3. `testscript_GLM_spatialStim.m` - fits GLM using two different parametrizations
     of stimulus kernel (linear vs. bilinear) 
4. `testscript_GLM_splineNlin.m` - fits GLM incorporating an arbitrary
     nonlinearity parametrized by cubic splines.
     - uses coordinate ascent of filter and nonlinearity params (not concave!)
     - illustrates code for plotting spike rasters


Notes
=====

- time is represented in units of "stimulus frames", which is
controlled by the global variable RefreshRate (Hz).  Thus, for
example, if RefreshRate=100, each unit of time is 10ms, and a
post-spike kernel discretized in time bins of width dt=.02 has time
bins of length 0.2 ms.

- fitting code relies on the matlab optimization toolbox ("fminunc",
"fmincon").


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
[Pillow et al, *Nature* 2008](http://pillowlab.princeton.edu/pubs/abs_Pillow08_nature.html)



Download
==========

1. Either clone the repository from github (```git clone git@github.com:pillowlab/GLMspiketools.git```) or
   [download as zip](https://github.com/pillowlab/GLMspiketools/archive/master.zip)
   and then unzip the archive.


Basic Usage
===

1. From the main code directory (e.g., `~/Downloads/GLMspiketools/`), run the `setpaths` script to add relevant
    sub-directories to the matlab path.
2. Examine demo scripts in sub-directory `demos/` to see simple
    scripts illustrating how to simulate and fit the GLM to spike
    train data.


**Demo Scripts:**

1. `demo1_GLM_temporalStim.m` - simulates and fits GLM  with 1D (purely temporal) stimulus.
2. `demo2_GLM_spatialStim.m` - simulates and fits GLM  with 2D (space
   x time) stimulus, and illustrates both linear and bilinear
   parametrization of stimulus filter.
3. `demo3_GLM_coupled.m` - simulates and fits GLM with two coupled neurons


Notes
=====

- The code allows for two discretizations of time: `dtStim` specifies
the size of time bins representing a single frame of the stimulus, and
`dtSp` specifies the size of time bins for spikes (both in units of
seconds).  The code requires `dtSp` to evenly divide `dtStim`. Thus,
for example, if the stimulus has a refresh rate of 100 Hz and spikes
are represented with 1ms precision, then `dtStim=.01` and `dtSp=.001`.

- fitting code relies on the matlab optimization toolbox function "fminunc".

- An older release of this code (now sitting in branch `old_v1`) had
  functionality that is no longer supported.  Namely: cubic spline
  parametrization of the nonlinearity, and a smart "chunking" of the
  design matrix that was more memory efficient (albeit slightly
  slower). If memory issues are a problem, due to large stimulus or
  coupling from many neurons, we suggest checking out version v1. (In
  the shell use`git checkout old_v1`, or download directly:
  [GLMspiketools-old\_v1.zip](https://github.com/pillowlab/GLMspiketools/archive/old_v1.zip)).


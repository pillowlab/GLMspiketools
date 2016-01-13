% SETPATHS.m - GLMspiketools code repository
%
% This simple script sets the path to include relevant directories for the
% GLMspiketools code package.  You must 'cd' into this directory in order
% to evaluate it.
%
% More info: http://pillowlab.princeton.edu/code_GLM.html
% Github page: https://github.com/pillowlab/GLMspiketools


basedir = pwd;  % The directory where this script lives

% Add a bunch sub-directories (with absoluate path names)
addpath([basedir '/glmtools_fitting/']);
addpath([basedir '/glmtools_misc/']);
addpath([basedir '/glmtools_mex/']);
addpath([basedir '/nlfuns/']);


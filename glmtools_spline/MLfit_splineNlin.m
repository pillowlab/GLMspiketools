function [gg,negloglival,pp] = MLfit_splineNlin(gg,Stim)
% [gg,neglogli,pp] = MLfit_splineNlin(gg,Stim,optimargs)
% 
%  Fit a nonlinear function in an GLM neuron with a (positive) spline defined using a set of knots 
% (discontinuities of a piecewise-polynomial function), using maximum likelihood
% loss
%
% Inputs: 
%   gg = glm param structure
%   Stim = stimulus
%
% Outputs: 
%   fun - function pointer to nonlinearity (uses 'ppval')
%   pp - piecewise polynomial structure 
%   Mspline - matrix for converting function params to spline coeffs
%   splinePrs -  coefs=Mspline*splinePrse

% Compute net GLM filter output
[~,~,~,Itot] = neglogli_GLM(gg,Stim);
Itot = Itot-gg.dc; % remove effect of dc
gg.dc = 0; % set dc to zero

% Do fitting
[fspline,pp,negloglival] = fit_tspline_poisson(Itot,gg.sps,gg.splineprs,gg.dtSp);

% Insert params into glm struct
gg.ppstruct = pp;  % piecewise polynomial parameters
gg.nlfun = fspline; % insert nonlinearity


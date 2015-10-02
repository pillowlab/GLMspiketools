function [gg,fval,pp] = MLfit_splineNlin(gg,Stim,optimargs)
% [gg,fval,pp] = MLfit_splineNlin(gg,Stim,optimargs)
% 
%  Fit a nonlinear function in an GLM neuron with a (strictly
%  positive) spline defined using a set of knots 
% (discontinuities of a piecewise-polynomial function), using MSE
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

global SPNDS  % integer spike times

% Compute net GLM filter output
[neglogli,rr,tt,Itot] = neglogli_GLM(gg,Stim);
rr=[];tt=[];

% Sorted filter outputs and spike times
xxsp = sparse(length(Itot),1);
xxsp(SPNDS) = 1;
if gg.tspi(1) == 1; istrt = 1;
else; istrt = SPNDS(gg.tspi(1)-1)+1;
end
Itot = Itot(istrt:end);
xxsp = xxsp(istrt:end);
[Isrt,xi] = sort(Itot);
xxsp = xxsp(xi);
spnds = find(xxsp);

% Set up spline params
pp = gg.ppstruct;
splinePrs = gg.splinePrs;
[prs0,Mspline] = extractSplineParams(pp,splinePrs);
knots = pp.breaks;
minval = splinePrs.minval;

%fprintf(1, '--- ML fitting spline (MLfit_splineNlin.m) ---\n');

% Set optimization params
opts = optimset('maxiter', 20, 'maxfunevals', 1e6, ...
    'Gradobj', 'on', 'Hessian', 'on');
if nargin > 2
    opts = optimset(opts, optimargs{:});
end

% Fit by maximize likelihood
[splinePrs,fval]=fminunc(@Loss_splineGLM_neglogli,prs0,opts,Isrt,xxsp,pp,Mspline,minval,gg.dt);

% Create function handle for resulting spline
[funptr0,pp] = makeSplinePos(knots, Mspline*splinePrs,minval); 
nlfun = @(x)ppfunpos(pp,x,minval);
gg.nlfun = nlfun;
gg.ppstruct = pp;

% ========= OPTIONAL: error-checking stuff ============
% % Check Deriv & Hessian
%neglogli_SplineGLM(prs0,Isrt,xxsp,pp,Mspline,minval,gg.dt)
%HessCheck(@neglogli_SplineGLM,prs0,opts,Isrt,xxsp,pp,Mspline,minval,gg.dt);




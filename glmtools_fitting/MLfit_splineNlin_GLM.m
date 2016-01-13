function [gg,fval,pp] = MLfit_splineNlin_GLM(gg,Stim,optimargs)
% [gg,fval,pp] = MLfit_splineNlin_GLM(gg,Stim,optimargs)
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
[neglogli0,NONE,NONE,Itot] = neglogli_GLM(gg,Stim);

% Compute mask times
rlen = length(Itot);
[iiSpk,iiLi] = initfit_mask(gg.mask,gg.dt,rlen);  % bins to use for likelihood calc

% Make spike vector
xxsp = sparse(rlen,1);
xxsp(iiSpk) = 1;

% Extract relevant bins and sort by Itot
Itot = Itot(iiLi);
xxsp = xxsp(iiLi);
[Isrt,xi] = sort(Itot);
xxsp = xxsp(xi);
spnds = find(xxsp);

% Set up spline params
pp = gg.nlprs.ppstruct;
splinePrs = gg.nlprs.splinePrs;
[prs0,Mspline] = extractSplineParams(pp,splinePrs);
knots = pp.breaks;
minval = splinePrs.minval;

%fprintf(1, '--- ML fitting spline (MLfit_splineNlin_GLM.m) ---\n');

% Set optimization params
opts = optimset('maxiter',100,'maxfunevals',1e6, ...
      'Gradobj', 'on', 'Hessian', 'on','largescale', 'off', ...
      'tolfun',1e-12, 'TolX', 1e-12);
if nargin > 2
    opts = optimset(opts, optimargs{:});
end

% Fit by maximize likelihood
[splinePrs,fval]=fminunc(@Loss_splineNlin_GLMlogli,prs0,opts,Isrt,xxsp,pp,Mspline,minval,gg.dt);

% Create function handle for resulting spline
[funptr0,pp] = makeSplinePos(knots, Mspline*splinePrs,minval); 
nlfun = @(x)ppfunpos(pp,x,minval);
gg.nlfun = nlfun;
gg.nlprs.ppstruct = pp;

% ========= OPTIONAL: error-checking stuff ============
% % Check Deriv & Hessian
%neglogli_SplineGLM(prs0,Isrt,xxsp,pp,Mspline,minval,gg.dt)
%HessCheck(@neglogli_SplineGLM,prs0,opts,Isrt,xxsp,pp,Mspline,minval,gg.dt);




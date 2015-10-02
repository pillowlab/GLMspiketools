function [fun,pp,fval] = fitSplineGLM(gg,Stim,SpTimes,splinePrs,optargs)
% [fun,pp,fval] = fitSplineLNP(knots, x, y, smoothness,
% extrapDeg);
% 
%  Fit a nonlinear function in an LNP neuron with a (strictly
%  positive) spline defined using a set of knots 
% (discontinuities of a piecewise-polynomial function), using MSE
% loss
%
% Inputs: 
%   knots - breaks between polynomial pieces
%   xx - filtered stimulus (input to nonlinearity)
%   spnds - indices of xx where spikes occurred
%   smoothness - derivs of 1 less are continuous 
%      (e.g., smoothness=3 -> 2nd derivs are continuous)
%   extrapDeg - degree polynomial for extrapolation on ends
%   prs0 - initial params (coefs0 = Mspline*prs0;)
%
% Outputs: 
%   fun - function pointer to nonlinearity (uses 'ppval')
%   pp - piecewise polynomial structure 
%   Mspline - matrix for converting function params to spline coeffs
%   splinePrs -  coefs=Mspline*splinePrse

global SPNDS

% Compute net GLM currents
[Istm,Ihist,Icpl] = compGLMcurrents(gg,Stim,SpTimes,[0,size(Stim,1)]);
Itot = Istm+Ihist+Icpl;
Ihist=[];Istm=[];Icpl=[];% Free memory;

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
[prs0,Mspline] = extractSplineParams(pp,splinePrs);
knots = pp.breaks;
minval = splinePrs.minval;

fprintf(1, '--- ML fitting spline (fitSplineGLM.m) ---\n');

% Set optimization params
opts = optimset('maxiter', 250, 'maxfunevals', 1e6, ...
    'display','iter', 'Gradobj', 'on', 'Hessian', 'on');
if exist('optargs')
    opts = optimset(opts, optargs{:});
end

% % Check values and Derivs
%neglogli_SplineGLM(prs0,Isrt,xxsp,pp,Mspline,minval)
%HessCheck(@neglogli_SplineGLM,prs0,opts,Isrt,xxsp,pp,Mspline,minval);


% Optimize
[splinePrs,fval]=fminunc(@neglogli_SplineGLM,prs0,opts,Isrt,xxsp,pp,Mspline,minval,gg.dt);

% Create function handle for resulting spline
[fun,pp] = makeSplinePos(knots, Mspline*splinePrs,minval); 

% % ============== Old loss func  =====================
% function err = glm_splineLogLi(prs,xx,spnds,pp,M,minval);
% % negative log-likelihood loss for spline fitting

% global RefreshRate DTsim;

% pp.coefs = fliplr(reshape(M*prs,4,[])');
% fx = max(minval,ppval(pp,xx))/RefreshRate;

% err = -sum(log(fx(spnds))) + sum(fx)*DTsim;





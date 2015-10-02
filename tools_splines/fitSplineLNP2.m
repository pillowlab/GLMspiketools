function [fun,pp,Mspline,splinePrs,fval] = fitSplineLNP(knots, xx, spnds, smoothness, extrapDeg,prs0,minval);
% [fun,pp,Mspline,splinePrs] = fitSplineLNP(knots, x, y, smoothness, extrapDeg);
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

if nargin < 7
    minval = [];
end


Mspline = splineParamMatrix(knots,smoothness,extrapDeg);

fprintf(1, '... ML fitting spline (fitSpline.m) ...\n');

% Do optimization
opts = optimset('maxiter', 250, 'maxfunevals', 1e6, ...
                'display','iter', 'Largescale', 'off'); 
[splinePrs,fval] = fminunc(@splineLogLi,prs0,opts,knots,xx,spnds,Mspline,minval);

% Create function handle for resulting spline
[fun,pp] = makeSplinePos(knots, Mspline*splinePrs,minval); 

% ================ LOSS FUNCTION =====================
function err = splineLogLi(prs, knots,xx,spnds,M,minval);
% MSE loss for spline fitting

global RefreshRate DTsim;

ff = makeSplinePos(knots, M*prs, minval);
fx = ff(xx)/RefreshRate;

err = -sum(log(fx(spnds))) + sum(fx)*DTsim;





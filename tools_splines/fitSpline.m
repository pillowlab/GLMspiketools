function [fun,pp,Mspline,splinePrs]=fitSpline(knots,x,y,smoothness,extrapDeg);
% [fun,pp,Mspline,splinePrs]=fitSpline(knots,x,y,smoothness,extrapDeg);
% 
%  Fit a function y = f(x) with a spline defined using a set of knots
% (discontinuities of a piecewise-polynomial function), using MSE loss
%
% Inputs:
%   knots - breaks between polynomial pieces
%   x,y - spline minimizes (f(x)-y).^2
%   smoothness - derivs of 1 less are continuous
%           (e.g., smoothness=3 -> 2nd derivs are continuous)
%   extrapDeg - degree polynomial for extrapolation on ends
%
% Outputs:
%   fun - function pointer to nonlinearity (uses 'ppval')
%   pp - piecewise polynomial structure
%   Mspline - matrix for converting params to spline coeffs
%   splinePrs -  coefs=Mspline*splinePrse

Mspline = splineParamMatrix(knots,smoothness,extrapDeg);
prs0 = zeros(size(Mspline,2),1);  % Initial params

opts = optimset('maxiter', 1e4, 'maxfunevals',1e6,'display','on'); 
[splinePrs,err] = fminunc(@splineMSE,prs0,opts,knots,x(:),y(:),Mspline);

% Create function handle for resulting spline
[fun,pp] = makeSplineFun(knots, Mspline*splinePrs); 

% ================ LOSS FUNCTION =====================
function err = splineMSE(prs, knots,x,y,M);
% MSE loss for spline fitting

ff = makeSplineFun(knots, M*prs);
err = sum((y-ff(x)).^2);




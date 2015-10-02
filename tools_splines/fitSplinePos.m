function [fun,pp,Mspline,splinePrs]=fitSplinePos(knots,x,y,smoothness,extrapDeg,minval);
% [fun,pp,Mspline,splinePrs]=fitSplinePos(knots,x,y,smoothness,extrapDeg);
% 
%  Fit a function y = f(x) with a (strictly positive) spline
%  defined using a set of knots 
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

if nargin < 6
    minval = [];
end

Mspline = splineParamMatrix(knots,smoothness,extrapDeg);
prs0 = zeros(size(Mspline,2),1);  % Initial params

fprintf(1, '... Least-squares fitting spline (fitSplinePos.m) ...\n');

opts = optimset('maxiter', 1e4, 'maxfunevals', 1e6,'display','on', ...
                'Largescale', 'off'); 
[splinePrs,err] = fminunc(@splineMSE,prs0,opts,knots,x(:),y(:),Mspline,minval);

% Create function handle for resulting spline
[fun,pp] = makeSplinePos(knots, Mspline*splinePrs,minval); 

% ================ LOSS FUNCTION =====================
function err = splineMSE(prs, knots,x,y,M,minval);
% MSE loss for spline fitting

ff = makeSplinePos(knots,M*prs,minval);
err = sum((y-ff(x)).^2);




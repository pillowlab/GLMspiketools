function [fun,Mspline,splinePrs] = fitSplinePos_LogErr(knots, x, y, smoothness, extrapDeg);
% [fun,Mspline,splinePrs] = fitSplinePos(knots, x, y, smoothness, extrapDeg);
% 
%  Fit a function y = f(x) with a (strictly positive) spline defined using a set of knots
% (discontinuities of a piecewise-polynomial function), using MSE loss

Mspline = splineParamMatrix(knots,smoothness,extrapDeg);
prs0 = .01*ones(size(Mspline,2),1);  % Initial params

opts = optimset('maxiter', 1e4, 'maxfunevals', 1e6,'display','on', 'Largescale', 'off'); 
fprintf(1, '... Least-squares fitting spline (fitSpline.m) ...\n');
[splinePrs,err] = fminunc(@splineMSE,prs0,opts,knots,x(:),y(:),Mspline);

% Create function handle for resulting spline
fun = makeSplinePos(knots, Mspline*splinePrs); 

% ================ LOSS FUNCTION =====================
function err = splineMSE(prs, knots,x,y,M);
% MSE loss for spline fitting

ff = makeSplinePos(knots, M*prs);
y = max(y,1e-20);
err = sum((log(y)-log(ff(x))).^2);




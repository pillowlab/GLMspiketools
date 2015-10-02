function [fun,pp] = makeSplineFun(knots,x);
% [fun,pp] = makeSplineFun(knots,x);
% 
% Make (minimum-dimensional) param vector x into a cubic spline.  Returns
% function handle to piecewise nonlinear function.

coeffs = fliplr(reshape(x,4,[])');
pp = mkpp(knots, coeffs);
fun = @(x)ppval(pp,x);


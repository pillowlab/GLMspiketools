function [fun,pp] = makeSplinePos(knots,x,MinVal);
% [fun,pp] = makeSplinePos(knots,x);
% 
% Make (minimum-dimensional) param vector x into a cubic spline.  Returns
% function handle to piecewise polynomial function, which is also
% constrained to be strictly positive. 
%
%  MinVal = 1e-20;

if nargin<3
    MinVal = 1e-20;
end

coeffs = fliplr(reshape(x,4,[])');
pp = mkpp(knots, coeffs);
fun = @(x)max(MinVal,ppval(pp,x));


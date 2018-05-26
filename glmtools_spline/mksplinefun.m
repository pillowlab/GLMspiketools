function [fun,pp] = mksplinefun(breaks,x)
% [fun,pp] = mksplinefun(breaks,x)
% 
% Make spline function handle out of breaks (knots) and coefficients in vector x. 
%
% Inputs:
%   breaks - points of discontinuity
%   x - vector that can be reshaped into N x 4 matrix of spline coefficients
%       for each segment
%
% Outputs:
%   fun - function pointer for spline function
%    pp - matlab piecewise polynomial structure (evaluate with 'ppval').

coeffs = reshape(x,4,[])';
pp = mkpp(breaks, coeffs);
fun = @(x)ppval(pp,x);


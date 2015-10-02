function [prs,Mspline] = extractSplineParams(pp,splineprs)
% [prs,Mspline] = extractSplineParams(pp,splineprs)
% 
% Extract the minimal description of the parameters needed to make
% the spline embedded in pp (piece-wise polynomial struct), using
% spline parameters splineprs.

knots = pp.breaks;
coefs = pp.coefs;
Mspline = splineParamMatrix(knots,splineprs.smthness, ...
                            splineprs.extrapDeg);

fullprs = reshape(fliplr(coefs)',[],1);
prs = (Mspline'*Mspline)*Mspline'*fullprs;




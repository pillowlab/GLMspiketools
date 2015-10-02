function [ggnew,lival] = initializeSplineNlin_GLM(gg,Stim,prs,optimargs);
%  [ggnew,lival] = initializeSplineNlin_GLM(gg,Stim,spline_prs);
% 
%  Initializes parameters for a spline nonlinearity using likelihood under
%  a GLM.
%
%  Inputs: 
%     gg = param struct
%     Stim = stimulus
%     spline_prs = structure parameters for spline nonlinearity
%        ---- Fields: ----
%        nknots = # of knots,
%        epprob = [low,high] quantiles for first and last knots
%        smthness = 1+degree of smoothness (3 -> 2nd order continuous)
%        extrapDeg = [lo,hi] degree of polynomials for extrapolation
%        minval = min value of funtion (rectification) 
%      
%  Outputs:
%     ggnew = new param struct (with ML params);
%     fval = negative log-likelihood at ML estimate

[neglogli,rr,tt,II] = neglogli_GLM(gg,Stim);
rr=[];tt=[];

% Set knots and insert piece-wise linear nonlinearity.
knotwin  = quantile(II,prs.epprob);
knots = knotwin(1):diff(knotwin)/prs.nknots:knotwin(2);

% Set up spline params & initialize using minimal MSE
xx1 = min(II):.1:max(knots);
yy1 = exp(xx1)+10;
fprintf('... Least-squares fitting spline (fitSpline.m) ...\n');
[nlfun0,pp0,Mspline0,nlprs0]=fitSpline(knots,xx1,yy1,prs.smthness,prs.extrapDeg);
[nlfun0,pp0] = makeSplinePos(knots,Mspline0*nlprs0,prs.minval);

% Examine logli with this new nonlinearity
% gg.nlfun = nlfun0;
 
% Now fit nonlinearity via ML
gg.ppstruct = pp0;
gg.splinePrs = prs;
fprintf('... ML fitting spline (MLfit_splineNlin.m) ...\n');
[ggnew,fval]= MLfit_splineNlin(gg,Stim,optimargs);





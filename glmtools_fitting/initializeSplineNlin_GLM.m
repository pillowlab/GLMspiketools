function [ggnew,neglogli] = initializeSplineNlin_GLM(gg,Stim,prs,optimargs)
%  [ggnew,neglogli] = initializeSplineNlin_GLM(gg,Stim,prs,optimargs)
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

[NONE,NONE,NONE,II] = neglogli_GLM(gg,Stim); % Get total currents from GLM

% Set knots and insert piece-wise linear nonlinearity.
knotwin  = quantile(II,prs.epprob);
knots = knotwin(1):diff(knotwin)/prs.nknots:knotwin(2);

% Initialize spline by minimizing MSE
xx1 = min(II):.1:max(knots);
yy1 = exp(xx1)+10;
fprintf('... Least-squares fitting spline (fitSpline.m) ...\n');
[NONE,NONE,Mspline0,nlprs0]=fitSpline(knots,xx1,yy1,prs.smthness,prs.extrapDeg);

% Make a strictly positive verstion of this spline
[zz,pp0] = makeSplinePos(knots,Mspline0*nlprs0,prs.minval);

% Now fit nonlinearity via maximum likelihood
gg.nlprs.ppstruct = pp0;
gg.nlprs.splinePrs = prs;
fprintf('... ML fitting spline (MLfit_splineNlin.m) ...\n');
[ggnew,neglogli]= MLfit_splineNlin_GLM(gg,Stim,optimargs);





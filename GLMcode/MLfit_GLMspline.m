function [gg,fval,H] = MLfit_GLMspline(gg,Stim,optimargs,npasses)
% [gg,fval,H] = MLfit_GLMspline(gg,Stim,optimargs,npasses)
% 
%  Fit GLM filters and spline nonlinearity using alternating 
%  ascent of filter and spline parameters
%
% Inputs: 
%   gg = glm param structure
%   Stim = stimulus
%   optimargs = optimization params
%   npasses (OPTIONAL) = # of passes of filter and nonlinear param fitting
%
% Outputs: 
%   gg = GLM struct;
%   fval = negative log-likelihood at maximum
%   H = Hessian (for filter params) at max

if nargin <= 3
    npasses = 3;
end

if isfield(gg, 'kxbas');
    fitFun = @MLfit_GLMbi;
else
    fitFun = @MLfit_GLM;
end

livalPrev = 0;
lival = neglogli_GLM(gg,Stim);
for j = 1:npasses
    % Fit filters
    fprintf('\n--- Pass %d: fitting filters (dLoss=%.3f)----', j,livalPrev-lival);
    [gg,lival] = fitFun(gg,Stim,optimargs);
    livalPrev = lival;
    
    % Fit spline nonlinearity
    fprintf('\n--- Pass %d: Fitting spline nlin -------', j);
    [gg,lival] = MLfit_splineNlin(gg,Stim,optimargs);
end    

fprintf('\n--- Final step: fitting filters (dLoss=%.3f)----',livalPrev-lival);
[gg,fval,H] = fitFun(gg,Stim,optimargs);



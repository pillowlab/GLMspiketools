function [gg,neglogli,H] = MLfit_GLMwithspline(gg,Stim,optimargs,npasses)
% [gg,neglogli,H] = MLfit_GLMwithspline(gg,Stim,optimargs,npasses)
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

% Select which GLM fitting function to use
if isfield(gg, 'kxbas')
    fitFun = @MLfit_GLMbi;
else
    fitFun = @MLfit_GLM;
end

neglogliPrev = 0;
neglogli = neglogli_GLM(gg,Stim);
for j = 1:npasses
    % Fit filters
    fprintf('\n--- Pass %d: fitting filters (dLoss=%.3f)----\n', j,neglogliPrev-neglogli);
    [gg,neglogli] = fitFun(gg,Stim,optimargs);
    neglogliPrev = neglogli;
    
    % Fit spline nonlinearity
    fprintf('--- Pass %d: Fitting spline nlin -------\n', j);
    [gg,neglogli] = MLfit_splineNlin(gg,Stim);
end    

% Filters one last time and return Hessian for filters (if desired)
fprintf('--- Final step: fitting filters (dLoss=%.3f)----',neglogliPrev-neglogli);
[gg,neglogli,H] = fitFun(gg,Stim,optimargs);

% %----------------------------------------------------
% Debugging code
% %----------------------------------------------------
%
% % ------ Check analytic gradients and Hessians -------
%  HessCheck(fitFun,prs0,opts);
%  HessCheck_Elts(@Loss_GLM_logli, [1 12],prs0,opts);
%  tic; [lival,J,H]=lfunc(prs0); toc;


function [gg,fval,H,Xstruct] = MLfit_GLMbi(gg,Stim,optimArgs)
%  [gg,fval,H] = MLfit_GLMbi(gg,Stim,optimArgs);
% 
%  Computes the ML estimate for GLM params, using grad and hessians.
%  Assumes bilinear parametrization of space-time filter.
%
%  Inputs: 
%     gg = param struct
%     Stim = stimulus
%     optimArgs = cell array of optimization params (optional)
%
%  Outputs:
%     ggnew = new param struct (with ML params);
%     fval = negative log-likelihood at ML estimate
%        H = Hessian of negative log-likelihood at ML estimate
%  Xstruct = structure with design matrices for spike-hist and stim terms

% Set optimization parameters 
if nargin > 2
    opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
else
    opts = optimset('Gradobj','on','Hessian','on','display','iter');
end

% --- Create design matrix extract initial params from gg ----------------
[prs0,Xstruct] = setupfitting_GLM(gg,Stim);

% minimize negative log likelihood
lfunc = @(prs)Loss_GLMbi_logli(prs,Xstruct); % loss function for exponential nonlinearity
[prs,fval] = fminunc(lfunc,prs0,opts); % optimize negative log-likelihood for prs

% Compute Hessian at maximum, if requested     
if nargout > 2
    [fval,~,H] = Loss_GLMbi_logli(prs);
end

% Put returned vals back into param structure ------
gg = reinsertFitPrs_GLMbi(gg,prs,Xstruct);


% %----------------------------------------------------
% Optional debugging code
% %----------------------------------------------------
%
% % ------ Check analytic gradients, Hessians -------
% DerivCheck(lfunc,prs0,opts);
% HessCheck_Elts(lfunc, [1 12],prs0,opts);
% tic; [lival,J,H]=lfunc(prs0,Xstruct); toc;


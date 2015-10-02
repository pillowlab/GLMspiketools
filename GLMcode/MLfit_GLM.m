function [gg, fval,H] = MLfit_GLM(gg,Stim,optimArgs);
%  [ggnew,fval,H] = MLfit_GLM(gg,Stim,optimArgs);
% 
%  Computes the ML estimate for GLM params, using grad and hessians.
%  Assumes basis for temporal dimensions of stim filter
%
%  Inputs: 
%     gg = param struct
%     Stim = stimulus
%     optimArgs = cell array of optimization params (optional)
%
%  Outputs:
%     ggnew = new param struct (with ML params);
%     fval = negative log-likelihood at ML estimate

MAXSIZE  = 1e7;  % Maximum amount to be held in memory at once;

% Set optimization parameters 
if nargin > 2
    opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
else
    opts = optimset('Gradobj','on','Hessian','on','display','iter');
end

% Set initial params
prs0 = extractFitPrs_GLM(gg,Stim,MAXSIZE);

% minimize negative log likelihood
[prs,fval] = fminunc(@Loss_GLM_logli,prs0,opts);
if nargout > 2 % Compute Hessian if desired
    [fv,gradval,H] = Loss_GLM_logli(prs);
end

% Put returned vals back into param structure ------
gg = reinsertFitPrs_GLM(gg,prs);

%----------------------------------------------------
% % ------ Check analytic gradients, Hessians -------
% HessCheck(@Loss_GLM_logli,prs0,opts);
% HessCheck_Elts(@Loss_GLM_logli, [1 12],prs0,opts);
% tic; [lival,J,H]=Loss_GLM_logli(prs0); toc;


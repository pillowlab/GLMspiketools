function [gg,fval,H] = MLfit_GLM(gg,Stim,optimArgs)
%  [ggnew,fval,H] = MLfit_GLM(gg,Stim,optimArgs)
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
%        H = Hessian of negative log-likelihood at ML estimate

% Set optimization parameters 
if nargin > 2
    opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
else
    opts = optimset('Gradobj','on','Hessian','on','display','iter');
end

% --- Create design matrix using bases and extract initial params from gg -------
[prs0,Xstruct] = setupfitting_GLM(gg,Stim);

Loss_GLM_logli(prs0,Xstruct)

% minimize negative log likelihood --------------------
lfunc = @(prs)Loss_GLM_logli(prs,Xstruct); % loss function
[prsML,fval] = fminunc(lfunc,prs0,opts); % find ML estimate of params

% Compute Hessian if desired
if nargout > 2 
    [fval,~,H] = Loss_GLM_logli(prsML,Xstruct);
end

% Put returned vals back into param structure ------
gg = reinsertFitPrs_GLM(gg,prsML,Xstruct);

% %----------------------------------------------------
% Optional debugging code
% %----------------------------------------------------
%
% % ------ Check analytic gradients and Hessians -------
%  HessCheck(lfunc,prs0,opts);
%  HessCheck_Elts(@Loss_GLM_logli, [1 12],prs0,opts);
%  tic; [lival,J,H]=Loss_GLM_logli(prs0); toc;


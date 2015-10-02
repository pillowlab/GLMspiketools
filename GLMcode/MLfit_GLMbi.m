function [gg, fval,H] = MLfit_GLMbi(gg,Stim,optimArgs);
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
%
%  Created: April 08.  (from maxli_GLM)


% Set optimization parameters 
if nargin > 2
    opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
else
    opts = optimset('Gradobj','on','Hessian','on','display','iter');
end

% Set initial params
MAXSIZE  = 1e7;  % Maximum amount to be held in memory at once;
prs0 = extractFitPrs_GLMbi(gg,Stim,MAXSIZE);

% minimize negative log likelihood
[prs,fval] = fminunc(@Loss_GLMbi_logli,prs0,opts);
if nargout > 2 % Compute Hessian if desired
    [fv,gradval,H] = Loss_GLMbi_logli(prs);
end

% Put returned vals back into param structure ------
gg = reinsertFitPrs_GLMbi(gg,prs);

%----------------------------------------------------
% % ------ Check analytic gradients, Hessians -------
% HessCheck(@Loss_GLMbi_logli,prs0,opts);
% HessCheck_Elts(@Loss_GLMbi_logli, [1 12],prs0,opts);
% tic; [lival,J,H]=Loss_GLMbi_logli(prs0); toc;


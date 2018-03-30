function [gg,neglogli,H,Xstruct,neglogp] = MAPfit_GLM(gg,Stim,Cinv,optimArgs)
%  [gg,neglogli,H,Xstruct,neglogp] = MAPfit_GLM(gg,C,Stim,optimArgs)
% 
%  Computes the MAP estimate for GLM params, using grad and hessians under
%  a zero-mean Gaussian prior with inverse covariance Cinv.
%
%  Minimizes negative log-likelihood plus a penalty of the form
%     0.5*x'*Cinv*x, where x is the parameter vector
%
%  Assumes basis for temporal dimensions of stim filter
%
%  Inputs: 
%  -------
%         gg = param struct
%       Stim = stimulus
%       Cinv = inverse of prior covariance matrix (gives quadratic penalty)
%  optimArgs = cell array of optimization params (optional)
%
%  Outputs:
%  -------
%     ggnew = new param struct (with MAP params);
%  neglogli = negative log-likelihood at MAP estimate
%         H = Hessian of negative log-likelihood at MAP estimate
%   Xstruct = structure with design matrices for spike-hist and stim terms
%   neglogp = negative log-posterior at MAP estimate

% Set optimization parameters 
algopts = getFminOptsForVersion(version);
if nargin > 3, opts = optimset(algopts{:}, optimArgs{:});
else, opts = optimset(algopts{:});
end

% --- Create design matrix extract initial params from gg ----------------
[prs0,Xstruct] = setupfitting_GLM(gg,Stim);

% --- Set log-likelihood function ----------------------------
if isequal(Xstruct.nlfun,@expfun) || isequal(Xstruct.nlfun,@exp)
    % loss function for exponential nonlinearity
    lfunc = @(prs)Loss_GLM_logli_exp(prs,Xstruct);
else
    lfunc = @(prs)Loss_GLM_logli(prs,Xstruct); 
    % loss function for all other nonlinearities
end

% set log-posterior function
lpost = @(prs)(neglogpost(prs,lfunc,Cinv));

% --- minimize negative log posterior --------------------
[prsMAP,neglogp] = fminunc(lpost,prs0,opts); % find MAP estimate of params

% Compute Hessian of log-likelihood, if desired
if nargout > 1 
    [neglogli,~,H] = lfunc(prsMAP);
end

% Put returned vals back into param structure
gg = reinsertFitPrs_GLM(gg,prsMAP,Xstruct);

end

%% =================================================
function [negLP,dLP,H] = neglogpost(prs,lfunc,Cinv)
% Compute log-posterior by adding quadratic penalty to log-likelihood

nprs = size(Cinv,1); % number of parameters in C
preg = prs(1:nprs);  % parameters being regularized.

switch nargout

    case 1  % evaluate function
        negLP = lfunc(prs) + .5*preg'*Cinv*preg;
    
    case 2  % evaluate function and gradient
        [negLP,dLP] = lfunc(prs);
        negLP = negLP + .5*preg'*Cinv*preg;        
        dLP(1:nprs) = dLP(1:nprs) + Cinv*preg;

    case 3  % evaluate function and gradient
        [negLP,dLP,H] = lfunc(prs);
        negLP = negLP + .5*preg'*Cinv*preg;        
        dLP(1:nprs) = dLP(1:nprs) + Cinv*preg;
        H(1:nprs,1:nprs) = H(1:nprs,1:nprs) + Cinv;
end

end


% %----------------------------------------------------
% Optional debugging code
% %----------------------------------------------------
%
% % ------ Check analytic gradients and Hessians -------
%  HessCheck(lfunc,prs0,opts);
%  HessCheck_Elts(@Loss_GLM_logli, [1 12],prs0,opts);
%  tic; [lival,J,H]=lfunc(prs0); toc;


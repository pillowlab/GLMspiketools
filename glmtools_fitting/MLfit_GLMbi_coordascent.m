function [gg,fval,H] = MLfit_GLMbi_coordascent(gg,Stim,maxiter,ftol,optimArgs)
%  [gg,fval,H] = MLfit_GLMbi_coordascent(gg,Stim,maxiter,ftol,optimArgs);
% 
%  Computes the ML estimate for GLM params, with gradient and hessians, via
%  coordinate ascent, for bilinear (low rank) parametrization of space-time filter.
%
%  Inputs: 
%     gg = param struct
%     Stim = stimulus
%     maxiter = maximum number of coordinate ascent steps 
%     ftol = tolerance for stopping based on function improvement 
%     optimArgs = cell array of optimization params (optional)
%
%  Outputs:
%     ggnew = new param struct (with ML params);
%     fval = negative log-likelihood at ML estimate
%        H = Hessian of negative log-likelihood at ML estimate

% Set maxiter if necessary
if (nargin <= 2) || isempty(maxiter)
    maxiter = 50;  % maximum number of coordinate ascent steps
end
if (nargin <=3) || isempty(ftol)
    ftol = .001;  % tolerance for stopping
end

% Set optimization parameters 
if nargin < 5
    optimArgs = [];
end

% create struct for spatial optimization
ggx = gg;
ggx.ktype = 'linear'; % set type to linear (from bilinear)
ggx.ktbas = 1; ggx.ktbasprs = []; % remove temporal basis
ggx.kt = gg.kx(:)';
ggx.k = ggx.kt;

% create struct for temporal optimization
ggt = gg;
ggt.ktype = 'linear'; % set type to linear (from bilinear)
ggt.kt = gg.kt;
ggt.k = ggt.ktbas*ggt.kt;

% Initialize spatial stimulus 
[nt,nkx] = size(Stim);
krank = gg.krank;
xStim = zeros(nt,krank*nkx); % initialize spatial stimulus

% compute initial log-likelihood
neglogli0 = neglogli_GLM(gg,Stim); % Compute logli of initial params
dlogli = inf;  % initialize change in logli
jjiter = 0;  % initialize counter

while (jjiter<maxiter) && dlogli>ftol
    
    % ---- Update temporal params -----
    fprintf('Iter #%d: Updating temporal params\n', jjiter);
    tStim = Stim*reshape(ggx.k',[],gg.krank);
    ggt.dc = ggx.dc;  % update dc param
    [ggt,tneglogli] = MLfit_GLM(ggt,tStim,optimArgs);
    fprintf('  dlogli = %.4f\n', neglogli0-tneglogli);
    
    % Convolve stimulus with temporal filters
    for irank = 1:krank
        for icol = 1:nkx
            xStim(:,icol+(irank-1)*nkx) = sameconv(Stim(:,icol),ggt.k(:,irank));
        end
    end
    
    % ---- Update spatial params ----
    fprintf('Iter #%d: Updating spatial params\n', jjiter);
    ggx.dc = ggt.dc; % update dc param
    [ggx,xneglogli] = MLfit_GLM(ggx,xStim,optimArgs);
    fprintf('  dlogli = %.4f\n', neglogli0-xneglogli);
    
    % Update iters
    jjiter = jjiter+1;  % counter
    dlogli = neglogli0-xneglogli; % change in log-likelihood
    neglogli0 = xneglogli;
    
end

fprintf('\nFinished coordinate ascent: %d iterations (dlogli=%.6f)\n',jjiter,dlogli);

% Compute conditional Hessians, if desired
if nargout > 2

    % Compute Hessian for time components
    tStim = Stim*(ggx.k');
    ggt.dc = ggx.dc;  % update dc param
    [ggt,~,Ht] = MLfit_GLM(ggt,tStim,optimArgs);

    % Compute Hessian for space components
    for irank = 1:krank
        for icol = 1:nkx
            xStim(:,icol+(irank-1)*nkx) = sameconv(Stim(:,icol),ggt.k(:,irank));
        end
    end    
    ggx.dc = ggt.dc; % update dc param
    [ggx,neglogli0,Hx] = MLfit_GLM(ggx,xStim,optimArgs);
    
    H = {Ht,Hx};  % Hessians for t and x components
end

% Update params of bilinear model
gg.dc = ggx.dc;
gg.kt = ggt.kt;
gg.kx = reshape(ggx.k',[],gg.krank);
gg.k = (gg.ktbas*gg.kt)*gg.kx';
fval = neglogli0; % value of loglikelihood


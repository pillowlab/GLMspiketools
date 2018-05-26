function [prs,fval,H] = fminNewton(fptr,prs,opts)
% [prs,fval,H] = fminNewton(fptr,prs,opts)
%
% Simple mplementation of Newton's method for minimizing a function
%
% Inputs:
%   fptr - function handle for loss function
%   prs -  initial value of parameters
%   opts - struct with fields: 'tolX', 'tolFun', 'maxIter', 'verbose'.
%          (verbose takes on values 0, 1 or 2)
%
% Outputs: 
%   prs - minimizer of function (or best guess)
%   fval - value of function at minimizer
%   H    - Hessian at minimizer
% 
% Derivation of Newton's method (refresher): 
% ------------------------------
% f = f0 + f1(x-x0) + 1/2 (x-x0)^T H (x-x0) %% 2nd order Taylor expansion 
% df = f1 + H(x-x0)                         %% derivative
% x = x0 - H^{-1} f1                        %% solve for x (new point)
% ------------------------------
% 
% $Id$

if nargin < 3
    opts = [];
end

% set options to defaults, if necessary
if isfield(opts,'tolX'),tolX= opts.tolX;          else tolX    = 1e-8; end
if isfield(opts,'tolFun'),tolFun= opts.tolFun;    else tolFun  = 1e-8; end
if isfield(opts,'maxIter'),maxIter= opts.maxIter; else maxIter = 50; end
if isfield(opts,'verbose'),verbose= opts.verbose; else verbose = 0; end

% Initialize 
dX = inf;
dF = inf;
iter = 0;
canDescend = true;
nprs = length(prs);

% % Change warning states so ill-conditioned matrix inverses aren't reported
% msgid1 = 'MATLAB:illConditionedMatrix';
% msgid2 = 'MATLAB:singularMatrix';
% msgid3 = 'MATLAB:nearlySingularMatrix';
% wstate = warning('off',msgid1);
% warning('off',msgid2);
% warning('off',msgid3);
% lastwarn('');

while (dX>tolX) && (dF>tolFun) && (iter < maxIter) && canDescend
    
    % evaluate function, gradient & Hessian 
    [f,df,H] = fptr(prs);  

    if any(isinf(H)) | any(isnan(H))
	error('fminNewton2: Infinity or Nan in user-supplied Hessian');
    end
    
    % Test for positive definite-ness
    try chol(H); REG=0;
    catch fprintf('[[not pos-def]]');REG=1;
    end
    
    % Check to see if H is ill-conditioned
    if ~REG && (cond(H)<1e6)
	% ---- Standard Newton's Method if well-conditioned ----
	dprs = -H\df;

% 	% Check if it's a descent direction
% 	if dprs'*df > 0 
% 	    dprs = -dprs; % flip sign otherwise
% 	    fprintf(1, 'FLIPPED\n');
% 	end
    else
	% ---- Regularize if ill-conditioned ----
	if verbose>=2
	    fprintf('(regularizing H)');
	end
	% Compute eigenvalues of H
	eigvals = eig(H); eigmin = min(eigvals);eigmax = max(eigvals); 
	if eigmin > 0
	    lam = eigmax/1e3;
	else
	    %lam = max(eigmax/1e6,-2*eigmin);
	    lam = eigmax/1e3-eigmin;
	end
	dprs = -(H+lam*speye(nprs))\(df);  % next attempt
	
	% Make sure *this* is a descent direction (if not, something is wrong)
	if dprs'*df > 0
	    error('Somehow we didn''t end up with a descent direction');
	end

    end
    
    % Evaluate function at new value
    fnew = fptr(prs+dprs);
    if fnew < f
	% Accept the position: update bookkeeping params
	dX = norm(dprs);
	dF = -(fnew-f);
	% Update params
	prs = prs+dprs;

    else
	% Line search in this search direction
	% (Note: will give errs if loss func doesn't accept matrix inputs)
	if verbose >= 2
	    fprintf('(line search 1)');
	end
	fracs = 3.^(-(1:5)); % fractions of gradient to move
	fnewvec = fptr(bsxfun(@plus,prs,+dprs*fracs));
	[fnew,imin] = nanmin(fnewvec);
	
	if ~(fnew<f)
	    % try again with even smaller steps
	    if verbose >= 2
		fprintf('(line search 2)');
	    end
	    fracs = 3.^(-(6:15)); % fractions of gradient to move
	    fnewvec = fptr(bsxfun(@plus,prs,+dprs*fracs));
	    [fnew,imin] = nanmin(fnewvec);
	end
	
	if isinf(fnew) || isnan(fnew)
	    error('fminNewton2: unavoidable infinity or Nan function value: try improving line search');
	end

	if ~(fnew<f)
	    canDescend = 0;
	    if verbose>=1
		fprintf('fminNewton: can''t descend anymore; minimum possible\n');
	    end
	else
	    % Accept the position: update bookkeeping params
	    dprs = dprs*fracs(imin);
	    dX = norm(dprs);
	    dF = -(fnew-f);
	    % Update params
	    prs = prs+dprs;

	end
	
    end
    
    % Print out information if requested
    if verbose>=2
	fprintf('iter %d: fun=%.4f\n', iter,fnew);
    end
    iter = iter+1;
    
    % % Check to make sure we're only decreasing function
    % if canDescend && (dF<0)
    %	warning('error in fminNewton: function increased somehow!'); % for debugging
    %	keyboard;
    % end
end

% If desired, Compute value and Hessian at final value
if nargout > 1
    [fval,~,H] = fptr(prs);
end

% Report stopping condition, if verbose >= 1
if verbose >=1
    if ~canDescend && iter>3
	fprintf('fminNewton could not descend further (iter=%d): flat local minimum likely\n',iter);
    elseif ~canDescend
	fprintf('fminNewton could not descend (iter=%d): quadratic approx may be bad\n',iter);
	fprintf('Try implementing a linesearch method or using better initialization\n');
    elseif dX<tolX
	fprintf('fminNewton finished (iter=%d): change in X less than tolX\n',iter);
    elseif dF<tolFun
	fprintf('fminNewton finished (iter=%d): change in fun less than tolFun\n',iter);
    else
	fprintf('fminNewton stopped: iter=%d exceeded maxIter',iter);
    end
end

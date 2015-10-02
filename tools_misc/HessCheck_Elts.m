function dH = HessCheck_Elts(funptr, inds, X0, opts, varargin);
% HessCheck_Elts(funptr, inds, X0, opts, varargin);
% Checks specific entries in the Hessian (2nd deriv) of a function
% using finite differencing

args = {funptr, X0, varargin{:}};
[f, JJ, HH] = feval(funptr, X0, varargin{:});

eps = .00001;
nprs = length(X0);
mask = zeros(nprs,2);
mask(inds(1),1) = eps;  
mask(inds(2),2) = eps;

vv(1,1) = feval(funptr, X0+mask*[-1;-1], varargin{:});
vv(2,1) = feval(funptr, X0+mask*[1;-1], varargin{:});
vv(1,2) = feval(funptr, X0+mask*[-1;1], varargin{:});
vv(2,2) = feval(funptr, X0+mask*[1;1], varargin{:});

val = diff(diff(vv))/eps.^2/4;

fprintf(1, 'Analytic Hess: %f, Finite Diffs: %f\n', HH(inds(1),inds(2)),val)

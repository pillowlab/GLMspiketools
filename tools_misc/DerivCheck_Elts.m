function dH = DerivCheck_Elts(funptr, ind, X0, opts, varargin);
% DerivCheck_Elts(funptr, inds, X0, opts, varargin);
% Checks specific entries in the Hessian (2nd deriv) of a function
% using finite differencing

args = {funptr, X0, varargin{:}};
[f, JJ] = feval(funptr, X0, varargin{:});

eps = 1e-5;
nprs = length(X0);
mask = zeros(nprs,1);
mask(ind) = eps;  

vv(1) = feval(funptr, X0-mask, varargin{:});
vv(2) = feval(funptr, X0+mask, varargin{:});

val = diff(vv)/eps./2;

fprintf(1, 'Analytic Deriv: %f, Finite Diffs: %f\n', JJ(ind),val)

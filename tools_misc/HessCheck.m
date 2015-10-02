function HessCheck(funptr, X0, opts, varargin);
% HessCheck(funptr, X0, opts, varargin);
% Checks the Hessian (2nd deriv) of a function 'funptr' for optimization purposes
%  Call the same as you would call the optimization routine.
args = {funptr, X0, varargin{:}};
[f, JJ, HH] = feval(args{:});

tol = 1e-6;  % Size of finite differencing step
rr = randn(length(X0),1).*tol;

X2 = X0+rr;
args2 = {funptr, X2, varargin{:}};
[f2, JJ2, HH2] = feval(args2{:});

normDerivs = [norm(JJ)-norm(JJ2)];  % If this is close to zero, then analytic deriv may not be useful;

% Print results
fprintf('Derivs: Analytic vs. Finite Diff = [%.4e, %.4e]\n', dot(rr, JJ), f2-f);
fprintf('Hess: Analytic vs. Finite Diff = [%.4e, %.4e]\n', sum(HH*rr), sum(JJ2-JJ));

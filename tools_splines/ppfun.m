function [f,df,ddf]=ppfun(pp,xx)
%  [f,df,ddf]=ppfun(pp,xx)
% 
%   Evaluate piecewise polynomial 'pp' at values 'xx'.  Optionally returns first and 
%   second derivs at these points.
%   
%   Inputs: 
%      pp = piecewise polynomial structure, given (eg) by makepp or spline
%      xx = vector of values at which to evaluate

if isstruct(xx) % flip input args if ppval(xx,pp) was called
   temp = xx; xx = pp; pp = temp;
end

% make into row vector
sz = size(xx);
lx = numel(xx); 
xs = reshape(xx,1,lx);

% if necessary, sort xs
if any(diff(xs)<0), [xs,ix] = sort(xs); end

% take apart PP
[breaks,coefs,l,k,dd]=unmkpp(pp);

% for each data point, compute its breakpoint interval
[ignored,index] = sort([breaks(1:l) xs]);
index = max([find(index>l)-(1:lx);ones(1,lx)]);

% now go to local coordinates
xs = xs-breaks(index);

% Recursively compute polynomial and its derivs
f = coefs(index,1);
for i=2:k
    f = xs(:).*f + coefs(index,i); 
end
if nargout > 1
    df = (k-1)*coefs(index,1);
    for i = 2:k-1
       df = xs(:).*df + (k-i)*coefs(index,i);
    end
end
if nargout > 2
    ddf = (k-2)*(k-1)*coefs(index,1);
    for i = 2:k-2
       ddf = xs(:).*ddf + (k-i)*(k-i-1)*coefs(index,i);
    end
end

if exist('ix'), 
    f(ix) = f; 
    if (nargout > 1), df(ix) = df;   end
    if (nargout > 2), ddf(ix) = ddf; end
end
f = reshape(f,sz);
if nargout > 1, df = reshape(df,sz); end
if nargout > 2, ddf = reshape(ddf,sz); end

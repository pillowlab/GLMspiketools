function [L,dL,ddL]=Loss_splineGLM_neglogli(prs,xx,yy,pp,M,minval,DTsim);
%  [L,dL,ddL]=neglogli_SplineGLM(pp,xx,spnds,pp,M,minval);
% 
%   Evaluates negative likelihood of GLM model with spline nonlinearity.
%   

global RefreshRate

% Extract spline params
[breaks,coefs0,nsegs,k,dd]=unmkpp(pp);
coefs = fliplr(reshape(M*prs,4,[])');
lx = numel(xx); 

% if necessary, sort xx
if any(diff(xx)<0)
    error('inputs should be sorted!');
end

% for each data point, compute its breakpoint interval
[ignored,index] = sort([breaks(1:nsegs) xx']);
index = max([find(index>nsegs)-(1:lx);ones(1,lx)])';

% now go to local coordinates
xx = xx-breaks(index)';

% Recursively compute polynomial and its derivs
f = coefs(index,1).*xx.^3+coefs(index,2).*xx.^2+ ...
    coefs(index,3).*xx + coefs(index,4);

% -----  Check for negative values -------
ii = find(f<minval);
if ~isempty(ii)
    f(ii) = minval;
end

% Compute Likelihood
L = -sum(log(f(find(yy)))) + sum(f)*DTsim/RefreshRate;

% Compute gradients
if nargout > 1
    fnonzero = (f>minval);
    trm0 = zeros(nsegs,4);
    trm1 = zeros(nsegs,4);
    H = zeros(nsegs*4,nsegs*4);
    for j = 1:nsegs
        ii = find((index == j) & fnonzero);
        if ~isempty(ii)
            xs = xx(ii); 
            nxs = length(xs);
            trm0(j,:) = [sum([xs.^3,xs.^2,xs],1),nxs];
            iisp = find(yy(ii));
            xxsp = xs(iisp);
            nsp = length(iisp);
            if nsp > 0
                Mf = spdiags(1./f(ii(iisp)),0,nsp,nsp);
                fprime = [xxsp.^3, xxsp.^2,xxsp,ones(nsp,1)];
                trm1(j,:) = sum(Mf*fprime,1);
                H((j-1)*4+1:j*4,(j-1)*4+1:j*4) = ...
                    fliplr(fprime)'*(Mf.^2)* ...
                    fliplr(fprime);
            end
        end
    end
    dL = trm0*DTsim/RefreshRate - trm1;
    dL = M'*reshape(fliplr(dL)',[],1);
    ddL = M'*H*M;
end


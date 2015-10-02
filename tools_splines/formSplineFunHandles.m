function [ff,splfuns] = formSplineFunHandles(prs,splineStruct);
% [ff,splfuns] = formSplineFunHandles(paramvec,splineStruct);
% 
% Re-insert reduced parameters into a piecewise polynomial struct (e.g.,
% after fitting a sum of splines to data).
%
% Inputs: Y - dependent variable (column vector)
%         X - indep variables (each column is a diff regressor)
%         splineStruct - structure with fields: "breaks", "smoothness", "extrapDeg"
%            - Use cell array if different params for each regressor
% 
% Outputs: ff - handle to function ff(x) = [f1(x(:,1), f2(x(:,2),  ... ]
%          splfuns - cell array of individual functions f1, f2, ....

nsplines = length(splineStruct);

% 1. Extract params for each spline structure 
for ispl = 1:nsplines
   
    % Get spline structure for this column (regressor)
    if iscell(splineStruct)
        spstr{ispl} = splineStruct{ispl};
    else
        spstr{ispl} = splineStruct;
    end
    ss = spstr{ispl};
    
    % Compute spline parametrization matrix
    Mspline{ispl} = splineParamMatrix(ss.breaks, ss.smoothness, ss.extrapDeg);
    nsegs(ispl) = length(ss.breaks)-1;  % number of segments
    breaks{ispl} = ss.breaks;
    nprs(ispl) = size(Mspline{ispl},2);

end

% 2. Compute piecewise polynomial coeffs and create func handle
for ispl = 1:nsplines
    iprs = sum(nprs(1:ispl-1));
    splprs  = prs(iprs+[1:nprs(ispl)]);
    coeffs = fliplr(reshape(Mspline{ispl}*splprs,4,[])');
    pp = mkpp(breaks{ispl}, coeffs); % make piecewise polynomial
    splfuns{ispl} = @(x)ppval(pp,x);
end

ff = @(x)splineVec(splfuns,x);


% ------------------------------
function y = splineVec(f,x)
% Evaluates a sum of functions 
%
%  Inputs:
%     f - cell array of function handles
%     x - column i is the input to the corresponding function f{i}
%
%  Output:
%     y = f{1}(x(:,1) + f{2}(x(:,2)) + .... + f{n}(x(:,n))
%
  nfuns = length(f);
  y = zeros(size(x,1),nfuns);
  for j = 1:nfuns
      y(:,j) = f{j}(x(:,j));
  end
  
      

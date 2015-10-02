function [ff,splfuns] = fitSumOfSplines(Y,X,splineStruct);
% [ff,splfuns] =  fitSumOfSplines(Y,X,breaks,smoothness,extrapDeg);
% 
% Fit parameters for spline functions f1, f2, f3, ....
% in order to fit:  Y = f1(X(:,1)) + f2(X(:,2)) + f3(X(:,3) + ...
% via least-squares regression
%
% Inputs: Y - dependent variable (column vector)
%         X - indep variables (each column is a diff regressor)
%         splineStruct - structure with fields: "breaks", "smoothness", "extrapDeg"
%            - Use cell array if different params for each regressor
% 
% Outputs: ff - handle to function ff(x) = f1(x(:,1) + f2(x(:,2) + ...
%          splfuns - cell array of individual functions f1, f2, ....

[nx,nsplines] = size(X);

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
    if ispl == 1  % Keep DC for first function
        Mspline{ispl} = splineParamMatrix(ss.breaks, ss.smoothness, ss.extrapDeg);
    else  % Remove DC of first segment for all subsequent splines
        Mspline{ispl} = splineParamMatrix(ss.breaks, ss.smoothness, ss.extrapDeg,1);
    end
    nsegs(ispl) = length(ss.breaks)-1;  % number of segments
    breaks{ispl} = ss.breaks;
    nprs(ispl) = size(Mspline{ispl},2);

end

% 2. Build up design matrix
MM = zeros(nx,sum(nprs));   % design matrix for full system
ntotprs = 0;  % cumulative # of spline params inserted into MM
for ispl = 1:nsplines
   
    % Set spline structure for this column (regressor)
    ss = spstr{ispl};

    tic;
    % Sort data by first column
    [xsrt,isrt] = sort(X(:,ispl));
    
    % for each data point, compute its breakpoint interval
    [ignored,index] = sort([ss.breaks(1:nsegs(ispl)) xsrt']);
    index = max([find(index>nsegs(ispl))-(1:nx);ones(1,nx)])';

    % convert x values to local coordinates
    xsrt = xsrt-ss.breaks(index)';

    % Insert data into design matrix

    % ---- Use this version if memory limitations occur
    %AA = sparse([],[],[],nx,nsegs(ispl)*4,nx*4);  
    
    % ---- Use this version if no memory problems (much faster!)
    AA = zeros(nx,nsegs(ispl)*4);

    % ---- Loop through each segment, inserting points into design matrix
    for j = 1:nsegs(ispl)
        ii = [index==j];
        ni = sum(ii);
        if ni > 0
            AA(ii,((j-1)*4+1:j*4)) = ...
                repmat(xsrt(ii),1,4).^repmat([0:3],ni,1);
        end
    end
    MM(isrt,ntotprs+(1:nprs(ispl))) = AA*Mspline{ispl};
    ntotprs = ntotprs+nprs(ispl);
end

% 3. Solve for spline prs & make spline functions
prs = MM\Y;

% 4. Compute piecewise polynomial coeffs and create func handle
for ispl = 1:nsplines
    iprs = sum(nprs(1:ispl-1));
    splprs  = prs(iprs+[1:nprs(ispl)]);
    coeffs = fliplr(reshape(Mspline{ispl}*splprs,4,[])');
    pp = mkpp(breaks{ispl}, coeffs); % make piecewise polynomial
    splfuns{ispl} = @(x)ppval(pp,x);
end

ff = @(x)splineSum(splfuns,x);


% ------------------------------
function y = splineSum(f,x)
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
  y = zeros(size(x,1),1);
  for j = 1:nfuns
      y = y+f{j}(x(:,j));
  end
  
      

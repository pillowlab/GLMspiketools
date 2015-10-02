function varargout = makeSplineDesignMats(X,splineStruct);
% [M1, M2, ..., Mspline] =  makeSplineDesignMats(X,splineStruct);
% 
% Fit parameters for spline functions f1, f2, f3, ....
% in order to fit:  Y = f1(X(:,1)) + f2(X(:,2)) + f3(X(:,3) + ...
% via least-squares regression
%
% Inputs: 
%   X - indep variables (each column is a diff regressor)
%   splineStruct - structure with fields: "breaks", "smoothness", "extrapDeg"
%                - Use cell array if different params for each regressor
% 
% Outputs: M1, M2, .... 
%          Mspline - cell array of matrices mapping spline params to
%                    piecewise polynomial params

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
    MM = zeros(nx,nprs(ispl));
    MM(isrt,:) = AA*Mspline{ispl};
    varargout{ispl} = MM;
    ntotprs = ntotprs+nprs(ispl);
end

varargout{end+1} = Mspline;

function [Xdesign,Ysrt,Mspline] = mksplineDesignMat(X,Y,ss,sortflag,dcflag)
% [Xdesign,Ysrt,Mspline] = mksplineDesignMat(X,Y,ss,sortflag)
% 
% Computes design matrix for least-squares fitting cubic splines, so that
% the problem can be written:
%      argmin_w  ||Ysrt - Xdesign*w||^2
%
% "Standard" spline coefficients can be recovered via
%      splinecoeffs = reshape(Mspline*w,4,[])';
%
% Inputs: 
%   X - input (indep variable)
%   Y - output (dep variable)
%   ss - spline structure with fields: "breaks", "smoothness", "extrapDeg"
%   sortflag - set to zero to obtain *unsorted* output (default=1)
%   dcflag - set to zero to constrain f at first knot to zero (default=1)
% 
% Outputs: 
%   Xdesign - design matrix (spline matrix projected onto X monomials)
%   Ysrt - sorted Y so rows of Ysrt correspond to rows of Xdesign
%   Mspline - spline parametrization matrix

if nargin<4
    sortflag=1;
end
if nargin<5
    dcflag = 1; % set to zero to have no DC component
end

breaks = ss.breaks;        % number of knots
nsegs = length(breaks)-1;  % number of segments 
nx = size(X,1);            % # elements in X and Y

% Make spline parametrization matrix
Mspline = mksplineParamMtx(ss,dcflag);

% for each data point, compute its breakpoint interval
[xsrt,isrt] = sort(X); % sort X values
[~,index] = sort([breaks(1:nsegs) xsrt']);
index = max([find(index>nsegs)-(1:nx);ones(1,nx)])'; % segment for each datapoint
    
% convert to local coordinates
xsrt = xsrt-breaks(index)';

% Insert data into design matrix
% - Use this version if no memory problems (much faster!)
AA = zeros(nx,nsegs*4);
% - Use this version if memory limitations occur
%AA = sparse([],[],[],nx,nsegs(ispl)*4,nx*4);

% Loop through each segment, inserting points into design matrix
for j = 1:nsegs
    ii = [index==j]; % points with this index
    ni = sum(ii);    % number of such points
    if ni > 0
	AA(ii,((j-1)*4+1:j*4)) = ...
	    repmat(xsrt(ii),1,4).^repmat([3:-1:0],ni,1);
    end
end

Xdesign = AA*Mspline;

if (sortflag==0)
    % Pass back unsorted Design matrix
    Xdesign(isrt,:) = Xdesign;
    Ysrt = Y;
elseif nargout>1
    % Sort Y values to match X values
    Ysrt = Y(isrt);
end

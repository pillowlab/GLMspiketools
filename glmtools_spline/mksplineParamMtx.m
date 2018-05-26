function Mspline = mksplineParamMtx(ss,dcflag)
% Mspline = mksplineParamMtx(ss,dcflag)
%
% Computes matrix parametrizing a cubic spline (3rd order piecewise
% polynomial) with given smoothness and end-point conditions 
% 
% Inputs:  
%   ss - spline struct with fields "breaks", "smoothness", and "extrapDeg"
%     fields:
%      .breaks - x values where polynomials are spliced together 
%      .smoothness - number of derivs of continuity.  (optional; default=3)
%         0 - not smooth (discontinuous)
%         1 - continuous (but not differentiable)
%         2 - differentiable (continuous 1st derivatives)
%         3 - twice differentiable (continuous 2nd derivs)
%             (3 is the standard condition for cubic splines)
%      .extrapDeg - degree of polynomials used in edge segments in [0,3]
%         0 - constant
%         1 - straight line
%         2 - quadratic final polynomials
%         3 - 3rd order polynomials (standard in cubic splines)
%         Note: use a 2-vector to specify diff polynomial degrees on left
%         and right segments.  (optional; default=3)
%   dcflag = [0 or 1].  Set to 0 to constrain avg function value at knots to 0.
%
% Output:  Mspline = matrix parametrizing the given class of splines
%     This matrix is given by the null space of a matrix enforcing boundary
%     constraints on the piecewise polynomials.
% See also: fitSpline.m, makeSplineFun.m, convrtSplinePrs.m 

if nargin==1
    dcflag=1;  % default: no constraint on DC of first segment
end

breaks = ss.breaks(:)';  % Convert to row vector;
nbreaks = length(breaks);
if nbreaks<3
    error('Must use at least 3 breaks. (2 breaks => single polynomial)');
end
nsegs = nbreaks-1;  % # of polynomial segments
ncoeffs = nsegs*4;
smoothness = ss.smoothness;
extrapDeg = ss.extrapDeg;

dd = diff(breaks(1:end-1))'; % spacing beteen breaks
v0 = zeros(nsegs-1,1); % vectors of zeros and ones
v1 = ones(nsegs-1,1);

% Generate row vectors needed to enforce smoothness constraints
dffs = [];
dffs(:,:,1) = [dd.^3   dd.^2 dd v1,  v0  v0  v0 -v1]; % continuity
dffs(:,:,2) = [3*dd.^2 2*dd  v1 v0,  v0  v0 -v1 v0]; % differentiability
dffs(:,:,3) = [3*dd    v1    v0 v0,  v0 -v1  v0 v0]; % twice differentiability

ConstrMat = zeros(smoothness*(nsegs-1), 4*nsegs);  % Constraint matrix

% Generate constraints for smoothness
for i = 1:smoothness
    nrow = (i-1)*(nsegs-1)+1;
    for j = 1:nsegs-1
        ncol = (j-1)*4;
        ConstrMat(nrow,ncol+[1:8]) = dffs(j,:,i);
        nrow=nrow+1;
    end
end

% Generate extra constraints for degree of last polynomial 
meye = eye(4);
nleftC = 3-extrapDeg(1);  % # additional constr on left side
nrightC = 3-extrapDeg(end); % # additional constr on right side
epMat1 = [meye(1:nleftC,:) zeros(nleftC,ncoeffs-4)];
epMat2 = [zeros(nrightC,ncoeffs-4), meye(1:nrightC,:)];
ConstrMat = [ConstrMat; epMat1; epMat2]; 

% Generate constraint on DC component of first segment (if desired)
if dcflag==0
    %dcconstr = [0 0 0 1, zeros(1,ncoeffs-4)];
    dcconstr = repmat([0 0 0 1],1,nsegs);
    ConstrMat = [ConstrMat; dcconstr];
end

% size(ConstrMat,1) = total # of constraints
% Param Matrix is given by the null space of the constraint matrix
Mspline = null(ConstrMat);

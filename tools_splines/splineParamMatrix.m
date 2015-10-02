function Mspline = splineParamMatrix(knots, smoothness, extrapDeg,dcConstr);
% Mspline = splineParamMatrix(knots, smoothness, extrapDeg, dcConstr);
%
% Computes matrix parametrizing a cubic spline (3rd order piecewise
% polynomial) with given smoothness and end-point conditions 
% 
% Inputs:  
%   knots - x values where polynomials are spliced together
%   smoothness - number of derivs of continuity.  (optional; default=3)
%        0 - not smooth (discontinuous)
%        1 - continuous (but not differentiable)
%        2 - differentiable (continuous 1st derivatives)
%        3 - twice differentiable (continuous 2nd derivs)
%            (3 is the standard condition for cubic splines)
%   extrapDeg - degree of polynomials used in edge segments ( >=0, <=3)
%        0 - constant
%        1 - straight line
%        2 - quadratic final polynomials
%        3 - 3rd order polynomials (standard in cubic splines)
%        Note: use a 2-vector to specify diff polynomial degrees on left
%        and right segments.  (optional; default=3)
%    dcConstr = 0 (default) or 1 (set first section to lack a dc term).
%
% Output:  Mspline = matrix parametrizing the given class of splines
%     This matrix is given by the null space of a matrix enforcing boundary
%     constraints on the piecewise polynomials.
%
% See also: fitSpline.m, makeSplineFun.m, convrtSplinePrs.m 

if nargin == 1
    smoothness = 3;
    extrapDeg = [3 3];
elseif nargin == 2
    extrapDeg = [3 3];
end

% Flag for setting DC of first piece to zero
if nargin <= 3
    dcConstr = 0;
else
    dcConstr = 1;
end

if length(extrapDeg) == 1
    extrapDeg = extrapDeg*[1 1];
end

knots = knots(:)';  % Convert to row vector;
nknots = length(knots);
nsegs = nknots-1;  % # of polynomial segments
ncoeffs = 4*nsegs;  % # parametrizing a collection of cubic funcs

dd = diff(knots(1:end-1))'; % spacing beteen knots
v0 = zeros(nsegs-1,1); % vectors of zeros and ones
v1 = ones(nsegs-1,1);

dffs = []; % values needed for enforcing smoothness constraints
dffs(:,:,1) = [v1 dd dd.^2 dd.^3 -v1 v0 v0];
dffs(:,:,2) = [v0 v1 2*dd 3*dd.^2 v0 -v1 v0];
dffs(:,:,3) = [v0 v0 v1 3*dd v0 v0 -v1];

ConstrMat = zeros(smoothness*(nsegs-1), 4*nsegs);  % Constraint matrix
% Generate constraints for smoothness
for i = 1:smoothness
    nrow = (i-1)*(nsegs-1);
    for j = 1:nsegs-1
        ncol = (j-1)*4;
        nrow=nrow+1;
        ConstrMat(nrow,ncol+[1:7]) = dffs(j,:,i);
    end
end
[nconstr,ncoeffs] = size(ConstrMat);

% Generate constraints for end-points
meye = fliplr(eye(4));
nleftC = 3-extrapDeg(1);  % # additional constr on left side
nrightC = 3-extrapDeg(2); % # additional constr on right side
epMat1 = [meye(1:nleftC,:) zeros(nleftC,ncoeffs-4)];
epMat2 = [zeros(nrightC,ncoeffs-4), meye(1:nrightC,:)];

ConstrMat = [ConstrMat; epMat1; epMat2];
if dcConstr
    ConstrMat(end+1,:) = [1,zeros(1,ncoeffs-1)];
end
% nconstr = size(ConstrMat,1);  % Total # of constraints

% Param Matrix given by the null space of the constraint matrix
Mspline = null(ConstrMat);  

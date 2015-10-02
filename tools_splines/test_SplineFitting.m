% test_SplineFitting1.m
%
% Simple script illustrating low-d parametrization of a cubic spline and
% MSE-based fitting


% 1. Set function to fit
f = @(x)sin(x*2*pi);

% 2. Generate some noisy data 
nsamps = 200;  % number of points
signse = .2;
x = rand(nsamps,1)*2-1;
y = f(x)+randn(nsamps,1)*signse;

subplot(221);
xvals = -1:.01:1;
plot(x, y, '.', xvals, f(xvals), 'k');
title('raw data');

% 3. Set up cubic spline parameters
dbreak = .25;
breaks = -1:dbreak:1; % points of discontinuity
smoothness = 3;  % 1+# of times differentiable at breaks
extrapDeg = [3,2]; % degree polynomial on each end segment 
Mspline = splineParamMatrix(breaks, smoothness, extrapDeg);


%% 4. Now do min-MSE fitting ====================

[xsrt,isrt] = sort(x);  % Sort data by x value
ysrt = y(isrt);

nx = length(xsrt);  % Number of elements in x.
nsegs = length(breaks)-1;  % number of segments 

% for each data point, compute its breakpoint interval
[ignored,index] = sort([breaks(1:nsegs) xsrt']);
index = max([find(index>nsegs)-(1:nx);ones(1,nx)])';
    
% convert go to local coordinates
xsrt = xsrt-breaks(index)';

% Insert data into design matrix
AA = sparse(nx,nsegs);
for j = 1:nsegs
    ii = find(index==j);
    ni = length(ii);
    AA(ii(1):ii(end),(j-1)*4+1:j*4) = repmat(xsrt(ii),1,4).^repmat([0:3],ni,1);
end

% Solve equation for spline params (in column space of Mspline).
splinePrs = (AA*Mspline)\ysrt;

%% 5. Put coeffs into a 'pp' (piecewise polynomial) structure and function
coeffs = fliplr(reshape(Mspline*splinePrs,4,[])');
pp = mkpp(breaks, coeffs); % make piecewise polynomial
splfun = @(x)ppval(pp,x);


% 6. Make Plots: examine fit.
subplot(222);
plot(xvals, f(xvals), 'k', xvals, splfun(xvals), 'r', ...
    breaks, splfun(breaks), 'ro');
title('true fun vs. spline fit');

xvals2 = -1.2:.05:1.2; % Examine extrapolation
subplot(224);
plot(xvals2, f(xvals2), 'k', xvals2, splfun(xvals2), 'r', ...
    breaks, splfun(breaks), 'ro');
title('Extrapolation behavior');
axis tight;

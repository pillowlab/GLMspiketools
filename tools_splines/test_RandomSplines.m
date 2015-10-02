% test code by making random splines

knots = [0:1:10];  % Set points where polynomials spliced together
xx = knots(1)-1:.01:knots(end)+1; % Fine mesh for examining interpolation

smoothness = 3;  % 1 + # of continuous derivs   (in [0,3])
extrapDeg = [2 0]; % degree of edge polynomials (in [0,3])
M = splineParamMatrix(knots, smoothness, extrapDeg); % generate param matrix

% Generate random cubic spline
nprs = size(M,2);  % # parameters needed to specify this spline
x = M*randn(nprs,1);
f = makeSplineFun(knots,x); % returns function handle to spline

% Plot knots and spline function
plot(knots, f(knots), 'o-', xx, f(xx), 'r');
title(sprintf('Smooth=%d;  L=deg%d, Right=deg%d',...
    smoothness,extrapDeg(1),extrapDeg(2)));


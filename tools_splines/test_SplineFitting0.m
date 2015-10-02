% Test fitting a known function, evaluated at a set of points, with a
% spline, by minimizing mean-squared error of the approximation

xxknots = -2:1:2;  % Points of discontinuity 
f0 = @exp;  % Function to approximate
yyknots = f0(xxknots); % True value of function at knots

xxtest = xxknots(1)-2:.01:xxknots(end)+2;  % test values
yytest = f0(xxtest); % true function values at test values

% Set parameters for spline
smoothness = 2;  % 1+# of times differentiable at knots
extrapDeg = [1,3]; % degree polynomial on each end segment 
   % ( left is straight line, right is 3rd deg polynomial).

% X and Y values for fitting 
xxfit = xxknots(1)-2:.1:xxknots(end)+2;
yyfit = f0(xxfit);
[fun,Mspline,prs] = fitSpline(xxknots,xxfit,yyfit,smoothness,extrapDeg); % Fit spline
yyPred1 = fun(xxtest);

% Compare matlab's cubic spline interpolation using the same knots
% (i.e., same # of polynomial pieces)
yyPred2 = spline(xxknots, yyknots, xxtest); % Matlab's spline fitting

% ------- Plot comparisons ------------
clf;
plot(xxknots, yyknots, 'o', xxtest, yytest, 'b', ...
    xxtest, yyPred1, 'r', xxtest, yyPred2, 'k');
legend('knots', 'true func',  'spline approx', 'Matlab spline',...
'location', 'northwest');
title('Function comparisons');
axis tight;

function coeffs = convrtSplinePrs(x);
% coeffs = convrtSplinePrs(x);
% 
% Convert spline params to those used by 'mkpp' (make piece-wise polynomial).

coeffs = fliplr(reshape(x,4,[])');


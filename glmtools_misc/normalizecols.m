function B = normalizecols(A);
%
%  B = normalizecols(A);
%
%  Normalizes the columns of a matrix, so each is a unit vector.

B = A./repmat(sqrt(sum(A.^2)), size(A,1), 1);
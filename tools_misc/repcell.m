function B = repcell(A, m)
% B = repcell(A, m)
%  
%  Like repmat, but places multiple copies of a matrix in a cell array;

for j = 1:m
    B{j} = A;
end

function M = repdiag(X,nreps);
% M = repdiag(X,nreps);
%
% Constructs a block-diagonal matrix using 'nreps' repeats of the matrix X.
% (Similar to blkdiag, but doesn't require passing X multiple times).
%
% Creates full (sparse) output if X is full (sparse)

[ny,nx] = size(X);


if issparse(X)
    M = X;
    for j = 2:nreps
        M = [[M;sparse(ny,nx*(j-1))] , [sparse(ny*(j-1),nx); X]];
    end
else
    M = zeros(ny*nreps,nx*nreps);
    for j = 1:nreps;
        M(ny*(j-1)+1:ny*j,nx*(j-1)+1:nx*j) = X;
    end
end


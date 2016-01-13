function B = blkdiagrep(A, mtimes)
% B = blkdiagrep(A, mtimes)
%  
% Forms a block-diagonal matrix with m blocks formed from copies of A 
% (uses kron)

B = kron(speye(mtimes),A);
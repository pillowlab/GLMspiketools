function Y = sameconv(X, F)
%  SAMECONV - causally filters X with F and returns output of same height
%
%  Y = sameconv(X, F)
%
%  F is not "flipped", so first element of Y is last row of F dot product with
%  first row of X. 
%  
%  Inputs: 
%     X [NxM] - tall matrix
%     F [RxM] - short matrix (same width as X)
%
%  Output: 
%     Y [Nx1] - X filtered causally with F.  
%               (first element is F(end,:)*X(1,:)' )
%
%  Attempts to use fastest method based on length and width of F. 
%  (i.e., using either matrix arithmetic or matlab's "conv2").
%
%  ex) sameconv(+- 0 0 -+, +- 1 4 -+) -->  +- 0 -+
%               |  1 0  |  |  2 5  |       |  3  |
%               |  0 0  |  +- 3 6 -+       |  2  |
%               |  0 0  |                  |  1  |
%               |  0 0  |                  |  0  |
%               |  0 1  |                  |  6  |
%               +- 0 0 -+                  +- 5 -+
%
%  (updated:  14/03/2011 JW Pillow)
% 
% Copyright 2011 Pillow Lab.
% $Id$

[nx, xwid] = size(X);
[nf, fwid] = size(F);
assert(xwid == fwid, 'X and F must have same second dimension');

% Decide which method to use based on size of F
if log2(nf)+3>2*log2(fwid)
    % Use CONV2 if F is tall and skinny
	
    % conv2 doesn't take sparse matrix
    if issparse(X); X = full(X); end
    if issparse(F); F = full(F); end

    Y = conv2([zeros(nf-1,xwid); X], rot90(F,2), 'valid');
else
    % Use Matrix Operations if F is short and fat.
    ny = nx+nf-1;  % # of columns in full convolution

    Y = zeros(ny,1);
    yy = X*F';
    for j = 1:nf
        Y(j:j+nx-1) = Y(j:j+nx-1) + yy(:,nf-j+1);
    end
    % Extract only the first part
    Y = Y(1:nx,:);
end

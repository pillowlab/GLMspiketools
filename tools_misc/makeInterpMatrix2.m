function M = makeInterpMatrix2(slen, nbn);
% M = makeInterpMatrix2(slen, dt);
%
%  Make (sparse matrix) for interpolation.  
%
% (second arg can be # bins or bin size).
%
%  Matrix performs the convolution equivalent to "fastinterp2.c";

if nbn<1
    nbn = round(1./nbn);
end

c1 = [1./nbn:1./nbn:1]';
cc = [c1; flipud(c1)-1/nbn];
M = spalloc(nbn*slen,slen, 2*nbn*slen);
phi = floor(nbn/2);
M(1:nbn+phi)= cc(phi+1:end);
for j = 2:slen-1
    M((j-1)*nbn-phi+1:(j+1)*nbn-phi,j) = cc;
end
M(nbn*slen-nbn-phi+1:nbn*slen,slen) = cc(1:nbn+phi);



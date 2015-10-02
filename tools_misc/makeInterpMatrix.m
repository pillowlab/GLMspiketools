function M = makeInterpMatrix(slen, nbn);
% M = makeInterpMatrix(slen, dt);
%
%  Make (sparse matrix) for interpolation.  
%
% (second arg can be # bins or bin size).

if nbn<1
    nbn = round(1./nbn);
end

c1 = [1./nbn:1./nbn:1]';
cc = [c1; flipud(c1)-1/nbn];
M = spalloc(nbn*slen,slen, 2*nbn*slen);
for j = 1:slen-1
    M((j-1)*nbn+1:(j+1)*nbn,j) = cc;
end
M(nbn*slen-nbn+1:nbn*slen,slen) = c1;



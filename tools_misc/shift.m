function xshift = shift(x, k)
%% SHIFT(x, k)
%%   shifts a matrix by k places (vertically downward) 
%% To induce horizontal (rightward) shift in matrices, use shift(x')'

[m, n] = size(x);  rowvect = 0;  
if (m==1)  % take transpose if we have a row vector
   x = x'; m = n; rowvect = 1;
end

s1 = mod(m-k,m); s2 = mod(k,m);
if (s1 ~= 0)   
  xshift(1:s2,:) = x((s1+1):m,:);
  xshift((s2+1):m,:) = x(1:s1,:);
else
  xshift = x;
end
  
if (rowvect), xshift = xshift';   end
   

% Old method:  make shifting matrix.  
%if (s1 ~= 0)
%  Mshift = diag(ones(s1,1),-s2)+diag(ones(s2,1),s1);
%  xshift = Mshift*x;
%else
%  xshift = x;
%end
   

function [ii, nds] = in(v, rnge)
% [ii,nds] = in(v, rnge)
% 
% returns indices of vector v within range rnge, inclusive of endpoints

if (nargout == 1)
    if length(rnge) == 1;
        ii = v(v<=vrnge);
    else
        ii = v((v>=rnge(1)) & (v<=rnge(2)));
    end

else
    if length(rnge) == 1;
        nds = find(v<=vrnge);
        ii = v(nds);
    else
        nds = find((v>=rnge(1)) & (v<=rnge(2)));
        ii = v(nds);
    end
end

function gg = reinsertFitPrs_GLMbi(gg,prs,Xstruct)
% gg = reinsertFitPrs_GLMbi(gg,prs);
%
% After fitting, reinsert params into param structure

% Put returned vals back into param structure ------
krank = Xstruct.krank;
nktprs = Xstruct.nkt*krank;
nkxprs = Xstruct.nkx*krank;
nktot = nktprs+nkxprs;
nh = Xstruct.nh; 
nh2 = Xstruct.nh2;

% Insert params into struct
gg.kt = reshape(prs(1:nktprs),[],krank);
gg.kx = reshape(prs(nktprs+1:nktprs+nkxprs),[],krank);
gg.dc = prs(nktot+1);
gg.ihw = reshape(prs(nktot+2:nktot+nh+1), nh,1);
gg.ihw2 = reshape(prs(nktot+nh+2:end), nh2, []);
gg.k = (gg.ktbas*gg.kt)*(gg.kxbas*gg.kx)';

if isempty(gg.ihw)
    gg.ihbas = [];
end
if isempty(gg.ihw2)
    gg.ihbas2 = [];
end
gg.ih = [gg.ihbas*gg.ihw, gg.ihbas2*gg.ihw2];


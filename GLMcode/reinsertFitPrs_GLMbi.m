function gg = reinsertFitPrs_GLMbi(gg,prs);
% gg = reinsertFitPrs_GLMbi(gg,prs);
%
% After fitting, reinsert params into param structure

global OPRS

% Put returned vals back into param structure ------
nkt = OPRS.nkt;
nkx = OPRS.nkx;
nktot = nkt+nkx;
nh = OPRS.nh; 
nh2 = OPRS.nh2;
krank = OPRS.krank;

% Insert params into struct
gg.kt = reshape(prs(1:nkt),[],krank);
gg.kx = reshape(prs(nkt+1:nkt+nkx),[],krank);
gg.dc = prs(nktot+1);
gg.ih = reshape(prs(nktot+2:nktot+nh+1), nh,1);
gg.ih2 = reshape(prs(nktot+nh+2:end), nh2, []);
gg.k = (gg.ktbas*gg.kt)*(gg.kxbas*gg.kx)';

function gg = reinsertFitPrs_GLM(gg,prs);
% gg = reinsertFitPrs_GLM(gg,prs);
%
% After fitting, reinsert params into param structure

global OPRS

nt = OPRS.nt;  
nx = OPRS.nx;   
nktot = nt*nx;
nh = OPRS.nh;  
nh2 = OPRS.nh2;

gg.kt = reshape(prs(1:nktot),nt,nx);
gg.k = gg.ktbas*gg.kt;
gg.dc = prs(nktot+1);
gg.ih = reshape(prs(nktot+2:nktot+nh+1), nh,1);
gg.ih2 = reshape(prs(nktot+nh+2:end), nh2, []);

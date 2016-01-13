function gg = reinsertFitPrs_GLMbimx(gg,prs);
% gg = reinsertFitPrs_GLMbimx(gg,prs);
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

% Find relevant indices for spatial filters
nx = OPRS.nx;
xwids = gg.xwids;
krank = gg.krank;
msk = zeros(nx,sum(krank));
icum = 0;
jcum = 0;
for jj = 1:length(krank)
    msk(icum+1:icum+xwids(jj),jcum+1:jcum+krank(jj)) = 1;
    icum = icum+xwids(jj);
    jcum = jcum+krank(jj);
end
ii = find(msk);

% Insert into param object -------------
gg.kt = reshape(prs(1:nkt),[],sum(krank));
gg.kx = msk*0;
gg.kx(ii) = prs(nkt+1:nkt+nkx);
gg.dc = prs(nktot+1);
gg.ih = reshape(prs(nktot+2:nktot+nh+1), nh,1);
gg.ih2 = reshape(prs(nktot+nh+2:end), nh2, []);
gg.k = (gg.ktbas*gg.kt)*(gg.kxbas*gg.kx)';

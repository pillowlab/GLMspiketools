function gg = reinsertFitPrs_GLM(gg,prs);
% gg = reinsertFitPrs_GLM(gg,prs);
%
% After fitting, reinsert params into param structure

global OPTprs

nkt = OPTprs.nkt;  
nkx = OPTprs.nkx;   
nktot = nkt*nkx;
nh = OPTprs.nh;  
nh2 = OPTprs.nh2;

gg.kt = reshape(prs(1:nktot),nkt,nkx);
gg.k = gg.ktbas*gg.kt;
gg.dc = prs(nktot+1);
gg.ihw = reshape(prs(nktot+2:nktot+nh+1), nh,1);
gg.ihw2 = reshape(prs(nktot+nh+2:end), nh2, []);

if isempty(gg.ihw)
    gg.ihbas = [];
end
if isempty(gg.ihw2)
    gg.ihbas2 = [];
end
gg.ih = [gg.ihbas*gg.ihw, gg.ihbas2*gg.ihw2];

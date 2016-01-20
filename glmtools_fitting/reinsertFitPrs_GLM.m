function gg = reinsertFitPrs_GLM(gg,prs,Xstruct)
% gg = reinsertFitPrs_GLM(gg,prs,Xstruct)
%
% After fitting, reinsert params into GLM param structure

% Extract relevant size information
nkt = Xstruct.nkt;  
nkx = Xstruct.nkx;   
nktot = nkt*nkx;
nh = Xstruct.nh;  
nh2 = Xstruct.nh2;

% Insert params
gg.kt = reshape(prs(1:nktot),nkt,nkx);
gg.k = gg.ktbas*gg.kt;
gg.dc = prs(nktot+1);
gg.ihw = reshape(prs(nktot+2:nktot+nh+1), nh,1);
gg.ihw2 = reshape(prs(nktot+nh+2:end), nh2, []);

% Ensure these bases are 'empty' if no spike history filters
if isempty(gg.ihw), gg.ihbas = [];
end
if isempty(gg.ihw2), gg.ihbas2 = [];
end

% Insert spike-history filters
gg.ih = [gg.ihbas*gg.ihw, gg.ihbas2*gg.ihw2];

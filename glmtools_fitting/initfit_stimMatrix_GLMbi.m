function initfit_stimMatrix_GLMbi(gg, Stim)
%  initfit_stimMatrix_GLMbi(gg, Stim);
%
% Initialize parameters relating to stimulus design matrix

global OPTprs

% ---- Set up filter and stim processing params ------------------- 
nkx = size(gg.kxbas,2); % # params per spatial vector
nkt = size(gg.ktbas,2); % # time params per stim pixel (# t params)
krank = gg.krank;     % rank of filter 
ncols = nkx*nkt;

[slen,swid] = size(Stim);
rlen = round(slen/gg.dt); % total length (fine bins)

% ---- Filter stimulus with spatial and temporal bases -----
xfltStim = Stim*gg.kxbas;
OPTprs.MSTM = zeros(slen,ncols);
for i = 1:nkx
    for j = 1:nkt
        OPTprs.MSTM(:,(i-1)*nkt+j) = sameconv(xfltStim(:,i),gg.ktbas(:,j));
    end
end

% ---- Set up stim processing params ----------
OPTprs.nkx = nkx;
OPTprs.nkt = nkt;
OPTprs.krank = krank;
OPTprs.ktbas = gg.ktbas;
OPTprs.kxbas = gg.kxbas;
OPTprs.slen = slen;  % Total stimulus length (course bins)
OPTprs.rlen = rlen;
OPTprs.dt = gg.dt;



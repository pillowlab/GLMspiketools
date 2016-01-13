function initfit_stimMatrix(gg,Stim)
% initfit_stimMatrix(gg,Stim)
%  
% Initialize parameters relating to stimulus design matrix 

global OPTprs 

% ---- Set up filter and stim processing params ------------------- 
nkx = size(gg.k,2);  % number stim pixels (# x params)
nkt = size(gg.ktbas,2); % # time params per stim pixel (# t params)
ncols = nkx*nkt;

[slen,swid] = size(Stim);
rlen = round(slen/gg.dt); % total length (fine bins)

% ---- Check size of filter and width of stimulus ----------
if (nkx ~= swid)
    error('Mismatch between stim width and kernel width');
end

% ---- Filter stimulus with spatial and temporal bases -----
OPTprs.MSTM = zeros(slen,ncols);
for i = 1:nkx
    for j = 1:nkt
        OPTprs.MSTM(:,(i-1)*nkt+j) = sameconv(Stim(:,i),gg.ktbas(:,j));
    end
end

% ---- Set fields of OPTprs -------------------------------------
OPTprs.nkx = nkx;
OPTprs.nkt = nkt;
OPTprs.ktbas = gg.ktbas;
OPTprs.slen = slen;  % Total stimulus length (course bins)
OPTprs.rlen = rlen;
OPTprs.dt = gg.dt;

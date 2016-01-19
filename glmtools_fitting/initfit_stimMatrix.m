function Xstruct = initfit_stimMatrix(gg,Stim)
% OPTprs = initfit_stimMatrix(gg,Stim)
%  
% Initialize parameters relating to stimulus design matrix 

% ---- Set up filter and stim processing params ------------------- 
nkx = size(gg.k,2);  % number stim pixels (# x params)
nkt = size(gg.ktbas,2); % # time params per stim pixel (# t params)
ncols = nkx*nkt; % total number of columns in stim design matrix

[slen,swid] = size(Stim);
upsampfactor = (gg.dtStim/gg.dtSp); % number of times by which spikes more finely sampled than stim
rlen = slen*upsampfactor;

% ---- Check size of filter and width of stimulus ----------
if (nkx ~= swid)
    error('Mismatch between stim width and kernel width');
end

% ---- Filter stimulus with spatial and temporal bases -----
Xstruct.Xstim = zeros(slen,ncols);
for i = 1:nkx
    for j = 1:nkt
        Xstruct.Xstim(:,(i-1)*nkt+j) = sameconv(Stim(:,i),gg.ktbas(:,j));
    end
end

% ---- Set fields of OPTprs -------------------------------------
Xstruct.nkx = nkx;
Xstruct.nkt = nkt;
Xstruct.ktbas = gg.ktbas;
Xstruct.slen = slen;  % Total stimulus length (course bins)
Xstruct.rlen = rlen;
Xstruct.upsampfactor = upsampfactor;
Xstruct.dtStim = gg.dtStim;
Xstruct.dtSp = gg.dtSp;


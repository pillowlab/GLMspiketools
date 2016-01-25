function Xstruct = initfit_stimDesignMat(gg,Stim)
% Xstruct = initfit_stimDesignMat(gg,Stim)
%  
% Initialize parameters relating to stimulus design matrix 

% ---- Set up filter and stim processing params ------------------- 
nkx = size(gg.k,2);  % number stim pixels (# x params)
nkt = size(gg.ktbas,2); % # time params per stim pixel (# t params)
ncols = nkx*nkt; % total number of columns in stim design matrix

[slen,swid] = size(Stim); % size of stimulus
upsampfactor = (gg.dtStim/gg.dtSp); % number of times by which spikes more finely sampled than stim
rlen = slen*upsampfactor; % length of spike train vector 

% ---- Check size of filter and width of stimulus ----------
assert(nkx == swid,'Mismatch between stim width and kernel width');

% ---- Convolve stimulus with spatial and temporal bases -----
Xstruct.Xstim = zeros(slen,ncols);
for i = 1:nkx
    for j = 1:nkt
        Xstruct.Xstim(:,(i-1)*nkt+j) = sameconv(Stim(:,i),gg.ktbas(:,j));
    end
end

% ---- Set fields of Xstruct -------------------------------------
Xstruct.nkx = nkx; % # stimulus spatial stimulus pixels 
Xstruct.nkt = nkt; % # time bins in stim filter
Xstruct.slen = slen;  % Total stimulus length (coarse bins)
Xstruct.rlen = rlen;  % Total spike-train bins (fine bins)
Xstruct.upsampfactor = upsampfactor; % rlen / slen
Xstruct.Minterp = kron(speye(slen),ones(upsampfactor,1)); % rlen x slen matrix for upsampling and downsampling 
Xstruct.dtStim = gg.dtStim; % time bin size for stimulus
Xstruct.dtSp = gg.dtSp; % time bins size for spike train


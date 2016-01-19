function [neglogli,rr,tt,Iinj,Istm,Ih,Icpl] = neglogli_GLM(gg,Stim)
% [neglogli,rr,tt,Iinj,Istm,Ih,Icpl] = neglogli_GLM(gg,Stim)
% 
% Compute glm model negative log-likelihood given the parameters in gg, 
%  
% Inputs: gg = param object
%         fields: k -  stimulus kernel
%                 ih - post-spike current
%                 dc - dc current injection
%                 nlfun - nonlinearity
%                 dt - time bin size
%         Stim = stimulus
%
% Outputs:
%     neglogli = negative log-likelihood of spike trains
%     rr = conditional intensity (in expected spikes /sec)
%     tt = time bins for for conditional intensity
%     Iinj = net linear input (sum of all filter outputs)
%     Istm = net linear input from stimulus
%     Ih = net linear input from own spike-history
%     Icpl = matrix of coupling current inputs 

% ---- Extract params from glm struct ---------------
k = gg.k;
dc = gg.dc;
dt = gg.dtSp;  % time bin size for spikes
upSampFactor = gg.dtStim/dt; % number of spike bins per Stim bin
slen = size(Stim,1);
rlen = size(gg.sps,1);

% Check that upSampFactor is an integer
assert(mod(upSampFactor,1) == 0, 'dtStim / dtSp must be an integer');
% Check factor relating size of stim and binned spike train
assert(slen*upSampFactor==rlen,'Spike train length must be an even multiple of stim length');

% ----  Compute filtered resp to stimulus -----------------------------
I0 = sameconv(Stim,k);
Iinj = kron(I0,ones(upSampFactor,1)) + dc;

% -------------- Compute net h current --------------------------------
spInds = find(gg.sps);

% Check if post-spike filters are present
if isempty(gg.ihw), gg.ihbas = []; end
if isempty(gg.ihw2), gg.ihbas2 = []; end

% Compute post-spike filter from basis and weights
ih = [gg.ihbas*gg.ihw, gg.ihbas2*gg.ihw2]; % h current
nCoupled = length(gg.couplednums); % # cells coupled to this one

if ~isempty(ih)
    % self-coupling filter
    Iinj = Iinj + spikefilt_mex(spInds,ih(:,1),[1,rlen]);
    % coupling filters from other neurons
    for j = 1:nCoupled
            Iinj = Iinj + spikefilt_mex(find(gg.sps2(:,j)),ih(:,j+1),[1,rlen]);
    end
end

rr = gg.nlfun(Iinj);  % Conditional intensity

% -- Compute negative log-likelihood from conditional intensity ------
bmask = initfit_mask(gg.mask,dt,rlen);  % bins to use for likelihood calc
trm1 = sum(rr(bmask))*dt;  % non-spiking term
trm2 = -sum(log(rr(gg.sps(bmask)>0)));  % Spiking term
neglogli = trm1 + trm2;


% ======  OPTIONAL OUTPUT ARGS ===================================
if nargout > 2 
    tt = (1:rlen)'*dt; % time indices 
end

if nargout > 4  % Return details of currents
    Istm = kron(I0,ones(upSampFactor,1)) + dc;
end

if nargout > 5
    if ~isempty(ih)
        Ih = spikefilt_mex(spInds,ih(:,1),[1,rlen]);
    else
        Ih = zeros(length(Iinj),1);
    end
end

if nargout > 6
    Icpl = zeros(length(Iinj),size(ih,2)-1);
    if ~isempty(ih)
        for j = 1:nCoupled
            Icpl(:,j) = spikefilt_mex(spInd2{j},ih(:,j+1),[1,rlen]);
        end
    end
end

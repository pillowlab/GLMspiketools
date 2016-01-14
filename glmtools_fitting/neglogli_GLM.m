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


% variables for spike times and frame rate
% global RefreshRate SPNDS SPNDS2;

% ---- Extract params from glm struct ---------------
k = gg.k;
dc = gg.dc;
dt = gg.dtSp;  % time bin size for spikes
upSampFactor = gg.dtStim/dt; % number of spike bins per Stim bin
assert(mod(upSampFactor,1) == 0, 'dtStim / dtSp must be an integer');


% ----  Compute filtered resp to stimulus -----------------------------
I0 = sameconv(Stim,k);
Iinj = kron(I0,ones(upSampFactor,1)) + dc;
rlen = length(Iinj);

% -------------- Compute net h current --------------------------------
[spInds,spInd_cpld] = initfit_spikeInds(gg); % get spike bin indices

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
            Iinj = Iinj + spikefilt_mex(spInd_cpld{j},ih(:,j+1),[1,rlen]);
    end
end

rr = gg.nlfun(Iinj);  % Conditional intensity
[iiSpk,iiLi] = initfit_mask(gg.mask,dt,rlen);  % bins to use for likelihood calc

% -- Compute negative log-likelihood from conditional intensity ------
trm1 = -sum(log(rr(iiSpk)));  % Spiking term
trm2 = sum(rr(iiLi))*dt;  % non-spiking term
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
            Icpl(:,j) = spikefilt_mex(spInd_cpld{j},ih(:,j+1),[1,rlen]);
        end
    end
end

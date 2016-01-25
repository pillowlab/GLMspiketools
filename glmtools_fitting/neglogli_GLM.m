function [neglogli,rr,tt,Itot,Istm,Ih,Icpl] = neglogli_GLM(gg,Stim)
% [neglogli,rr,tt,Itot,Istm,Ih,Icpl] = neglogli_GLM(gg,Stim)
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
%   neglogli = negative log-likelihood of spike trains
%         rr = conditional intensity (in expected spikes /sec)
%         tt = time bins for for conditional intensity
%       Itot = sum of all filter outputs
%       Istm = sum of stim filter output and dc term
%         Ih =   output of neuron's own spike-history filter
%       Icpl = matrix of coupling current inputs 

% ---- Extract params from glm struct ---------------
k = gg.k;
dc = gg.dc;
ih = gg.ih;
dt = gg.dtSp;  % time bin size for spikes
upsampfactor = gg.dtStim/dt; % number of spike bins per Stim bin
slen = size(Stim,1);
rlen = size(gg.sps,1);

% Check that upSampFactor is an integer
assert(mod(upsampfactor,1) == 0, 'dtStim / dtSp must be an integer');
% Check factor relating size of stim and binned spike train
assert(slen*upsampfactor==rlen,'Spike train length must be an even multiple of stim length');

% ----  Compute filtered resp to stimulus -----------------------------
I0 = sameconv(Stim,k);
Itot = kron(I0,ones(upsampfactor,1)) + dc;

% -------------- Compute net h current --------------------------------

% Check if post-spike filters are present
if isempty(gg.ihw), gg.ihbas = []; end
if isempty(gg.ihw2), gg.ihbas2 = []; end

% Compute convolution of post-spike filters with spike history
nCoupled = length(gg.couplednums); % # cells coupled to this one
if ~isempty(ih)
    Itot = Itot + spikefilt(gg.sps,ih(:,1));  % self-coupling filter
    for j = 1:nCoupled   % coupling filters from other neurons
            Itot = Itot + spikefilt(gg.sps2(:,j),ih(:,j+1));
    end
end

rr = gg.nlfun(Itot);  % Conditional intensity

% -- Compute negative log-likelihood from conditional intensity ------
bmask = initfit_mask(gg.mask,dt,rlen);  % bins to use for likelihood calc
trm1 = sum(rr(bmask))*dt;  % non-spiking term
trm2 = -sum(log(rr(gg.sps(bmask)>0)));  % Spiking term
neglogli = trm1 + trm2;

% ======  OPTIONAL OUTPUT ARGS ===================================

% time indices 
if nargout > 2, tt = (1:rlen)'*dt; 
end

% Stim filter output
if nargout > 4, Istm = kron(I0,ones(upsampfactor,1)) + dc;
end

% Spike-history filter output
if nargout > 5
    Ih = zeros(length(Itot),1);
    if ~isempty(ih), Ih = spikefilt(gg.sps,ih(:,1));
    end
end

% Coupling filter outputs
if nargout > 6 
    Icpl = zeros(length(Itot),size(ih,2)-1);
    if ~isempty(ih)
        for j = 1:nCoupled
            Icpl(:,j) = spikefilt(gg.sps2(:,j),ih(:,j+1));
        end
    end
end

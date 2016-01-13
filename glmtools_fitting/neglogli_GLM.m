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
global RefreshRate SPNDS SPNDS2;

initfit_SPNDS(gg);  % set SPNDS and SPNDS2

% ---- Extract params from gg ---------------
k = gg.k;
dc = gg.dc;
dt = gg.dt;

% ----  Compute filtered resp to stimulus -----------------------------
I0 = sameconv(Stim,k);
Iinj = fastinterp2(I0,round(1/dt)) + dc;
rlen = length(Iinj);

% -------------- Compute net h current --------------------------------
if isempty(gg.ihw)
    gg.ihbas = [];
end
if isempty(gg.ihw2);
    gg.ihbas2 = [];
end
ih = [gg.ihbas*gg.ihw, gg.ihbas2*gg.ihw2]; % h current
ih = interp_spikeFilt(ih,gg.iht,dt);   % interpolated h current
nCoupled = length(gg.couplednums); % # cells coupled to this one

if ~isempty(ih)
    Iinj = Iinj + spikefilt_mex(SPNDS,ih(:,1),[1,rlen]);
    for j = 1:nCoupled
            Iinj = Iinj + spikefilt_mex(SPNDS2{j},ih(:,j+1),[1,rlen]);
    end
end

rr = gg.nlfun(Iinj);  % Conditional intensity
[iiSpk,iiLi] = initfit_mask(gg.mask,dt,rlen);  % bins to use for likelihood calc

% ---- Compute negative log-likelihood from conditional intensity ------
trm1 = -sum(log(rr(iiSpk)));  % Spiking term
trm2 = sum(rr(iiLi))*dt/RefreshRate;  % non-spiking term
neglogli = trm1 + trm2;


% ======  OPTIONAL OUTPUT ARGS ===================================
if nargout > 2 
    tt = [dt:dt:size(Stim,1)]'/RefreshRate;  % time indices 
end

if nargout > 4  % Return details of currents
    Istm = fastinterp2(I0,round(1/dt)) + dc;
end

if nargout > 5
    if ~isempty(ih)
        Ih = spikefilt_mex(SPNDS,ih(:,1),[1,rlen]);
    else
        Ih = zeros(length(Iinj),1);
    end
end

if nargout > 6
    Icpl = zeros(length(Iinj),size(ih,2)-1);
    if ~isempty(ih)
        for j = 1:nCoupled
            Icpl(:,j) = spikefilt_mex(SPNDS2{j},ih(:,j+1),[1,rlen]);
        end
    end
end

function [neglogli,rr,tt,Iinj,Istm,Ih,Icpl] = neglogli_GLM(gg,Stim);
% [neglogli,rr,tt,Iinj,Istm,Ih,Icpl] = neglogli_GLM(gg,Stim);
% 
% Compute glm model negative log-likelihood given the parameters in gg, 
%  
% Inputs: gg = param object, must have fields 
%              k -  stimulus kernel
%              ih - post-spike current
%              dc - dc current injection
%              dt - time bin size on which likelihood computed
%         Stim = stimulus
%
% Outputs:
%   neglogli = negative log-likelihood of spike trains
%   rr = conditional intensity (in spikes /sec)
%   tt = time binning for for conditional intensity
%   Iinj = net linear input (sum of all filter outputs)
%   Istm = net linear input from stimulus
%   Ih = net linear input from own spike-history
%   Icpl = matrix of coupling current inputs 

% variables for linearly interpolating stim, spike currents
global MSP SPNDS SPNDS2 RefreshRate;

setSPNDS(gg);  % set global vars

% ---- Extract params from gg ---------------
k = gg.k;
dc = gg.dc;
dt = gg.dt;
ih = [];
if ~isempty(gg.ih)
    ih = gg.ihbas*gg.ih; % h current, interpolated correctly
end
if ~isempty(gg.ih2)
    if isempty(gg.ihbas2)
        % Use the self-coupling basis for other inputs as well
        ih = [ih gg.ihbas*gg.ih2];
    else
        ih = [ih gg.ihbas2*gg.ih2];
    end
end
if ~isempty(ih)
    ih = MSP*ih; % Interpolate h current to correct time sampling
end
if isfield(gg, 'nlfun')
    nlfun = gg.nlfun;
else
    nlfun = @exp;
end
nsp = length(gg.tsp);
nCoupled = length(SPNDS2);

% ----  Compute filtered resp to signal ----------------
I0 = sameconv(Stim,k);
Iinj = fastinterp2(I0,round(1/dt)) + dc;
rlen = length(Iinj);

% -------------- Compute net h current ---------------------------
if ~isempty(ih)
    Iinj = Iinj + spikeconv_mex(SPNDS,ih(:,1),[1,rlen]);
    for j = 1:nCoupled
            Iinj = Iinj + spikeconv_mex(SPNDS2{j},ih(:,j+1),[1,rlen]);
    end
end

% ---- Nonlinearity ------
rr = nlfun(Iinj);  % Conditional intensity

% ---- Compute negative logli from binned Cond Intensity ------
if (gg.tspi(1) == 1), istrt = 1;
else 
    istrt = SPNDS(gg.tspi(1)-1)+1;
end
trm1 = -sum(log(rr(SPNDS(gg.tspi(1):end))));  % Spiking term
trm2 = sum(rr(istrt:end))*dt/RefreshRate;  % non-spiking term
neglogli = trm1 + trm2;

if nargout > 2 
    tt = [dt:dt:size(Stim,1)]'/RefreshRate;  % time indices 
end

% ======  OPTIONAL OUTPUT ARGS ===================================
if nargout > 4  % Return details of currents
    Istm = fastinterp2(I0,round(1/dt)) + dc;
end
if nargout > 5
    if ~isempty(ih)
        Ih = spikeconv_mex(SPNDS,ih(:,1),[1,rlen]);
    else
        Ih = zeros(length(Iinj),1);
    end
end
if nargout > 6
    Icpl = zeros(length(Iinj),size(ih,2)-1);
    if ~isempty(ih)
        for j = 1:nCoupled
            Icpl(:,j) = spikeconv_mex(SPNDS2{j},ih(:,j+1),[1,rlen]);
        end
    end
end

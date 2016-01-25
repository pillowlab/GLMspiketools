function [tsp,sps,Itot,Istm] = simGLM(glmprs,Stim)
% [tsp,sps,Itot,Ispk] = simGLM(glmprs,Stim)
% 
% Compute response of glm to stimulus Stim.
%
% Uses time rescaling instead of Bernouli approximation to conditionally
% Poisson process
%
% Dynamics:  Filters the Stimulus with glmprs.k, passes this through a
% nonlinearity to obtain the point-process conditional intensity.  Add a
% post-spike current to the linear input after every spike.
%
% Input: 
%   glmprs - struct with GLM params, has fields 'k', 'h','dc' for params
%              and 'dtStim', 'dtSp' for discrete time bin size for stimulus
%              and spike train (in s).
%     Stim - stimulus matrix, with time running vertically and each
%              column corresponding to a different pixel / regressor.
% Output:
%   tsp - list of spike times (in s)
%   sps - binary matrix with spike times (at resolution dtSp).
%  Itot - summed filter outputs 
%  Istm - just the spike-history filter output

% Check nonlinearity (default is exponential)
if ~isfield(glmprs, 'nlfun'), 
    glmprs.nlfun = @exp;
end

upSampFactor = glmprs.dtStim/glmprs.dtSp; % number of spike bins per Stim bin
assert(mod(upSampFactor,1) == 0, 'dtStim / dtSp must be an integer');

% Determine which version to run 
if size(glmprs.k,3) > 1  % Run "coupled" GLM model if multiple cells present
    [tsp,sps,Itot,Istm] = simGLMcpl(glmprs,Stim,upSampFactor);
else   % Single cell simulation
    [tsp,sps,Itot,Istm] = simGLMsingle(glmprs,Stim,upSampFactor);
end


% ======================================================================
function [tsp,sps,Itot,Istm] = simGLMsingle(glmprs,Stim,upSampFactor)
% Sub-function within simGLM.m
%
% Simulates the GLM point process model for a single (uncoupled) neuron
% using time-rescaling
    
% --------------- Check Inputs ----------------------------------
nbinsPerEval = 100;  % Default number of bins to update for each spike
dt = glmprs.dtSp; % bin size for simulation

slen = size(Stim,1); % length of stimulus
rlen = slen*upSampFactor;  % length of binned spike response
hlen = size(glmprs.ih,1); % length of post-spike filter

% -------------  Compute filtered resp to signal ----------------
Istm = sameconv(Stim,glmprs.k);  % filter stimulus with k
Istm = kron(Istm,ones(upSampFactor,1)) + glmprs.dc; % upsample stim current
Itot = Istm; % total filter output
    
% --------------- Set up simulation dynamics variables ---------------
nsp = 0; % number of spikes
sps = zeros(rlen,1); % sparse spike time matrix
jbin = 1; % current time bin
tspnext = exprnd(1);  % time of next spike (in rescaled time)
rprev = 0;  % Integrated rescaled time up to current point

% --------------- Run dynamics ---------------------------------------
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
    rrnxt = glmprs.nlfun(Itot(iinxt))*dt; % Cond Intensity
    rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
    if (tspnext >= rrcum(end)) % No spike in this window
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
    else   % Spike!
        ispk = iinxt(find(rrcum>=tspnext, 1, 'first')); % time bin where spike occurred
        nsp = nsp+1;
        sps(ispk) = 1; % spike time
        mxi = min(rlen, ispk+hlen); % max time affected by post-spike kernel
        iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
        if ~isempty(iiPostSpk)
            Itot(iiPostSpk) = Itot(iiPostSpk)+glmprs.ih(1:mxi-ispk);
        end
        tspnext = exprnd(1);  % draw next spike time
        rprev = 0; % reset integrated intensity
        jbin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = jbin/nsp;
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
end
tsp = find(sps>0)*dt;

% ========================================================================
function [tsp,sps,Itot,Istm] = simGLMcpl(glmprs,Stim,upSampFactor)
% [tsp,sps,Itot,Istm] = simGLMcpl(glmprs,Stim,upSampFactor)
% 
% Compute response of (multi-cell) coupled-glm to stimulus Stim.
%
% Uses time rescaling to sample conditionally Poisson process
%
% Dynamics:  Filters the Stimulus with glmprs.k, passes this through a
% nonlinearity to obtain the point-process conditional intensity.  Add a
% post-spike current to the linear input after every spike.
%
% Sub-function within simGLM.m

% --------------- Check Inputs ----------------------------------
ncells = size(glmprs.k,3);
nbinsPerEval = 100;  % Default number of bins to update for each spike
dt = glmprs.dtSp;
% Check that dt evenly divides 1
if mod(1,dt) ~= 0
    dt = 1/round(1/dt);
    fprintf(1, 'glmrunmod: reset dtsim = %.3f\n', dt);
end
slen = size(Stim,1); % length of stimulus
rlen = slen*upSampFactor;  % length of binned response
ih = permute(glmprs.ih,[1 3 2]); % permute so all outgoing filters in a single plane
hlen = size(ih,1); % length of h filter

% -------------  Compute filtered resp to signal ----------------
Istm = zeros(slen,ncells);
for jj = 1:ncells
    Istm(:,jj) = sameconv(Stim,glmprs.k(:,:,jj));  % filter stimulus with k
end
Istm = bsxfun(@plus,Istm,glmprs.dc); % add DC term

% Compute finely sampled version
Istm = kron(Istm,ones(upSampFactor,1)); % upsample stim current
Itot = Istm; % total filter output

% --------------- Set up simulation dynamics variables ---------------
tsp(1,1:ncells) = {zeros(round(slen/25),1)};  % allocate space for spike times
nsp = zeros(1,ncells);
jbin = 1;  % counter ?
tspnext = exprnd(1,1,ncells); % time of next spike (in rescaled time) 
rprev = zeros(1,ncells); % Integrated rescaled time up to current point

% --------------- Run dynamics ---------------------------------------
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen); % Bins to update in this iteration
    nii = length(iinxt);  % Number of bins
    rrnxt = glmprs.nlfun(Itot(iinxt,:))*dt; % Cond Intensity
    rrcum = cumsum(rrnxt+[rprev;zeros(nii-1,ncells)],1);  % Cumulative intensity
    if all(tspnext >= rrcum(end,:)) % No spike in this window
        jbin = iinxt(end)+1;
        rprev = rrcum(end,:);
    else   % Spike!
        [ispks,jspks] =  find(rrcum>=repmat(tspnext,nii,1));
        spcells = unique(jspks(ispks == min(ispks))); % cell number(s)
        ispk = iinxt(min(ispks)); % time bin of spike(s)
        rprev = rrcum(min(ispks),:); % grab accumulated history to here

        % Record this spike
        mxi = min(rlen, ispk+hlen); % determine bins for adding h current
        iiPostSpk = ispk+1:mxi;
        for ic = 1:length(spcells)
            icell = spcells(ic);
            nsp(icell) = nsp(icell)+1;
            tsp{icell}(nsp(icell),1) = ispk*dt;
            if ~isempty(iiPostSpk)
                Itot(iiPostSpk,:) = Itot(iiPostSpk,:)+ih(1:mxi-ispk,:,icell);
            end
            rprev(icell) = 0;  % reset this cell's integral
            tspnext(icell) = exprnd(1); % draw RV for next spike in this cell
        end
        jbin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = jbin/(sum(nsp));
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
end

% Remove any extra bins from tsp and compute binned spike train 'sps'
sps = zeros(rlen,ncells);
for jj = 1:ncells
    tsp{jj} = tsp{jj}(1:nsp(jj));
    if ~isempty(tsp{jj})
        sps(round(tsp{jj}/dt),jj) = 1;
    end
end


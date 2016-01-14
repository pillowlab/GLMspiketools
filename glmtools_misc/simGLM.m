function [tsp,Vmem,Ispk] = simGLM(glmprs,Stim)
% [tsp, Vmem,Ispk] = simGLM(glmprs,Stim);
% 
% Compute response of glm to stimulus Stim.
%
% Uses time rescaling instead of Bernouli approximation to conditionally
% Poisson process
%
% Dynamics:  Filters the Stimulus with glmprs.k, passes this through a
% nonlinearity to obtain the point-process conditional intensity.  Add a
% post-spike current to the linear input after every spike.

% Check nonlinearity (default is exponential)
if ~isfield(glmprs, 'nlfun'), 
    glmprs.nlfun = @exp;
end

upSampFactor = glmprs.dtStim/glmprs.dtSp; % number of spike bins per Stim bin
assert(mod(upSampFactor,1) == 0, 'dtStim / dtSp must be an integer');

% Run "simple" or "coupled" simulation code
if size(glmprs.k,3) > 1  % Run "coupled" GLM model if multiple cells present
    if nargout <= 2
        [tsp,Vmem] = simGLMcpl(glmprs,Stim,upSampFactor);
    else
        [tsp,Vmem,Ispk] = simGLMcpl(glmprs,Stim,upSampFactor);
    end

else   % Single cell simulation
    if nargout <= 2
        [tsp,Vmem] = simGLMsingle(glmprs,Stim,upSampFactor);
    else
        [tsp,Vmem,Ispk] = simGLMsingle(glmprs,Stim,upSampFactor);
    end
end


% ======================================================================
function [tsp,Vmem,Ispk] = simGLMsingle(glmprs,Stim,upSampFactor)
% Sub-function within simGLM.m
% Simulates the GLM point process model for a single (uncoupled) neuron
% using time-rescaling
    
% --------------- Check Inputs ----------------------------------
nbinsPerEval = 100;  % Default number of bins to update for each spike
dt = glmprs.dtSp; % bin size for simulation

slen = size(Stim,1); % length of stimulus
rlen = slen*upSampFactor;  % length of binned spike response

% -------------  Compute filtered resp to signal ----------------
Istm = sameconv(Stim,glmprs.k);  % filter stimulus with k
Vmem = kron(Istm,ones(upSampFactor,1)) + glmprs.dc; % upsample stim current

% -------------- Compute interpolated h current ----------------------
ih = glmprs.ih;
hlen = length(ih);
    
% --------------- Set up simulation dynamics variables ---------------
nsp = 0; % number of spikes
tsp = zeros(round(slen/4),1);  % pre-allocate space for spike times
if (nargout > 2), 
    Ispk = Vmem*0;
end
jbin = 1; % current time bin
tspnext = exprnd(1);  % time of next spike (in rescaled time)
rprev = 0;  % Integrated rescaled time up to current point

% --------------- Run dynamics ---------------------------------------
while jbin <= rlen
    iinxt = jbin:min(jbin+nbinsPerEval-1,rlen);
    rrnxt = glmprs.nlfun(Vmem(iinxt))*dt; % Cond Intensity
    rrcum = cumsum(rrnxt)+rprev; % integrated cond intensity
    if (tspnext >= rrcum(end)) % No spike in this window
        jbin = iinxt(end)+1;
        rprev = rrcum(end);
    else   % Spike!
        ispk = iinxt(find(rrcum>=tspnext, 1, 'first')); % time bin where spike occurred
        nsp = nsp+1;
        tsp(nsp,1) = ispk*dt; % spike time
        mxi = min(rlen, ispk+hlen); % max time affected by post-spike kernel
        iiPostSpk = ispk+1:mxi; % time bins affected by post-spike kernel
        if ~isempty(iiPostSpk)
            Vmem(iiPostSpk) = Vmem(iiPostSpk)+ih(1:mxi-ispk);
            if nargout == 3  % Record post-spike current
                Ispk(iiPostSpk) = Ispk(iiPostSpk)+ih(1:mxi-ispk);
            end
        end
        tspnext = exprnd(1);  % draw next spike time
        rprev = 0; % reset integrated intensity
        jbin = ispk+1;  % Move to next bin
        % --  Update # of samples per iter ---
        muISI = jbin/nsp;
        nbinsPerEval = max(20, round(1.5*muISI)); 
    end
end
tsp = tsp(1:nsp); % prune extra entries, if necessary


% ========================================================================
function [tsp,Vmem,Ispk] = simGLMcpl(glmprs,Stim,upSampFactor)
% [tsp,Vmem,Ispk] = simGLMcpl(glmprs,Stim,upSampFactor)
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
dt = glmprs.dt;
% Check that dt evenly divides 1
if mod(1,dt) ~= 0
    dt = 1/round(1/dt);
    fprintf(1, 'glmrunmod: reset dtsim = %.3f\n', dt);
end
slen = size(Stim,1); % length of stimulus
rlen = slen*upSampFactor;  % length of binned response

% -------------  Compute filtered resp to signal ----------------
Istm = zeros(slen,ncells);
for jj = 1:ncells
    Istm(:,jj) = sameconv(Stim,glmprs.k(:,:,jj));  % filter stimulus with k
end
% Compute finely sampled version
Vmem = kron(Istm,ones(upSampFactor,1)) + glmprs.dc; % upsample stim current

% -------------- Compute interpolated h current ----------------------
ih = reshape(interp_spikeFilt(glmprs.ih, glmprs.iht, dt),[],ncells,ncells);
ih = permute(ih,[1 3 2]);
hlen = length(ih);

% -------------  Static nonlinearity & spiking -------------------
if nargout > 2
    Ispk = Vmem*0;
end

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
    rrnxt = glmprs.nlfun(Vmem(iinxt,:))*dt/RefreshRate; % Cond Intensity
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
                Vmem(iiPostSpk,:) = Vmem(iiPostSpk,:)+ih(1:mxi-ispk,:,icell);
                if nargout == 3  % Record post-spike current
                    Ispk(iiPostSpk,:)=Ispk(iiPostSpk,:)+ih(1:mxi-ispk,:,icell);
                end
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

% Remove any extra bins from cell array of spike times
for j = 1:ncells
    tsp{j} = tsp{j}(1:nsp(j));
end


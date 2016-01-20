function [logli, dL, H] = Loss_GLM_logli_exp(prs,Xstruct)
% [neglogli, dL, H] = Loss_GLM_logli_exp(prs)
%
% Compute negative log-likelihood of data undr the GLM model with
% exponential nonlinearity (with standard linear parametrization of stim filter) 
%
% Inputs:
%   prs = [kprs - weights for stimulus kernel
%          dc   - dc current injection
%          ihprs - weights on post-spike current];
% Outputs:
%      logli = negative log likelihood of spike train
%      dL = gradient with respect to prs
%      H = hessian

% Extract some vals from Xstruct (Opt Prs);
nktot = Xstruct.nkx*Xstruct.nkt;   % total # params for k
dt = Xstruct.dtSp;           % absolute bin size for spike train (in sec)

% Unpack GLM prs;
kprs = prs(1:nktot);
dc = prs(nktot+1);
ihprs = prs(nktot+2:end);

% Extract some other stuff we'll use a lot
XXstm = Xstruct.Xstim; % stimulus design matrix
XXsp = Xstruct.Xsp;    % spike history design matrix
bsps = Xstruct.bsps;   % binary spike vector
M = Xstruct.Minterp;   % matrix for interpolating from stimulus bins to spike train bins
ihflag = Xstruct.ihflag;  % flag
rlen = Xstruct.rlen;   % number of bins in spike train vector
nsp = sum(bsps);     % number of spikes

% -------- Compute sum of filter reponses -----------------------
if Xstruct.ihflag
    Itot = M*(XXstm*kprs) + XXsp*ihprs + dc; % stim-dependent + spikehist-dependent inputs
else
    Itot = M*(XXstm*kprs) + dc;  % stim-dependent input only
end

% ---------  Compute output of nonlinearity  ------------------------
rr = exp(Itot);

% ---------  Compute log-likelihood ---------------------------------
Trm1 = sum(rr)*dt;  % non-spike term
Trm2 = -sum(Itot(bsps)); % spike term
logli = Trm1 + Trm2;

% ---------  Compute Gradient -----------------
if (nargout > 1)
    
    % Non-spiking terms (Term 1)
    dLdk0 = (rr'*M*XXstm)';
    dLdb0 = sum(rr);
    if ihflag, dLdh0 = (rr'*XXsp)'; 
    end
    
    % Spiking terms (Term 2)
    Msp = M(bsps,:); % interpolation matrix just for spike bins
    dLdk1 = (sum(Msp*XXstm))';
    dLdb1 = nsp;
    if ihflag, dLdh1 = sum(XXsp(bsps,:),1)';
    end

    % Combine terms
    dLdk = dLdk0*dt - dLdk1;
    dLdb = dLdb0*dt - dLdb1;
    if ihflag, dLdh = dLdh0*dt - dLdh1;
    else dLdh = [];
    end
    dL = [dLdk; dLdb; dLdh];
    
end

% ---------  Compute Hessian -----------------
if nargout > 2
    % --- Non-spiking terms -----
    
    % multiply each row of M with drr
    ddrrdiag = spdiags(rr,0,rlen,rlen); 
    ddqqIntrp = ddrrdiag*M; % this is MUCH faster than using bsxfun, due to sparsity!

    % k and b terms
    Hk = (XXstm'*(M'*ddqqIntrp)*XXstm)*dt; % Hkk (k filter)
    Hb = dLdb0*dt;                       % Hbb (constant b)
    Hkb = (sum(ddqqIntrp,1)*XXstm)'*dt;    % Hkb (cross-term)
    if ihflag  % h terms
        Hh = (XXsp'*(bsxfun(@times,XXsp,rr)))*dt;  % Hh (h filter)
        % (here bsxfun is faster than diagonal multiplication)
        
        Hkh = ((XXsp'*ddqqIntrp)*XXstm*dt)';         % Hhk (cross-term)
        Hhb = (rr'*XXsp)'*dt;                      % Hhb (cross-term)
    else
        Hkh=[]; Hhb=[]; Hh=[];
    end

    H = [[Hk Hkb Hkh]; [Hkb' Hb Hhb']; [Hkh' Hhb Hh]];
end


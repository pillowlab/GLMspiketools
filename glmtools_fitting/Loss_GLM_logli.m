function [logli, dL, H] = Loss_GLM_logli(prs,Xstruct)
% [logli, dL, H] = Loss_GLM_logli(prs,Xstruct)
%
% Compute negative log-likelihood of data undr the GLM model
% (with standard linear parametrization of stimulus kernel);
%
% Uses arbitrary nonlinearity 'nlfun' instead of exponential
%
% Inputs:
%   prs = [kprs - weights for stimulus kernel
%          dc   - dc current injection
%          ihprs - weights on post-spike current]
%
% Outputs:
%      logli = negative log likelihood of spike train
%      dL = gradient with respect to prs
%      H = hessian


% Unpack GLM prs;
nktot = Xstruct.nkx*Xstruct.nkt;   % total # params for k
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
nsp = (sum(bsps));     % number of spikes
dt = Xstruct.dtSp;  % absolute bin size for spike train (in sec)

% -------- Compute sum of filter reponses -----------------------
if Xstruct.ihflag
    Itot = M*(XXstm*kprs) + XXsp*ihprs + dc; % stim-dependent + spikehist-dependent inputs
else
    Itot = M*(XXstm*kprs) + dc;  % stim-dependent input only
end

% ---------  Compute output of nonlinearity  ------------------------
[rr,drr,ddrr] = Xstruct.nlfun(Itot);

% ---------  Compute log-likelihood ---------------------------------
Trm1 = sum(rr)*dt;  % non-spike term
Trm2 = -sum(log(rr(bsps))); % spike term
logli = Trm1 + Trm2;

% ---------  Compute Gradient -----------------
if (nargout > 1)
    
    % Non-spiking terms (Term 1)
    dLdk0 = (drr'*M*XXstm)';
    dLdb0 = sum(drr);
    if ihflag, dLdh0 = (drr'*XXsp)'; 
    end
    
    % Spiking terms (Term 2)
    Msp = M(bsps,:); % interpolation matrix just for spike bins
    frac1 = drr(bsps)./rr(bsps);
    dLdk1 = ((frac1'*Msp)*XXstm)';
    dLdb1 = sum(frac1);
    if ihflag, dLdh1 = (frac1'*XXsp(bsps,:))';
    end

    % Combine Term 1 and Term 2
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
    ddrrdiag = spdiags(ddrr,0,rlen,rlen); 
    ddqqIntrp = ddrrdiag*M; % this is MUCH faster than using bsxfun, due to sparsity!

    % k and b terms
    Hk = (XXstm'*(M'*ddqqIntrp)*XXstm)*dt; % Hkk (k filter)
    Hb = sum(ddrr)*dt;                   % Hbb (constant b)
    Hkb = (sum(ddqqIntrp,1)*XXstm)'*dt;  % Hkb (cross-term)
    if ihflag  % h terms
        Hh = (XXsp'*(bsxfun(@times,XXsp,ddrr)))*dt;  % Hh (h filter)
        % (here bsxfun is faster than diagonal multiplication)
                
        Hkh = ((XXsp'*ddqqIntrp)*XXstm*dt)';         % Hhk (cross-term)
        Hhb = (ddrr'*XXsp)'*dt;                      % Hhb (cross-term)
    else
        Hkh=[]; Hhb=[]; Hh=[];
    end

    % --- Add in spiking terms ----
    frac2 = (rr(bsps).*ddrr(bsps) - drr(bsps).^2)./rr(bsps).^2; % needed weights from derivation of Hessian
    fr2Interp = spdiags(frac2,0,nsp,nsp)*Msp; % rows of Msp re-weighted by these weights

    % Spiking terms, k and b
    Hk= Hk - XXstm'*(Msp'*fr2Interp)*XXstm; % Hkk (k filter)
    Hb =  Hb-sum(frac2);           % Hbb (constant b)
    Hb1 = sum(fr2Interp,1)*XXstm;  % Spiking term, k and const
    Hkb = Hkb - Hb1';
    if Xstruct.ihflag
        XXrrsp = bsxfun(@times,XXsp(bsps,:),frac2);
        Hh= Hh - XXsp(bsps,:)'*XXrrsp;
        % Const by h
        Hhb1 = sum(XXrrsp,1)';
        Hhb = Hhb - Hhb1;
        % k by h term
        Hkh0 = XXsp(bsps,:)'*(fr2Interp*XXstm);
        Hkh = Hkh-Hkh0';
    end
    H = [[Hk Hkb Hkh]; [Hkb' Hb Hhb']; [Hkh' Hhb Hh]];
end


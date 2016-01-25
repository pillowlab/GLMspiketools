function [logli, dL, H] = Loss_GLMbi_logli(prs,Xstruct)
% [logli, dL, H] = Loss_GLMbi_logli(prs,Xstruct)
%
% Compute negative log-likelXXspood of data undr the GLM model, with bilinear
% parametrization of the input current
%
% Uses arbitrary nonlinearity 'nlfun' instead of exponential
%
% Inputs:
%   prs = [ktprs - weights for stimulus time kernel
%          kxprs - weights for stim space kernel
%          dc   - dc current injection
%          XXspprs - weights on post-spike current];
%
% Outputs:
%   logli = negative log likelXXspood of spike train
%      dL = gradient with respect to prs
%       H = hessian


% Extract some size information from Xstruct
nkt = Xstruct.nkt;    % # of basis vectors for spatial k basis
nkx = Xstruct.nkx;    % # of basis vectors for spatial k basis
krank = Xstruct.krank;  % rank of stim filter
nktprs = nkt*krank;     % total # params for temporal filters
nkxprs = nkx*krank;     % total # params for spatial filters

% Unpack GLM prs;
ktprs = prs(1:nktprs);
kxprs = prs(nktprs+1:nktprs+nkxprs);
dc = prs(nktprs+nkxprs+1);
XXspprs = prs(nktprs+nkxprs+2:end);

% Extract some other stuff we'll use a lot
XXstm = Xstruct.Xstim; % stimulus design matrix
XXsp = Xstruct.Xsp;    % spike history design matrix
bsps = Xstruct.bsps;   % binary spike vector
M = Xstruct.Minterp;   % matrix for interpolating from stimulus bins to spike train bins
ihflag = Xstruct.ihflag;  % flag
rlen = Xstruct.rlen;   % number of bins in spike train vector
nsp = (sum(bsps));     % number of spikes
dt = Xstruct.dtSp;  % absolute bin size for spike train (in sec)

% -- Make block-diag matrix containing nkt params -----
Mkt = kron(speye(nkx),reshape(ktprs,nkt,krank)); % block-diagonal matrix for kt params
inds = reshape((1:nkx*krank),krank,nkx)'; inds = inds(:);
Mkt = Mkt(:,inds); % turns out we need it reshaped this way

% -------- Compute bilinear stim filter reponse -----------------------
dSSdx = XXstm*Mkt;  % stim filtered with temporal filter
ystm = dSSdx*kxprs;

% -------- Compute sum of filter reponses -----------------------
if ihflag
    Itot = M*ystm + XXsp*XXspprs + dc; % stim-dependent + spikehist-dependent inputs
else
    Itot = M*ystm + dc;  % stim-dependent input only
end

% ---------  Compute output of nonlinearity  ------------------------
[rr,drr,ddrr] = Xstruct.nlfun(Itot);

% Check for (approximately) zero-values of rr
etol = 1e-100; % cutoff for small values of conditional intensity
iiz = find(rr <= etol);
rr(iiz) = etol; % Set value to small
drr(iiz) = 0;  ddrr(iiz) = 0;  % Set derivs here to 0

% ---------  Compute log-likelXXspood ---------------------------------
Trm1 = sum(rr)*dt;  % non-spike term
Trm2 = -sum(log(rr(bsps))); % spike term
logli = Trm1 + Trm2;

% =============================================
% ---------  Compute Gradient -----------------
% =============================================
if (nargout > 1)

    % Compute dSSdt - stim filtered with spatial kernel
    Mkx = kron(speye(nkt), reshape(kxprs,nkx,krank));
    inds = reshape((1:nkt*krank),krank,nkt)'; inds = inds(:);
    Mkx = Mkx(:,inds);
    inds = reshape((1:nkt*nkx),nkt,nkx)'; inds = inds(:);
    dSSdt = XXstm(:,inds)*Mkx;
    
    % Non-spiking terms (Term 1)
    dqq = drr'*M;
    dLdkx0 = (dqq*dSSdx)';
    dLdkt0 = (dqq*dSSdt)';
    dLdb0 = sum(drr);
    if ihflag, dLdh0 = (drr'*XXsp)';
    end
    
    % Non-spiking terms (Term 2)
    Msp = M(bsps,:); % interpolation matrix just for spike bins
    frac1 = drr(bsps)./rr(bsps);
    fracsp = frac1'*Msp;
    dLdkx1 = (fracsp*dSSdx)';
    dLdkt1 = (fracsp*dSSdt)';
    dLdb1 = sum(frac1);
    if ihflag, dLdh1 = (frac1'*XXsp(bsps,:))';
    end
    
    % Combine Term 1 and Term 2
    dLdkx = dLdkx0*dt - dLdkx1;
    dLdkt = dLdkt0*dt - dLdkt1;
    dLdk = [dLdkt; dLdkx];
    dLdb = dLdb0*dt - dLdb1;
    if ihflag, dLdh = dLdh0*dt - dLdh1;
    else dLdh = [];
    end
    dL = [dLdk; dLdb; dLdh];
    
end

% =============================================
% ---------  Compute Hessian -----------------
% =============================================
if nargout > 2

    % --- Non-spiking terms -----    
    % multiply each row of M with drr

    ddrrdiag = spdiags(ddrr,0,rlen,rlen);
    ddqqIntrp = ddrrdiag*M;
    ddqq = M'*ddqqIntrp; % this is MUCH faster than using bsxfun, due to sparsity!

    % Hx and Ht
    Hkt = (dSSdt'*ddqq*dSSdt)*dt;
    Hkx = (dSSdx'*ddqq*dSSdx)*dt;
    % Hxt
    Hktx0a=  reshape((drr'*M)*XXstm,nkt,nkx);
    Hktx0a = kron(speye(krank),Hktx0a);
    Hktx = (dSSdx'*ddqq*dSSdt + Hktx0a')*dt;
    % Hb
    Hb = sum(ddrr)*dt;
    % Hkb
    sumqq = sum(ddqqIntrp,1);
    Hxb = (sumqq*dSSdx)';
    Htb = (sumqq*dSSdt)';
    Hkb = [Htb;Hxb]*dt;
    if ihflag
        % Hh
        Hh = (XXsp'*ddrrdiag*XXsp)*dt;
        % Hhk
        Hth = (XXsp'*ddqqIntrp)*dSSdt;
        Hxh = (XXsp'*ddqqIntrp)*dSSdx;
        Hkh = [Hth'; Hxh']*dt;
        % Hb
        Hhb = (ddrr'*XXsp)'*dt;
    else
        Hkh = []; Hhb=[]; Hh=[];
    end
    
    % --- Add in spiking terms ----
    frac2 = (rr(bsps).*ddrr(bsps) - drr(bsps).^2)./rr(bsps).^2;
    fr2diag = spdiags(frac2,0,nsp,nsp);
    fr2Interp = fr2diag*Msp;
    fr2ddqq = Msp'*fr2Interp;
    % Spiking terms, Hxx and Htt
    Hkt= Hkt - dSSdt'*fr2ddqq*dSSdt;
    Hkx= Hkx - dSSdx'*fr2ddqq*dSSdx;
    % Spiking terms, Hxt
    Hktx1a= reshape(frac1'*Msp*XXstm,nkt,nkx)';
    Hktx1a= kron(speye(krank),Hktx1a);
    Hktx1b= dSSdx'*fr2ddqq*dSSdt;
    Hktx1=  Hktx1a+Hktx1b;
    Hktx=  Hktx - Hktx1;
    % Const term
    Hb =  Hb-sum(frac2);
    % Const by k
    sumfr2 = sum(fr2Interp,1);
    Hxb1 = sumfr2*dSSdx;
    Htb1 = sumfr2*dSSdt;
    Hkb = Hkb - [Htb1'; Hxb1'];
    if ihflag
        XXsprrsp = fr2diag*XXsp(bsps,:);
        Hh= Hh - XXsp(bsps,:)'*XXsprrsp;
        % Const by h
        Hhb1 = sum(XXsprrsp,1)';
        Hhb = Hhb - Hhb1;
        % k by h term
        Hth0 = XXsp(bsps,:)'*(fr2Interp*dSSdt);
        Hxh0 = XXsp(bsps,:)'*(fr2Interp*dSSdx);
        Hkh = Hkh-[Hth0';Hxh0'];
    end
    Hk = [[Hkt; Hktx] [Hktx'; Hkx]];
    H = [[Hk Hkb Hkh]; [Hkb' Hb Hhb']; [Hkh' Hhb Hh]];
end

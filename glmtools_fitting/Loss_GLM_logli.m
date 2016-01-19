function [logli, dL, H] = Loss_GLM_logli(prs,XX)
% [neglogli, dL, H] = Loss_GLM_logli(prs)
%
% Compute negative log-likelihood of data undr the GLM model
% (with standard linear parametrization of stimulus kernel);
%
% Uses arbitrary nonlinearity 'nlfun' instead of exponential
%
% Inputs:
%   prs = [kprs - weights for stimulus kernel
%          dc   - dc current injection
%          ihprs - weights on post-spike current];
% Outputs:
%      logli = negative log likelihood of spike train
%      dL = gradient with respect to prs
%      H = hessian


etol = 1e-100;  % fudge factor for avoiding log(0)

% Extract some vals from Xstruct (Opt Prs);
nktot = Xstruct.nkx*Xstruct.nkt;   % total # params for k
nh = Xstruct.nh;             % # basis vectors for history current 
nh2 = Xstruct.nh2;           % # basis vectors for coupling currents
nCoupled = Xstruct.ncoupled; % number of neurons coupled to this one
nhtot = nh+nh2*nCoupled;     % total # of params for history/coupling kernels
dt = Xstruct.dtSp;           % absolute bin size for spike train (in sec)

% Unpack GLM prs;
kprs = prs(1:nktot);
dc = prs(nktot+1);
ihprs = prs(nktot+2:end);
if length(ihprs)~= nhtot
    error('Problem with number of params passed in (Loss_GLM_logli.m)');
end

% -------- Compute stim filter reponse -----------------------
SS = Xstruct.MSTM(jwin(1)+1:jwin(2),:);  % Relevant stim
MM = Xstruct.MMintrp(ii(1):ii(2),jj(1):jj(2)); % Interpolation filter for this chunk
Istm = Xstruct.Xstim*kprs+dc;  % filtered stimulus

% -------- Compute net h filter response -----------------------
if Xstruct.ihflag
    Iinj = kron(Istm,ones(Xstruct.upsampfactor,1)) + Xstruct.Xsp*hprs;
else
    Iinj = kron(Istm,ones(Xstruct.upsampfactor,1));
end

% ---------  Compute output of nonlinearity  ------------------------
[rr,drr,ddrr] = Xstruct.nlfun(Iinj);

% ---------  Compute log-likelihood ---------------------------------
Trm1 = sum(rr)*dt/RefreshRate;  % non-spike term
Trm2 = -sum(log(rr(Xstruct.spInds))); % spike term
logli = Trm1 + Trm2;

% ---------  Compute Gradient -----------------
if (nargout > 1)
    
    % Non-spiking terms
    dLdk0 = ((drr'*MM)*SS)';
    dLdb0 = sum(drr);
    if Xstruct.ihflag
        dLdh0 = (drr'*Ih)';
    end
    if spflag
        MMsp = MM(i_sp,:);
        %isprel = i_sp+ii(1)-1;
        frac1 = drr(i_sp)./rr(i_sp);
        dLdk1 = ((frac1'*MMsp)*SS)';
        dLdb1 = sum(frac1);
        if Xstruct.ihflag
            dLdh1 = (frac1'*Ih(i_sp,:))';
        end
    else
        dLdk1=0; dLdb1=0; dLdh1=0;
    end
    dLdk = dLdk0*dt/RefreshRate - dLdk1;
    dLdb = dLdb0*dt/RefreshRate - dLdb1;
    if Xstruct.ihflag
        dLdh = dLdh0*dt/RefreshRate - dLdh1;
    else
        dLdh = [];
    end
    dL = dL + [dLdk; dLdb; dLdh];
    
end

% ---------  Compute Hessian -----------------
if nargout > 2
    ddrrdiag = spdiags(ddrr,0,ilen,ilen); % make diagonal matrix
    ddqqIntrp = ddrrdiag*MM;  % multiply each row of MM by ddrr (faster than bsxfun!)
    ddqq = MM'*ddqqIntrp;
    %Hxx and Htt
    Hk = (SS'*ddqq*SS)*dt/RefreshRate;
    % Hb
    Hb = sum(ddrr)*dt/RefreshRate;
    % Hkb
    sumqq = sum(ddqqIntrp,1);
    Hkb = (sumqq*SS)'*dt/RefreshRate;
    if Xstruct.ihflag
        % Hh
        Hh = (Ih'*(ddrrdiag*Ih))*dt/RefreshRate;
        % Hhk
        Hkh = ((Ih'*ddqqIntrp)*SS*dt/RefreshRate)';
        % Hb
        Hhb = (ddrr'*Ih)'*dt/RefreshRate;
    else
        Hkh = []; Hhb=[]; Hh=[];
    end
    if spflag  % -------  If spikes in window ------
        frac2 = (rr(i_sp).*ddrr(i_sp) - drr(i_sp).^2)./rr(i_sp).^2;
        fr2diag = spdiags(frac2,0,length(i_sp),length(i_sp));
        fr2Interp = fr2diag*MMsp;
        fr2ddqq = MMsp'*fr2Interp;
        % Spiking terms, Hk
        Hk= Hk - SS'*fr2ddqq*SS;
        % Spiking term, const
        Hb =  Hb-sum(frac2);
        % Spiking term, k and const
        sumfr2 = sum(fr2Interp,1);
        Hb1 = sumfr2*SS;
        Hkb = Hkb - Hb1';
        if Xstruct.ihflag
            Ihrrsp = fr2diag*Ih(i_sp,:);
            Hh= Hh - Ih(i_sp,:)'*Ihrrsp;
            % Const by h
            Hhb1 = sum(Ihrrsp,1)';
            Hhb = Hhb - Hhb1;
            % k by h term
            Hkh0 = Ih(i_sp,:)'*(fr2Interp*SS);
            Hkh = Hkh-Hkh0';
        end
    end
    H = H + [[Hk Hkb Hkh]; [Hkb' Hb Hhb']; [Hkh' Hhb Hh]];
end


function [logli, dL, H] = Loss_GLMbi_logli(prs)
% [neglogli, dL, H] = Loss_GLMbi_logli(prs)
%
% Compute negative log-likelihood of data undr the GLM model, with bilinear
% parametrization of the input current
%
% Uses arbitrary nonlinearity 'nlfun' instead of exponential
%
% Inputs:
%   prs = [ktprs - weights for stimulus time kernel
%          kxprs - weights for stim space kernel
%          dc   - dc current injection
%          ihprs - weights on post-spike current];
% Outputs:
%      logli = negative log likelihood of spike train
%      dL = gradient with respect to prs
%      H = hessian


etol = 1e-100;

% global vars for optimization
global MSTM MMntrp SPNDS SPNDS2 OPRS RefreshRate 

% Extract some vals from OPRS (Opt Prs);
nprs = length(prs);
nt = OPRS.nt;nx = OPRS.nx; 
krank = OPRS.krank;
nh = OPRS.nh; nh2 = OPRS.nh2;
nsp = length(SPNDS);
nCoupled = length(SPNDS2);
nhtot = nh+nh2*nCoupled;
dt = OPRS.dt; ndt = round(1/dt);
slen = size(MSTM,1); rlen = slen*ndt;

% Unpack GLM prs;
ktprs = prs(1:OPRS.nkt);
kxprs = prs(OPRS.nkt+1:OPRS.nkt+OPRS.nkx);
dc = prs(OPRS.nkt+OPRS.nkx+1);
ihprs = prs(OPRS.nkt+OPRS.nkx+2:end);

% -- Make block-diag matrix containing kt params -----
ktprs2 = repcell(sparse(reshape(ktprs,nt,krank)), nx);
Mkt = blkdiag(ktprs2{:});
inds = reshape([1:nx*krank],krank,nx)'; inds = inds(:);
Mkt = Mkt(:,inds);

% Initialize likelihood, gradient and Hessian -----------
logli = 0;
dL = zeros(nprs,1);
H = zeros(nprs,nprs);

for jch = 1:OPRS.nchunks
    %------- Compute indices relevant to this chunk ---------------
    iwin = OPRS.ichunk([jch, jch+1]); % abs fine bin win
    iwin1 = [iwin(1)+1,iwin(2)];
    jwin = OPRS.jchunk(jch,:);   % abs window of stimulus, coarse
    jj = [1,diff(jwin)];
    ilen = diff(iwin);  
    jlen = diff(jwin);
    ii = [(iwin(1)-jwin(1)*ndt+1),(iwin(1)-jwin(1)*ndt+ilen)];
    %   if (iwin(1)/ndt-1)  == jwin(1)  %--- evenly-dividing bins ----
    %   ii = [ndt+1, ndt+ilen];
    
    % -------- Compute stim filter reponse -----------------------
    SS = MSTM(jwin(1)+1:jwin(2),:);  % Relevant stim
    dSSdx = SS*Mkt;  % stim filtered with temporal filter
    MM = MMntrp(ii(1):ii(2),jj(1):jj(2));
    ystm = dSSdx*kxprs;
    ystmhi = MM*ystm + dc;

    % -------------- Compute net h current -----------------------
    if OPRS.ihflag
        Ih = zeros(ilen,nh+nh2*nCoupled);
        Ih(:,1:nh) = spikeconv_mex(SPNDS, OPRS.ihbas, iwin1);
        for jcpl = 1:nCoupled
            Ih(:,nh+nh2*(jcpl-1)+1:nh+nh2*jcpl) = ...
                spikeconv_mex(SPNDS2{jcpl}, OPRS.ihbas2, iwin1);
        end
        yh = Ih*ihprs;
        Iinj = ystmhi+yh;
    else
        Iinj = ystmhi;
    end
    
    % ---------  Compute likelihood itself  ------------------------
    i_sp = in(SPNDS, iwin+.5)-iwin(1);
    [rr,drr,ddrr] = OPRS.nlfun(Iinj);
    spflag = ~isempty(i_sp);

    % Check for zero-values of rr
    iiz = find(rr <= etol);
    rr(iiz) = etol; % Set value to small
    drr(iiz) = 0;  ddrr(iiz) = 0;  % Set derivs here to 0

    Trm1 = sum(rr)*dt/RefreshRate;  % non-spike term
    if spflag, 
        Trm2 = -sum(log(rr(i_sp))); % spike term
    else, Trm2 = 0;
    end    
    logli = logli + (Trm1 + Trm2);

    % ---------  Compute Gradient -----------------
    if (nargout > 1)
        % Compute dSSdt - stim filtered with spatial kernel 
        Mkx = repcell(reshape(kxprs,nx,krank), nt);
        Mkx = blkdiag(Mkx{:});
        inds = reshape([1:nt*krank],krank,nt)'; inds = inds(:);
        Mkx = Mkx(:,inds);
        inds = reshape([1:nt*nx],nt,nx)'; inds = inds(:);
        dSSdt = SS(:,inds)*Mkx;

        % Non-spiking terms
        dqq = drr'*MM;
        dLdkx0 = (dqq*dSSdx)'; 
        dLdkt0 = (dqq*dSSdt)';
        dLdb0 = sum(drr);
        if OPRS.ihflag
            dLdh0 = (drr'*Ih)';
        end
        if spflag
            MMsp = MM(i_sp,:);
            %isprel = i_sp+ii(1)-1;
            frac1 = drr(i_sp)./rr(i_sp);
            fracsp = frac1'*MMsp;
            dLdkx1 = (fracsp*dSSdx)';
            dLdkt1 = (fracsp*dSSdt)';
            dLdb1 = sum(frac1);
            if OPRS.ihflag
                dLdh1 = (frac1'*Ih(i_sp,:))';
            end
        else
            dLdkx1=0; dLdkt1=0; dLdb1=0; dLdh1=0;
        end
        dLdkx = dLdkx0*dt/RefreshRate - dLdkx1;
        dLdkt = dLdkt0*dt/RefreshRate - dLdkt1;
        dLdk = [dLdkt; dLdkx];
        dLdb = dLdb0*dt/RefreshRate - dLdb1;
        if OPRS.ihflag
            dLdh = dLdh0*dt/RefreshRate - dLdh1;
        else
            dLdh = [];
        end
        dL = dL + [dLdk; dLdb; dLdh];

    end

    % ---------  Compute Hessian -----------------
    if nargout > 2
        ddrrdiag = spdiags(ddrr,0,ilen,ilen);
        ddqqIntrp = ddrrdiag*MM;
        ddqq = MM'*ddqqIntrp;
        %Hxx and Htt
        Hkt = (dSSdt'*ddqq*dSSdt)*dt/RefreshRate;
        Hkx = (dSSdx'*ddqq*dSSdx)*dt/RefreshRate;
        %Hxt
        Hktx0a=  reshape((drr'*MM)*SS,nt,nx);
        Hktx0a = repcell(Hktx0a,krank);
        Hktx0a = blkdiag(Hktx0a{:});
        Hktx = (dSSdx'*ddqq*dSSdt + Hktx0a')*dt/RefreshRate;
        % Hb
        Hb = sum(ddrr)*dt/RefreshRate;
        % Hkb
        sumqq = sum(ddqqIntrp,1);
        Hxb = (sumqq*dSSdx)';
        Htb = (sumqq*dSSdt)';
        Hkb = [Htb;Hxb]*dt/RefreshRate;
        if OPRS.ihflag
            % Hh
            Hh = (Ih'*ddrrdiag*Ih)*dt/RefreshRate;
            % Hhk
            Hth = (Ih'*ddqqIntrp)*dSSdt;
            Hxh = (Ih'*ddqqIntrp)*dSSdx;
            Hkh = [Hth'; Hxh']*dt/RefreshRate;
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
            % Spiking terms, Hxx and Htt
            Hkt= Hkt - dSSdt'*fr2ddqq*dSSdt;
            Hkx= Hkx - dSSdx'*fr2ddqq*dSSdx;
            % Spiking terms, Hxt
            Hktx1a= reshape(frac1'*MMsp*SS,nt,nx);
            Hktx1a= repcell(Hktx1a,krank);
            Hktx1a= blkdiag(Hktx1a{:})';
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
            if OPRS.ihflag
                Ihrrsp = fr2diag*Ih(i_sp,:);
                Hh= Hh - Ih(i_sp,:)'*Ihrrsp;
                % Const by h
                Hhb1 = sum(Ihrrsp,1)';
                Hhb = Hhb - Hhb1;
                % k by h term
                Hth0 = Ih(i_sp,:)'*(fr2Interp*dSSdt);
                Hxh0 = Ih(i_sp,:)'*(fr2Interp*dSSdx);
                Hkh = Hkh-[Hth0';Hxh0'];
            end
        end
        Hk = [[Hkt; Hktx] [Hktx'; Hkx]];
        H = H + [[Hk Hkb Hkh]; [Hkb' Hb Hhb']; [Hkh' Hhb Hh]];
    end
    
end

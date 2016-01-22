
% ---------------
% Make param object with "true" params (if desired).
% ---------------
Filter_rank = 1;
ggTrue = makeFittingStruct_GLMbi(ggsim.k,dtSp,Filter_rank,ggsim);
% Set true kernel in this struct (i.e. not represented by default basis).
[u,s,v] = svd(ggsim.k);  
ggTrue.k = ggsim.k;
ggTrue.ktbas = eye(nkt);
ggTrue.kt = u(:,1);  
ggTrue.kx = v(:,1)*s(1,1);
% Insert spike times
ggTrue.tsp = tsp;

% ---------------
% Check that conditional intensity calc is correct 
% (if desired, compare to vmem returned by simGLM above)
[logliTrue, rrTrue,tt] = neglogli_GLM(ggTrue,Stim);
subplot(211); plot(tt,vmem,tt,log(rrTrue));
title('total filter output (computed 2 ways)');
subplot(212); plot(tt,log(rrTrue)-vmem);
title('difference');

% ---------------
[logliTrue, rrT,tt,Iinj,Istm,Ih] = neglogli_GLM(ggTrue,Stim);

%% 5. Unit tests

% On total log-conditional intensity
assert(max(abs(vmem-Iinj))<1e-8,'Unit failed: log-conditional intensity consistency');

% On just the spike-history component of the conditional intensity
assert(max(abs(ispk-Ih))<1e-8,'Unit failed: spike-history filter output consistency');

fprintf('Unit passed: condIntensityConsistency\n');



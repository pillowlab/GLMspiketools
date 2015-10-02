function [nl_est,xx,nl_std,II] = compEmpiricalNlin_GLM(gg,Stim,nbins);
% [nl_est,xx,nl_std,II] = compEmpiricalNlin_GLM(gg,Stim,nbins);
% 
% Computes a histogram estimate of the nonlinearity for a GLM by regressing
% (binned) filter outputs against the probability of spiking
%
% Inputs:  gg = glm param object (which includes spike times)
%          Stim = stimulus
%          bins = number of bins to use, or bin centers (optional)
%          
% Outputs: nl_est = y values of nonlinearity;  
%          xx = x values of nonlinearity;
%          nl_std = stdev of nonlinearity estimate at each xx 
%          II = net linear input vs. time (before nonlinearity)


global RefreshRate SPNDS;

if nargin < 3
    nbins = 200;
end

% Compute conditional filter output
[logli,CondInstensity,tt,II] = neglogli_GLM(gg,Stim);
II = sum(II,2);
slen = length(II);
IIsp = II(SPNDS);
II = sort(II);

% Initialize variables
spdist = zeros(nbins,1);
rawdist = zeros(nbins,1);
xx = zeros(nbins,1);

% Set up bins using quantiles
dbin = slen/nbins;
binedges = round([0:dbin:slen]);

% Compute spiking distributions
spdist = histc(IIsp, [II(binedges(1:end-1)+1);inf]);
spdist = spdist(1:nbins);

% Compute raw distribution and xx bin centers
for jj = 1:nbins
    rawdist(jj) = binedges(jj+1)-binedges(jj);
    xx(jj) = mean(II([binedges(jj)+1:binedges(jj+1)]));
end

% Compute nonlinearity and stdev
psp = spdist./rawdist;
nl_est = psp./gg.dt*RefreshRate;
nl_std = sqrt(psp.*(1-psp))./sqrt(rawdist)/gg.dt*RefreshRate;


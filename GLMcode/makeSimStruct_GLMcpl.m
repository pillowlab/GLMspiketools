function S = makeSimStruct_GLMcpl(varargin);
%  S = makeSimStruct_GLMcpl(gg1, gg2, gg3, ...);
%
%  Creates a structure for simulating coupled GLM from multiple
%  single-neuron GLMs.
%
%  Input: structures params of glms to combine

ncells = nargin;

gg = varargin{1};
[klen,kwid] = size(gg.k);
kk = zeros(klen,kwid,ncells);  % Stimulus filters
nh = length(gg.iht);
ih = zeros(nh,ncells,ncells);  % post-spike & coupling kernels

ihflag = 0;  % Check for spike-history kernel
if ~isempty(gg.ih), ihflag = 1; end

ih2flag = 0; % Check for presence of coupling kernels
if size(gg.ih,2)>1
    ih2flag = 1; end
if isfield(gg, 'ih2'); if ~isempty(gg.ih2);
        ih2flag = 1; end; end

ihbasflag = 0;  % Basis for history kernel?
if isfield(gg, 'ihbas'); if ~isempty(gg.ihbas)
        ihbasflag = 1; end; end

ihbas2flag=0; % Distinct basis for coupling kernels?
if isfield(gg, 'ihbas2'); if ~isempty(gg.ihbas2)
        ihbas2flag = 1; end; end

% Add cells to single param structure 
nltestvals = [0, 1, 2];
for j = 1:ncells
  gg = varargin{j};
  kk(:,:,j) = gg.k;
  dc(j) = gg.dc;
  if isfield(gg, 'nlfun')
      nlfuns{j} = gg.nlfun;
  else
      nlfuns{j} = @exp;
  end
  NLtst(j,:) = nlfuns{j}(nltestvals);

  if ihflag  % Post-spike kernel
      if ihbasflag
          ih(:,j,j) = gg.ihbas*gg.ih(:,1);
      else
          ih(:,j,j) = gg.ih(:,1);
      end
  end
 
  if ih2flag % Coupling kernels
      if ihbas2flag
          ih(:,setdiff(1:ncells,j),j) = gg.ihbas2*gg.ih2;
      elseif ihbasflag
          ih(:,setdiff(1:ncells,j),j) = gg.ihbas*gg.ih(:,2:end);
      else
          ih(:,setdiff(1:ncells,j),j) = gg.ih(:,2:end);
      end
  end
  
  
end

S.k = kk;
S.dc = dc;
S.ih = ih;
S.iht = varargin{1}.iht;
S.dt = varargin{1}.dt;

if all(min(NLtst)==max(NLtst))
    S.nlfun = nlfuns{1};
else
    S.nlfun = @(x)multiNLfun(x,nlfuns);
end

if isfield(gg,'ihbasprs')
    S.ihbasprs = gg.ihbasprs;
end
if isfield(gg,'ihbasprs2')
    S.ihbasprs2 = gg.ihbasprs2;
end



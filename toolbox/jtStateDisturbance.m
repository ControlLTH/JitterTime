function N = jtStateDisturbance(N,sysid,m,P)
% N = jtStateDisturbance(N,sysid,m,P)
%
% Add a deterministic and/or random disturbance to the state of a
% continuous- or discrete-time state-space system. Can be called at any time
% after jtCalcDynamics.
%
% Arguments:
% N        The JitterTime system.
% sysid    The system ID.
%
% Optional arguments (default to zero if empty):
% m        The mean value of the disturbance. Must be of size n x 1 for
%          continuous systems and size (n + p) x 1 for discrete systems,
%          where n is the number of states and p is the number of outputs.
% P        The covariance of the disturbance. Must be of size n x n for
%          continuous systems and size (n + p) x (n + p) for discrete systems,
%          where n is the number of states and p is the number of outputs.

if ~isfield(N,'idtoindex')
  error('jtStateDisturbance cannot be called before jtCalcDynamics')
end
  
k = N.idtoindex(sysid);

% Find state indexes 
if N.sys{k}.type == 1     % Continuous
  indexes = N.sys{k}.stateindex;
elseif N.sys{k}.type == 2 % Discrete
  indexes = [N.sys{k}.stateindex N.sys{k}.outputindex];
end
n = length(indexes);

if ~isempty(m)
  if ~isequal(size(m), [n 1])
    error('Disturbance mean value must have size %d x 1',n)
  end
  N.m(indexes) = N.m(indexes) + m;
end

if ~(nargin < 4 || isempty(P))
  if ~isequal(size(P), [n n])
    error('Disturbance covariance must have size %d x %d',n,n)
  end
  if ~issymmetric(P) || min(eig(P)) < -eps
    error('Disturbance covariance must be positive semidefinite')
  end
  N.P(indexes,indexes) =  N.P(indexes,indexes) + P;
end

if N.periodicAnalysis
  if ~isempty(m)
    N.dtot(indexes) = N.dtot(indexes) + m;
  end
  if ~(nargin < 4 || isempty(P))
    N.Rtot(indexes,indexes) =  N.Rtot(indexes,indexes) + P;
  end
end

function N = jtBeginPeriodicAnalysis(N)
% N = jtBeginPeriodicAnalysis(N)
%
% Start a periodic covariance analysis. Must be called after jtCalcDynamics.
%
% Return values:
% N         The JitterTime system which must be passed to all other functions.

if ~isfield(N,'idtoindex')
  error('jtBeginPeriodicAnalysis cannot be called before jtCalcDynamics')
end

N.periodicAnalysis = true;
N.Atot = eye(size(N.Ac));
N.Rtot = zeros(size(N.Ac));
N.dtot = zeros(size(N.Ac,1),1);

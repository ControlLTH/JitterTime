function N = jtEndPeriodicAnalysis(N)
% N = jtEndPeriodicAnalysis(N)
%
% Start a periodic covariance analysis. Must be called after jtCalcDynamics.
%
% Return values:
% N         The JitterTime system which must be passed to all other functions.

if ~N.periodicAnalysis
  error('Periodic analysis not started')
end

if max(abs(eig(N.Atot))) > 1-sqrt(eps)
  error('Periodic analysis failed due to unstable system mode(s)')
end

N.Pperiodic = dlyap(N.Atot,N.Rtot);             % solve P = Atot*P*Atot' + Rtot
N.mperiodic = (eye(size(N.Ac))-N.Atot)\N.dtot;  % solve m = Atot*m + dtot
N.periodicAnalysis = false;

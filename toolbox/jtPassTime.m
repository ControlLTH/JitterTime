function N = jtPassTime(N,T)
% N = jtPassTime(N,T)
%
% Let time pass by T units, running all continuous-time systems and
% accumulating cost.
%
% Arguments:
% N      The JitterTime system.
% T      The time interval.

if ~isfield(N,'idtoindex')
  error('jtPassTime cannot be called before jtCalcDynamics')
end

[Ad,Rd,Qd,Qconst] = calcc2d(N.Ac,N.Rc,N.Qc,T);

costm = N.m' * Qd * N.m;               % Cost due to mean value 
costP = trace(Qd * N.P) + Qconst;      % Cost due to variance
N.J = N.J + costm + costP;
N.dJ = (costm + costP) / T;            % Average cost increase during interval

N.m = Ad * N.m;                        % Update mean value
N.P = Ad * N.P * Ad' + Rd;             % Update variance
if ~isreal(N.P)
  error('Covariance became imaginary')
end
N.Tsim = N.Tsim + T;

if N.periodicAnalysis
  N.Atot = Ad * N.Atot;
  N.Rtot = Ad * N.Rtot * Ad' + Rd;
  N.dtot = Ad * N.dtot;
end


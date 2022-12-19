function N = jtResetSys(N,sysid,ver)
% N = jtResetSys(N,sysid,ver)
%
% Reset the continuous-time system k. The reset dynamics must have been
% previously specified using jtAddResetDynamics. The system state
% is momentarily updated as x(t+) = Phi * x(t) [+ Gamma u(t)].
%
% Arguments:
% N        The JitterTime system.
% sysid    ID of the contiuous-time system to be reset. 
% ver      (optional) What reset dynamics to use (default = 1)

if ~isfield(N,'idtoindex')
  error('jtResetSys cannot be called before jtCalcDynamics')
end
  
k = N.idtoindex(sysid);

if ~(N.sys{k}.type == 1 && N.sys{k}.resetDynamics)
	error('No reset dynamics have been specified for the system')
end

if (nargin < 3) || isempty(ver)
	ver = 1;
end

Ad = N.sys{k}.Ad{ver};

N.m = Ad * N.m;
N.P = Ad * N.P * Ad';    % Rd = 0 for reset updates

if N.periodicAnalysis
  N.Atot = Ad * N.Atot;
  N.Rtot = Ad * N.Rtot * Ad';
end

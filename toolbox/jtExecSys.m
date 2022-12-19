function N = jtExecSys(N,sysid,ver)
% N = jtExecSys(N,sysid,ver)
%
% Execute the discrete-time system k
%
% Arguments:
% N        The JitterTime system.
% sysid    ID of the discrete-time system to be updated.
% ver      (optional) What version to execute (default = 1)

if ~isfield(N,'idtoindex')
  error('jtExecSys cannot be called before jtCalcDynamics')
end

k = N.idtoindex(sysid);

if N.sys{k}.type ~= 2
	error('Can only execute discrete-time systems')
end

if (nargin < 3) || isempty(ver)
	ver = 1;
end

Ad = N.sys{k}.Ad{ver};
Rd = N.sys{k}.Rd{ver};

N.m = Ad * N.m;
N.P = Ad * N.P * Ad' + Rd;
if ~isreal(N.P)
  error('N.P became imaginary')
end

if N.periodicAnalysis
  N.Atot = Ad * N.Atot;
  N.Rtot = Ad * N.Rtot * Ad' + Rd;
end

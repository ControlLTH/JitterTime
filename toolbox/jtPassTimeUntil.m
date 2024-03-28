function N = jtPassTimeUntil(N,time)
% N = jtPassTimeUntil(N,time)
%
% Let time pass until N.Tsim = time, running all continuous-time systems
% and accumulating cost.
%
% Arguments:
% N      The JitterTime system.
% time   The target time. Must be greater than or equal to N.Tsim.

T = time - N.Tsim;
if T < 0
	warning('target time smaller than the current time, ignoring the command')
	return
end

N = jtPassTime(N,T)
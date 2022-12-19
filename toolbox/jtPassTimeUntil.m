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

[A,R,Q,Qconst] = calcc2d(N.Ac,N.Rc,N.Qc,T);
N.J = N.J + trace(Q * N.P) + Qconst;
N.P = A * N.P * A' + R;
N.Tsim = time;

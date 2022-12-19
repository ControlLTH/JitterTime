% Process
s = zpk('s'); z = zpk('z');
Pc = 0.101*(200-s)/((s+0.1)*(s+1.01)*(s+200));
Rc = 1;
Qc = diag([1 1]); % J = E (y^2 + (u+w)^2)

% PID controller
K = 7070; Ti = 0.12; Td = 0.058; N = 100;
Cc = K*(1+1/(Ti*s)+Td*s/(1+s*Td/N));
h = 0.007493883058611;
Cd = -c2d(Cc,h,'foh') * (z-1)/z;  % Controller outputs delta u
Rd = 0.01;

% Integrator 
Hold = ss(0,[],1,[]); 
Rint = 1;
PhiReset = 1;
GammaReset = 1;

N = jtInit;
N = jtAddContSys(N,1,Pc,3,Rc,Qc);
N = jtAddDiscSys(N,2,Cd,1,Rd);
N = jtAddContSys(N,3,Hold,[],Rint);
N = jtAddResetDynamics(N,3,PhiReset,GammaReset,2);
N = jtCalcDynamics(N);

% Perform periodic analysis
N = jtBeginPeriodicAnalysis(N);
N = jtPassTime(N,h);
N = jtExecSys(N,2);  % Execute controller
N = jtResetSys(N,3); % Actuate through the hold
N = jtEndPeriodicAnalysis(N);


N = jtPassTime(N,h);
N = jtExecSys(N,2);
N = jtResetSys(N,3);


%%% ﻿Example 1: Simple Control Loop with Controller Activation %%%

% Define process, controller and sampler
P = ss(0,1,1,0);  % Integrator process
Qc = diag([1 0]); % State and control cost
R1c = 1;          % State noise
h = 1;            % Sampling period
CA = -1/h*(3+sqrt(3))/(2+sqrt(3)); % Minimum-variance controller
S = 1;            % Sampler

% Define the JitterTime model
N = jtInit;
N = jtAddContSys(N,1,P,3,R1c,Qc); % Add sys 1 (P)
N = jtAddDiscSys(N,2,S,1);        % Add sys 2 (S)
N = jtAddDiscSys(N,3,CA,2);       % Add sys 3 (CA)
N = jtCalcDynamics(N);

% Simulate the system and log the results
Nsteps = 6;       % Large time steps (control periods)
dt = h/10;        % Small time steps (for plotting)
l = 0;
tvec = [];pvec = []; Jvec = [];
for k = 1:Nsteps
  for j = 1:10  
    l = l+1;
    tvec(l) = N.Tsim; pvec(l) = N.P(1,1); Jvec(l) = N.J;
    N = jtPassTime(N,dt);
  end
  l = l+1;
  tvec(l) = N.Tsim; pvec(l) = N.P(1,1);	Jvec(l) = N.J;
  if k >= 3       % Activate controller after 3 samples
    N = jtExecSys(N,2);
    N = jtExecSys(N,3);
  end
end

% Plot the results
subplot(211)
plot(tvec,pvec)
xlabel('Time')
ylabel('Process variance')
subplot(212)
plot(tvec,Jvec)
xlabel('Time')
ylabel('Accumulated cost')


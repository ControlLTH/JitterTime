%%% ﻿Example 2: Asynchronous Kalman Filtering with Two Sensors %%%

% Define the process (inverted pendulum) and design objective
A = [0 1; 1 0];
B = [0; 1];
C = eye(2);
P = ss(A,B,C,0);           % Process dynamics
S1 = [1 0];                % Sensor 1 (position)
S2 = [0 1];                % Sensor 2 (velocity)
Qc = diag([1 0.1 0.1]);    % Cost function
R1c = B*B';                % Process input noise
R21 = 0.1;                 % Measurement noise on sensor 1
R22 = 0.01;                % Measurement noise on sensor 2

% Design state feedback and the two Kalman filters
h1 = 0.15;
[~,L,~,~,K1,sysd1] = lqgdesign(P(1),Qc,R1c,R21,h1);
h2bar = 0.5;
[~,~,~,~,K2,sysd2] = lqgdesign(P(2),Qc,R1c,R22,h2bar);

% Formulate Kalman filter based on y1
[Phi1e,Gam1e,C1e] = ssdata(sysd1);
IK1C1 = eye(size(Phi1e)) - K1*C1e;
Obs1 = ss(IK1C1*Phi1e,[IK1C1*Gam1e K1 0*K2],IK1C1*Phi1e,[IK1C1*Gam1e K1 0*K2],-1);

% Formulate Kalman filter based on y2
[Phi2e,Gam2e,C2e] = ssdata(sysd2);
IK2C2 = eye(size(Phi2e)) - K2*C2e;
Obs2 = ss(IK2C2,[0*Gam2e 0*K1 K2],IK2C2,[0*Gam2e 0*K1 K2],-1);

% Define the JitterTime model
N = jtInit;
N = jtAddContSys(N,1,P,5,R1c,Qc);          % Add sys 1 (Process)
N = jtAddDiscSys(N,2,S1,1,diag([R21 0]));  % Add sys 2 (Sensor 1)
N = jtAddDiscSys(N,3,S2,1,diag([0 R22]));  % Add sys 3 (Sensor 2)
N = jtAddDiscSys(N,4,{Obs1,Obs2},[5 2 3]); % Add sys 4 (Observers)
N = jtAddDiscSys(N,5,-L,4);                % Add sys 5 (Feedback)
N = jtCalcDynamics(N);

% Simulate the system
l = 0;
tvec = []; p1vec = []; p2vec = []; Jvec = [];
for k = 1:12/h1  % Simulate for 12 time units
  for j = 1:10   % Small time steps (for plotting)
    l = l+1;
    N = jtPassTime(N,h1/10);
    tvec(l) = N.Tsim; p1vec(l) = N.P(1,1); p2vec(l) = N.P(2,2); Jvec(l) = N.J;
    if rand < 0.05
      N = jtExecSys(N,3);    % Run the aperiodic sampler
      N = jtExecSys(N,4,2);  % Execute the aperiodic observer update
    end
  end
  N = jtExecSys(N,2);     % Run the periodic sampler
  N = jtExecSys(N,4,1);   % Execute the regular observer
  N = jtExecSys(N,5);     % Execute the state feedback
end

% Plot the result
subplot(211); plot(tvec,p1vec,tvec,p2vec)
xlabel('Time'); ylabel('State variance')
subplot(212); plot(tvec,Jvec)
xlabel('Time'); ylabel('Accumulated cost')
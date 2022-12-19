%%% DC motor LQG example

A = [0 1; 0 -1];
B = [0; 1];
C = [1 0];

R1c = B*B';
R2c = 0.1^2;
Q1c = diag([1.1 0.3]);
Q2c = 0.7;

L = lqr(A,B,Q1c,Q2c);
K = lqe(A,eye(2),C,R1c,R2c);

%%% 1. Calculate the continuous-time LQG cost

proc = ss(A,B,C,0);
ctrl = reg(proc,L,K);

N = jtInit;
N = jtAddContSys(N,1,proc,2,R1c,blkdiag(Q1c,Q2c));
N = jtAddContSys(N,2,ctrl,1,ctrl.b*R2c*ctrl.b');
N = jtCalcDynamics(N);

% Perform a periodic analysis to find stationary variance
T = 1; % Any time step works
N = jtBeginPeriodicAnalysis(N);
N = jtPassTime(N,T);  
N = jtEndPeriodicAnalysis(N);
N.P = N.Pperiodic;  % Set initial variance to stationary variance
N.J = 0;            % Set initial cost to zero
N = jtPassTime(N,T);
J_cont = N.J / T    % Calculate cost per time unit (1.4583)

%%% 1a. Baseline LQG controller with integral sampling (without resets, which
%%% leads to an unstable mode and prohibits periodic analysis)

hvec = [0.001 0.025:0.025:0.5];
J1vec = [];
J2vec = [];
for h = hvec 
  % Define extended plant with integral output
  Ae = [A zeros(2,1); C -sqrt(eps)];
  Be = [B; 0];
  Ce = [zeros(1,2) 1];
  R1e = blkdiag(R1c,R2c);
  R2e = 0;
  proce = ss(Ae,Be,Ce,0);
  
  %%% 2. ZOH-sampled LQG controller with reset integral sampling
  
  % Design Kalman filter and LQR
  [Phie,Re] = calcc2d(Ae,R1e,zeros(size(Ae)),h);
  Phie(end,end) = 0; % Reset output integrator in each sample
  [~,Game] = c2d(Ae,[B;0],h);
  Pe = dare(Phie',Ce',Re,0);
  Kf = Pe*Ce'/(Ce*Pe*Ce');  % Kalman filter gain (filter form)
  Ld = lqrd(A,B,Q1c,Q2c,h); % State feedback gain
  Le = [Ld 0];
  [Ar,Br,Cr,Dr] = dreg(Phie,Game,Ce,0,Le,Kf); % Form LQG controller
  ctrlr = -ss(Ar,Br,Cr,Dr,h);
  AReset = blkdiag(eye(size(A)),0); % Reset the sampler state
  
  N = jtInit;
  N = jtAddContSys(N,1,proce,2,blkdiag(R1c,R2c),blkdiag(Q1c,0,Q2c));
  N = jtAddResetDynamics(N,1,AReset);
  N = jtAddDiscSys(N,2,ctrlr,1);
  N = jtCalcDynamics(N);
  
  N = jtBeginPeriodicAnalysis(N);
  N = jtResetSys(N,1);  % Reset sampler integrator
  N = jtPassTime(N,h);  % Pass time h
  N = jtExecSys(N,2);   % Execute controller
  N = jtEndPeriodicAnalysis(N);
  
  N.P = N.Pperiodic;
  N.J = 0;
  N = jtResetSys(N,1);  % Reset sampler integrator
  N = jtPassTime(N,h);  % Pass time h
  J1b = N.dJ
  J1vec = [J1vec J1b];
  
  %%% 3. Mirkin's optimal sampled-data controller
  
  L = lqr(A,B,Q1c,Q2c);
  K = lqe(A,eye(2),C,R1c,R2c);
  
  Am = [A-B*L-K*C -B*L; K*C A];
  Bm = [K; -K];
  Cm = [-L -L];
  ctrlm = ss(Am,Bm,Cm,0);
  PhiReset = blkdiag(eye(size(A)),zeros(size(A)));
  
  N = jtInit;
  N = jtAddContSys(N,1,proc,2,R1c,blkdiag(Q1c,Q2c));
  N = jtAddContSys(N,2,ctrlm,1,ctrlm.b*R2c*ctrlm.b');
  N = jtAddResetDynamics(N,2,PhiReset);
  N = jtCalcDynamics(N);
  
  N = jtBeginPeriodicAnalysis(N);
  N = jtResetSys(N,2);
  N = jtPassTime(N,h);
  N = jtEndPeriodicAnalysis(N);
  
  N.P = N.Pperiodic;
  N.J = 0;
  N = jtResetSys(N,2);
  N = jtPassTime(N,h);
  J2 = N.dJ
  J2vec = [J2vec J2];

end

clf
plot(hvec,J1vec/J_cont)
hold on
plot(hvec,J2vec/J_cont)
xlabel('Sampling interval')
ylabel('Cost relative to continuous LQG')
legend('ZOH-sampled LQG control with reset sampling','Mirkin''s optimal sampled-data LQG')





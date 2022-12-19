function scheduling_init(Tsamp)

ttInitKernel('prioFP')

% Create two high-priority dummy tasks
data.E = 0.027 - sqrt(eps);
ttCreatePeriodicTask('dummy_task1', 0, 0.080, 'scheduling_dummycode', data);
ttSetPriority(1,'dummy_task1');
data.E = 0.044 - sqrt(eps);
ttCreatePeriodicTask('dummy_task2', 0, 0.140, 'scheduling_dummycode', data);
ttSetPriority(2,'dummy_task2');

% Create periodic timer and sampling handler
h = Tsamp;                          % Sampling period from init argument
ttCreateHandler('samp_handler', 1, 'scheduling_sampcode')
ttCreatePeriodicTimer('samp_timer', 0, h, 'samp_handler')

% Create low-priority controller task
data = [];
s = zpk('s');
P = 1/(s^2-1);                     % The process (inverted pendulum)
R1c = 1;                           % Continuous-time input noise
R2 = 0.01;                         % Discrete-time measurement noise
Qc = blkdiag(1,0.01);              % Continuous cost J = E(y^2 + 0.01*u^2)
tau = 0.025;                       % Assumed input-output delay (= E)
S = ss(1);                         % Sampler system
A = ss(1);                         % Actuator system
C = zpk(ss(lqgdesign(P,Qc,R1c,R2,h,tau),'min'));

global N
N = jtInit;                          % Initialize Jittersim
N = jtAddContSys(N,1,P,4,R1c,Qc);    % Add sys 1 (P), input from sys 4
N = jtAddDiscSys(N,2,S,1,R2);        % Add sys 2 (S), input from 1
N = jtAddDiscSys(N,3,C,2);           % Add sys 3 (C), input from 2
N = jtAddDiscSys(N,4,A,3);           % Add sys 4 (A), input from 3
N = jtCalcDynamics(N);               % Calculate the internal dynamics
N.samp = 0;
data.E = 0.025 - sqrt(eps);
ttCreateTask('ctrl_task', h, 'scheduling_ctrlcode', data)
ttSetPriority(3,'ctrl_task');

% J = 0.386

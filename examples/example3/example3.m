%%% Example 3 - Control Task Period Selection %%%

% This example requires the MatlabTrueTime and Jitterbug to be run. 

Tvec = 0.050:0.001:0.200; % Vector of periods to try
Jvec = [];

Tsim = 1000;
for Tsamp = Tvec          % Period of task 3
  Tsamp
	sim('scheduling',Tsim)  % Simulate the system for Tsim seconds
	J = JP.signals(1).values(end) / Tsim  % Calculate average cost
	Jvec = [Jvec J];
end

% Plot the results
plot(Tvec,Jvec)
xlabel('Task period')
ylabel('Controller cost')

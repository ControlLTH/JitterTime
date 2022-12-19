function [exectime,data] = scheduling_ctrlcode(seg,data)

global N

switch seg
	case 1
		% Pass time since last execution
		N = jtPassTimeUntil(N,ttCurrentTime);
 		ttAnalogOut(1,N.J);
    ttAnalogOut(2,trace(N.P));
      
		% Execute sampler
		N = jtExecSys(N,2);

% 		% Read sample from input buffer
% 		N.samp = N.samp - 1;
% 		if N.samp > 0
% 			disp(['Sample buffer overflow at ' num2str(ttCurrentTime)])
% 		end
		
		% Execute controller
		N = jtExecSys(N,3);
	
		exectime = data.E;
	case 2
		% Pass time since last execution
		N = jtPassTimeUntil(N,ttCurrentTime);
		ttAnalogOut(1,N.J);
    ttAnalogOut(2,trace(N.P));
		
		% Actuation - execute sys 4
		N = jtExecSys(N,4);
		exectime = -1;
end


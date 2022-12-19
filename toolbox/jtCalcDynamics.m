function N = jtCalcDynamics(N)
% N = jtCalcDynamics(N)
%
% Calculate the total system dynamics for the JitterTime system N. 
% The resulting system description is stored in the array N.sys.
%
% Arguments:
% N      The JitterTime system.

%
% Definition of N.sys
%
% State dynamics:
% Let the state of the total system (including all subsystems) be x.
%
% Continous time evolution:
% dx(t) = N.Ac x(t)dt + dvc(t)  where vc is continuous white noise with
%                               intensity N.Rc
%
% Discrete-event evolution when executing system S:
% x+(t_k) = S.Ad x(t_k) + vd(t_k) where vd is white Gaussian noise with
%                                 variance S.Rd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Count states and inputs and build a system-to-state
% index mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = zeros(2,length(N.sys)); % Indices into state vector for the subsystems
totstates = 0;   % Total number of states in JitterTime model
idtoindex = [];  % Array from subsys ID to index
for s = 1:length(N.sys)
	if (N.sys{s}.id <= 0)
		error('System ID is not a positive number.');
	end
	if (length(idtoindex) >= N.sys{s}.id && idtoindex(N.sys{s}.id) ~= 0)
		error('System ID %d defined twice.',N.sys{s}.id);
	end
	idtoindex(N.sys{s}.id) = s;
  n = N.sys{s}.n;
  ind(:,s) = [totstates+1; totstates+n];
  N.sys{s}.stateindex = ind(1,s):ind(2,s);
  totstates = totstates+n;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Check that the number of outputs and inputs match for
% all system connections (including reset dynamics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:length(N.sys)
  if (N.sys{s}.type == 1) % Continuous-time system
    ninputs = 0;          % Count total number of inputs
    for i = 1:length(N.sys{s}.inputid)
      in = idtoindex(N.sys{s}.inputid(i));
      ninputs = ninputs + N.sys{in}.p;
    end
    if ninputs > N.sys{s}.r
      error('System ID %d has more inputs in total (%d) than expected (%d)', ...
        N.sys{s}.id, ninputs, N.sys{s}.r)
    elseif ninputs < N.sys{s}.r
      warning(['System ID %d has fewer inputs in total (%d) than expected' ...
        '(%d)'], N.sys{s}.id, ninputs, N.sys{s}.r)
    end
    if N.sys{s}.resetDynamics
      for v = N.sys{s}.resetVersions
        ninputs = 0;          % Count total number of inputs
        for i = 1:length(N.sys{s}.resetInputid{v})
          in = idtoindex(N.sys{s}.resetInputid{v}(i));
          ninputs = ninputs + N.sys{in}.p;
        end
        resetp = size(N.sys{s}.resetB{v},2);
        if ninputs > resetp
          error(['System ID %d reset dynamics version %d has more inputs in'
            ' total (%d) than expected (%d)'], N.sys{s}.id, v, ninputs, resetp)
        elseif ninputs < resetp
          warning(['System ID %d reset dynamics version %d has fewer inputs' ...
          ' in total (%d) than expected (%d)'], N.sys{s}.id, v, ninputs, resetp)
        end
      end
    end
  elseif N.sys{s}.type == 2  % Discrete-time system
    for v = N.sys{s}.versions
      ninputs = 0;          % Count total number of inputs
      for i = 1:length(N.sys{s}.inputid{v})
        in = idtoindex(N.sys{s}.inputid{v}(i));
        ninputs = ninputs + N.sys{in}.p;
      end
      if ninputs > N.sys{s}.r
        error(['System ID %d version %d has more inputs in total (%d) than' ...
          ' expected (%d)'], N.sys{s}.id, v, ninputs, N.sys{s}.r)
      elseif ninputs < N.sys{s}.r
        warning(['System ID %d version %d has fewer inputs in total (%d)' ...
          ' than expected (%d)'], N.sys{s}.id, v, ninputs, N.sys{s}.r)
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Formulate the total continuous-time dynamics,
% continuous noise, and continuous cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N.Ac = zeros(totstates);  % total continuous state matrix
N.Rc = zeros(totstates);  % total continuous noise matrix
N.Qc = zeros(totstates);  % total continuous cost matrix

for s = 1:length(N.sys)
  % Construct the Ac matrix
  if (N.sys{s}.type == 1) % Continuous
    N.Ac(ind(1,s):ind(2,s),(ind(1,s):ind(2,s))) = N.sys{s}.A;
    bix = 1; % Index into B matrix
    for i = 1:length(N.sys{s}.inputid)
      in = idtoindex(N.sys{s}.inputid(i));
      Bpart = N.sys{s}.B(:,bix:bix+N.sys{in}.p-1); % Extract part of B matrix
      N.Ac(ind(1,s):ind(2,s),(ind(1,in):ind(2,in))) = Bpart*N.sys{in}.C;
      bix = bix + N.sys{in}.p;
    end
  end
  % Construct the Rc matrix
  if (N.sys{s}.type == 1) % Continuous
    N.Rc(ind(1,s):ind(2,s),(ind(1,s):ind(2,s))) = N.sys{s}.Rc;
  end
  % Construct the Qc matrix
  xtou = zeros(N.sys{s}.n,totstates);
  xtou(:,ind(1,s):ind(2,s)) = eye(N.sys{s}.n);
  if (N.sys{s}.type == 1) % 
    for i = 1:length(N.sys{s}.inputid)
      in = idtoindex(N.sys{s}.inputid(i));
      xtoy = zeros(N.sys{in}.p,totstates);
      xtoy(:,ind(1,in):ind(2,in)) = N.sys{in}.C;
      xtou = [xtou; xtoy];
    end
  end
  N.Qc = N.Qc + xtou'*N.sys{s}.Qc*xtou; % Combined cost of state, outputs and inputs 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Formulate all versions of the discrete-time dynamics,
% including system resets, and discrete-time noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:length(N.sys)
  if (N.sys{s}.type == 2) % Discrete-time system
    for v = 1:N.sys{s}.versions
      N.sys{s}.Ad{v} = eye(totstates);
      N.sys{s}.Ad{v}(ind(1,s):ind(2,s),ind(1,s):ind(2,s)) = N.sys{s}.A{v};
      N.sys{s}.Rd{v} = zeros(totstates);    
      N.sys{s}.Rd{v}(ind(1,s):ind(2,s),ind(1,s):ind(2,s)) = N.sys{s}.R{v};      
      bix = 1; % Index into B matrix
      for i = 1:length(N.sys{s}.inputid{v})
        in = idtoindex(N.sys{s}.inputid{v}(i));
        Bpart = N.sys{s}.B{v}(:,bix:bix+N.sys{in}.p-1); % Extract part of B
        N.sys{s}.Ad{v}(ind(1,s):ind(2,s),(ind(1,in):ind(2,in))) = ...
          Bpart*N.sys{in}.C;
        bix = bix + N.sys{in}.p;
      end
    end
  elseif (N.sys{s}.type == 1 && N.sys{s}.resetDynamics) % Cont. reset system
    for v = 1:N.sys{s}.resetVersions
      N.sys{s}.Ad{v} = eye(totstates);
      N.sys{s}.Ad{v}(ind(1,s):ind(2,s),ind(1,s):ind(2,s)) = N.sys{s}.resetA{v};    
      bix = 1; % Index into B matrix
      for i = 1:length(N.sys{s}.resetInputid{v})
        in = idtoindex(N.sys{s}.resetInputid{v}(i));
        Bpart = N.sys{s}.resetB{v}(:,bix:bix+N.sys{in}.p-1); % Extract part of B
        N.sys{s}.Ad{v}(ind(1,s):ind(2,s),(ind(1,in):ind(2,in))) = ...
          Bpart*N.sys{in}.C;
        bix = bix + N.sys{in}.p;   
      end
    end    
  end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 4: Initialize some data fields in N
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N.idtoindex = idtoindex;
N.m = zeros(length(N.Ac),1); % Initialize mean m to zero
N.P = zeros(size(N.Ac));     % Initialize variance P to zero
N.J = 0;
N.Tsim = 0;
N.periodicAnalysis = false;

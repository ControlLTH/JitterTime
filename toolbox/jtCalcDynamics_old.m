function N = jtCalcDynamics(N)
% N = jtCalcDynamics(N)
%
% Calculate the total system dynamics for the JitterTime system N. 
% The resulting system description is stored in N.nodes.
%
% Arguments:
% N      The JitterTime system.

% This function builds the timing node structure N.nodes which is
% usAd by calccost(). The starting point is the N.systems which is
% a set of linear continuous-time and discrete-time systems which
% are interconnected.
%
% Definition of N.sys
% ==============================================================
%
% State dynamics:
% Let the state of the system (including all subsystems) be x.
%
% Continous time evolution:
% dx(t) = N.Ac x(t)dt + dvc(t)  where vc is continuous white noise with
%                               intensity N.Rc
%
% Discrete-event evolution when executing system S:
% x+(t_k) = S.Ad*x(t_k) + vd(t_k) where vd is white Gaussian noise with
%                                 variance S.Rd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Count states and inputs and build a system-to-state
% index mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Indices into Ac for the subsystems
ind = zeros(2,length(N.sys));
% Indices into continuous outputs
indcontout = zeros(2,length(N.sys));
states = 0; % total # of states
contoutputs = 0;
Rc = [];    % total continuous noise matrix
idtoindex = []; % Array from subsys ID to index
for s = 1:length(N.sys)
	if (N.sys{s}.id < 1)
		error(sprintf('System ID %d is not >= 1.',N.sys{s}.id));
	end
	if ((length(idtoindex) >= N.sys{s}.id) && ...
			idtoindex(N.sys{s}.id) ~= 0)
		error(sprintf('System ID %d defined twice.', ...
			N.sys{s}.id));
	end
	idtoindex(N.sys{s}.id) = s;
	if (N.sys{s}.type == 1) % Continuous states
		n = size(N.sys{s}.A{1},1);
		ind(:,s) = [states+1; states+n];
		indcontout(:,s) = [contoutputs+1;contoutputs+size(N.sys{s}.C{1},1)];
		N.sys{s}.states = n;
		N.sys{s}.stateindex = ind(1,s):ind(2,s);
		states = states + n;
		contoutputs = contoutputs + size(N.sys{s}.C{1},1);
		Rc = blkdiag(Rc,N.sys{s}.Rc);
	end
	if (N.sys{s}.type == 2) % Discrete states = states + held outputs
		n = size(N.sys{s}.A{1},1)+size(N.sys{s}.C{1},1);
		ind(:,s) = [states+1; states+n];
		N.sys{s}.states = n;
		N.sys{s}.stateindex = [states+1:states+size(N.sys{s}.A{1},1)];
		N.sys{s}.outputindex = [states+size(N.sys{s}.A{1},1)+1:states+n];
		states = states + n;
		% Discrete systems have no continuous noise, so pad with zeros
		Rc = blkdiag(Rc,zeros(n));
	end
end

Ac = zeros(states);  % total continuous system matrix
Qc = zeros(states);  % total continuous cost matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Check that inputs matches outputs and build cell arrays
% for multiple inputs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:length(N.sys) % Check and fix interconnections
	for k = 1:N.sys{s}.versions  % Iterate over (discrete-time) versions
		insys = N.sys{s}.inputid{k};
		totinputs = 0; % Count total number of inputs
		cellB = {};
		cellD = {};
		% Map whole state vector to states in Q (state, output, input)
		xtou = zeros(N.sys{s}.states,states);
		xtou(:,ind(1,s):ind(2,s)) = eye(N.sys{s}.states);
    if (isempty(insys) || isequal(insys,0))
      % No input system – do nothing
    else
      for t = 1:length(insys)
        % No of outputs in input system
        if (insys(t) < 1 || insys(t) > length(idtoindex) || idtoindex(insys(t)) == 0)
          error(sprintf(['System ID %d as referred by system %d '...
            'not exist.']',insys(t),N.sys{s}.id));
        end
        inputs = size(N.sys{idtoindex(insys(t))}.C{1},1);
        % Split B and D into cell arrays of Bs and Ds for input systems
        if (size(N.sys{s}.B{k},2) >= totinputs+inputs)
          B = N.sys{s}.B{k}(:,totinputs+1:totinputs+inputs,:);
          D = N.sys{s}.D{k}(:,totinputs+1:totinputs+inputs,:);
          cellB = [cellB(:)' {B}];
          cellD = [cellD(:)' {D}];
          in = idtoindex(insys(t));
          if (N.sys{in}.type == 1)  % input sys is continuous
            xtoy = zeros(size(N.sys{in}.C{1},1),states);
            xtoy(:,ind(1,in):ind(2,in)) = N.sys{in}.C{1};
          end
          if (N.sys{in}.type == 2)  % input sys is discrete
            xtoy = zeros(N.sys{in}.outputs,states);
            xtoy(:,ind(1,in):ind(2,in)) = ...
              [zeros(N.sys{in}.outputs,size(N.sys{in}.A{1},1)) ...
              eye(N.sys{in}.outputs)];
          end
          xtou = [xtou; xtoy];
        end
        totinputs = totinputs + inputs;
      end  
      if (size(N.sys{s}.B{k},2) ~= totinputs)
        error(sprintf(['System ID %d: number of inputs (%d) does not equal total number of' ...
          ' ouputs in input systems (%d).'], N.sys{s}.id, size(N.sys{s}.B{k},2), totinputs));
      end
    end
    Qc = Qc + xtou'*N.sys{s}.Qc{k}*xtou; % Combined cost of state, outputs and inputs
		N.sys{s}.B{k} = cellB;
		N.sys{s}.D{k} = cellD;
	end
end

for s = 1:length(N.sys)
	if (N.sys{s}.type == 1) % Continuous
		Ac(ind(1,s):ind(2,s),(ind(1,s):ind(2,s))) = N.sys{s}.A{1};
		insys = N.sys{s}.inputid{1};
		for t = 1:length(insys)
      if (isempty(insys) || isequal(insys,0))
				% No input system - do nothing
			else
				in = idtoindex(insys(t));
				if (N.sys{in}.type == 1) % Continuous input
					Ac(ind(1,s):ind(2,s),(ind(1,in):ind(2,in))) = N.sys{s}.B{1}{t}*N.sys{in}.C{1};
				else % Discrete input (use output-state)
					Ac(ind(1,s):ind(2,s),(ind(1,in)+size(N.sys{in}.A{1},1):ind(2,in))) = N.sys{s}.B{1}{t};
				end
			end
		end
	end
end

% Save the total continuous-time dynamics, noise and cost matrices
N.Ac = Ac;
N.Rc = Rc;
N.Qc = Qc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Construct Ad matrices for discrete system execution,
% and corresponding noises Rd.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:length(N.sys)
	if N.sys{n}.type == 2
		for k = 1:N.sys{n}.versions  % Iterate over discrete-time versions
			N.sys{n}.Ad{k} = zeros(size(Ac,1));
			N.sys{n}.Rd{k} = zeros(size(Ac,1));
			
			for s = 1:length(N.sys)
				% Go through all discrete-time systems and set Ad and Rd matrices.
				if (s == n)
					insys = N.sys{s}.inputid{k};
					% Set Ad matrices
					N.sys{n}.Ad{k}(ind(1,s):ind(2,s),ind(1,s):ind(2,s)) = ...
						N.sys{n}.Ad{k}(ind(1,s):ind(2,s),ind(1,s):ind(2,s)) + ...
						[N.sys{s}.A{k} zeros(size(N.sys{s}.A{k},1),N.sys{s}.outputs); ...
						N.sys{s}.C{k} zeros(N.sys{s}.outputs)];
					% Set Rd
					N.sys{n}.Rd{k}(ind(1,s):ind(2,s),ind(1,s):ind(2,s)) = ...
						N.sys{n}.Rd{k}(ind(1,s):ind(2,s),ind(1,s):ind(2,s)) + ...
						N.sys{s}.R{k};
					
					% Now, go through all inputs to the system and alter the
					% Ad matrices accordingly
          for t = 1:length(insys)
            if (isempty(insys) || isequal(insys,0))
              % No input system – do nothing
            else
              in = idtoindex(insys(t));
              if (N.sys{in}.type == 1) % Continuous-time input system
                N.sys{s}.Ad{k}(ind(1,s):ind(2,s),(ind(1,in):ind(2,in))) ...
                  = N.sys{s}.Ad{k}(ind(1,s):ind(2,s),(ind(1,in):ind(2,in))) ...
                  +[N.sys{s}.B{k}{t} * N.sys{in}.C{1}; ...
                  N.sys{s}.D{k}{t} * N.sys{in}.C{1}];
              else % Discrete-time input system (use output-state)
                N.sys{n}.Ad{k}(ind(1,s):ind(2,s),(ind(1,in)+size(N.sys{in}.A{1},1):ind(2,in))) = ...
                  N.sys{n}.Ad{k}(ind(1,s):ind(2,s),(ind(1,in)+size(N.sys{in}.A{1},1):ind(2,in))) + ...
                  [N.sys{s}.B{k}{t};N.sys{s}.D{k}{t}];
              end
            end
          end
        else
          
					% Another subsystem than the executed one, set Ad = I
					N.sys{n}.Ad{k}(ind(1,s):ind(2,s),ind(1,s):ind(2,s)) = ...
						eye(N.sys{s}.states);
        end
        
      end
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Construct Ad matrices for continuous system resets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = 1:length(N.sys)
	if (N.sys{n}.type == 1 && N.sys{n}.resetDynamics)
		for k = 1:N.sys{n}.resetVersions  % Iterate over reset versions
			N.sys{n}.Ad{k} = eye(size(Ac,1));
			N.sys{n}.Rd{k} = zeros(size(Ac,1));		
      insys = N.sys{n}.inputid{k};
      % Set Ad matrices
      N.sys{n}.Ad{k}(ind(1,n):ind(2,n),ind(1,n):ind(2,n)) = ...
        N.sys{n}.resetPhi{k};
      if ~isempty(N.sys{n}.resetGamma) 
        if length(N.sys{n}.resetGamma) ~= 1
          error('Multiple input systems not yet implemented')
        end
        insys = N.sys{n}.resetInputid{k};
        if ~isempty(insys) 
          in = idtoindex(insys(1)); % Only one in system allowed
          if (N.sys{in}.type == 1) % Continuous-time input system
            N.sys{n}.Ad{k}(ind(1,n):ind(2,n),(ind(1,in):ind(2,in))) ...
              = N.sys{n}.Ad{k}(ind(1,n):ind(2,n),(ind(1,in):ind(2,in))) ...
              + N.sys{n}.resetGamma{k} * N.sys{in}.C{1};
          else % Discrete-time input system (use output-state)
            N.sys{n}.Ad{k}(ind(1,n):ind(2,n),(ind(1,in)+size(N.sys{in}.A{1},1):ind(2,in))) = ...
              N.sys{n}.Ad{k}(ind(1,n):ind(2,n),(ind(1,in)+size(N.sys{in}.A{1},1):ind(2,in))) + ...
              N.sys{n}.resetGamma{k};
          end
        end
      end     
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Initialize data fields in N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N.idtoindex = idtoindex;
N.m = zeros(length(N.Ac),1); % Initialize mean m to zero vector of correct size
N.P = zeros(size(N.Ac));     % Initialize variance P to zero matrix of correct size
N.J = 0;
N.Tsim = 0;
N.periodicAnalysis = false;

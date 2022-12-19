function N = jtAddDiscSys(N,sysid,sys,inputid,Rd,Qc)
% N = jtAddDiscSys(N,sysid,sys,inputid,Rd,Qc)
%
% Add a discrete-time linear system "sys" to the JitterTime system N.
% u(k) --> sys --> y(k)
%          x(k)           
% The state/output initially have mean E([x;y])=0 and covariance V([x;y])=0.
%
% Arguments:
% N        The JitterTime system to add this discrete-time system to.
% sysid    A unique positive ID number for this system (pick any). Used when
%          referred to from other systems. 
% sys      A discrete-time LTI system in state-space or transfer function form,
%          or a double/matrix (interpreted as a static gain transfer function).
%          Internally, the system is converted to state-space form, where the
%          held outputs are treated as additional states.
% inputid  A vector of input system IDs. The outputs of the corresponding
%          systems will be used as inputs to this system. The number of inputs
%          in this system must equal the total number of outputs in the input
%          systems.  An empty vector (or zero) indicates that the system inputs
%          are unconnected.
%
% Optional arguments (assumed zero if missing/empty):
% Rd       The discrete state/output/input noise covariance matrix. The noise
%          vector is assumed to have the same size as [x;y] (for state-
%          space systems) or [y;u] (for transfer-function systems). Noise is
%          added each time the system is executed.
% Qc       The continuous cost function is E(Int [x(t);y(t)]'*Qc*[x(t);y(t)] dt)
%          (for state-space systems) or E(Int y(t)'*Qc*y(t) dt) (for
%          transfer-function/zpk systems). Note that both x(t) and y(t) are
%          held constant between system executions.
%
% For multiple-version systems, sys, inputid, and R can be given as
% cell arrays (which must all have the same length, or length 1). Only
% state-space systems as allowed for multiple-version systems.

% Sanity checks
if (nargin < 4)
  error('To few arguments to function: N = jsAddDiscSys(N,sysid,sys,inputid)');
end

if iscell(sys)
	if ~isrow(sys)
		error('SYS must be a row')
	end 
	versions = length(sys);
	for v = 1:versions
		if ~isa(sys{v},'ss')
			error('System version %d is not in state-space form',v);
		end
	end
	origclass = 'ss';
else
	versions = 1;
	origclass = class(sys);
	switch origclass
		case 'ss'
		case 'tf'
		case 'zpk'
		case 'double'
			sys = tf(sys); % convert constant gain matrix to transfer function
			sys.Ts = -1;
			origclass = 'tf';
		otherwise
			error(['System class ' origclass ' not supported. SYS must be' ...
				' either ss, tf, zpk, or double.']);
	end
	sys = {sys};
end

for v = 1:versions
	sysk = sys{v};
	if ~isdt(sysk)
		error('System version %d is not discrete time.',v);
	end
	if ~isproper(sysk)
		error('System version %d is not proper.',v);
	end
	if hasdelay(sysk)
		error('System version %d is not delay free.',v);
	end
	if ~isa(sysk,'ss')
		sysk = ss(sysk,'min'); % Convert to state-space form if needed
  end
  [A,B,C,D] = ssdata(sysk);
	if v == 1
		n = size(A,2); % nbr of states
		p = size(C,1); % nbr of outputs
		r = size(B,2); % nbr of inputs
	else % Check size consistensy of all system versions 
		if (size(A,2) ~= n)
			error('System version %d has inconsistent number of states',v);
		elseif (size(C,1) ~= p)
			error('System version %d has inconsistent number of outputs',v);
		elseif (size(B,2) ~= r)
			error('System version %d has inconsistent number of inputs',v);
		end
  end
  Aarray{v} = [A zeros(n,p); C zeros(p)];
  Barray{v} = [B; D];
end

% Check inputid argument
if inputid == 0
  inputid = [];
end
if iscell(inputid)
	if length(inputid) ~= versions
		error('Length of INPUTID array must match length of SYS array')
  end
else % Copy inputid vector to each version
	inputid = squeeze(num2cell(repmat(inputid,[1 1 versions 1]),[1 2]))';
end
for v = 1:length(versions)
  if ismember(sysid,inputid{v})
    error('System cannot be connected to itself.')
  end
end

% Check Rd argument
switch origclass
	case 'ss'
		if (nargin < 5 || isempty(Rd))
			Rd = zeros(n+p);
		end
	case {'tf','zpk'}
		if (nargin < 5 || isempty(Rd))
			Rd = zeros(r+p);
		end
end

if iscell(Rd)
	if length(Rd) ~= versions
		error('Length of R array must match length of SYS array')
	end
else % Copy R matrix to each version
	Rd = squeeze(num2cell(repmat(Rd,[1 1 versions 1]),[1 2]))';
end
for v = 1:versions
	switch origclass
		case 'ss'
			if ~isequal(size(Rd{v}),[n+p n+p])
				error(['For state-space systems, the noise R should be a' ...
					' (n+p)*(n+p) matrix, where n is the number of states' ...
					' and p is the number of outputs.'])
			end
   case {'tf','zpk'} 	
		if (size(Rd{v},1) == r)
			Rd{v} = blkdiag(Rd{v},zeros(p));
		end
		if ~isequal(size(Rd{v}),[r+p r+p])
			error(['For transfer function systems, the noise R should be either' ...
				' a (r+p)*(r+p) matrix, where r is the number of inputs' ...
				' and p is the number of outputs; or a r*r matrix,' ...
				' where r is the input of inputs (assuming zero output noise)'])
		else
			Rd{v} = [Barray{v}]*Rd{v}(1:r,1:r)*[Barray{v}]'; 
		end
	end
end

% Check Qc argument
if nargin < 6 || isempty(Qc)
  switch origclass
    case 'ss'
      Qc = zeros(n+p);
    case {'tf','zpk'}
      Qc = zeros(p);
  end
end
if iscell(Qc)
  error('Qc is defined in continuous time and cannot have multiple versions')
end

switch origclass
  case 'ss'
    if ~isequal(size(Qc),[n+p n+p])
      error(['For state-space systems, the cost Qc should be a matrix' ...
        ' punishing all states and outputs: [x;y]^T*Qc*[x;y].'])
    end
  case {'tf','zpk'}
    if ~isequal(size(Qc),[p p])
      error(['For transfer function systems, the cost Qc should be' ...
        ' matrix punishing all outputs: y^T*Qc*y.'])
    end
    Qc = blkdiag(zeros(n),Qc);
end

S = struct('id',sysid,'type',2);
S.A = Aarray;
S.B = Barray;
S.C = [zeros(p,n) eye(p)];
S.R = Rd;
S.Qc = Qc;
S.inputid = inputid;
S.n = n+p;
S.r = r;
S.p = p;
S.origclass = origclass;
S.versions = versions;

N.sys = [N.sys(:)' {S}];

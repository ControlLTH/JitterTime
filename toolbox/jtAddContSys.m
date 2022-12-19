function N = jtAddContSys(N,sysid,sys,inputid,Rc,Qc)
% N = jtAddContSys(N,sysid,sys,inputid,Rc,Qc)
%
% Add a continuous-time linear system "sys" to the JitterTime system N.
% u(t) --> sys --> y(t)
%          x(t)
% The state initially has mean value E(x)=0 and covariance V(x)=0.
%
% Arguments:
% N        The JitterTime system to add this continuous-time system to.
% sysid    A unique positive ID number for this system (pick any). Used when
%          referred to from other systems. 
% sys      A strictly proper, delay-free continuous-time LTI system in
%          state-space or transfer function (or zpk) form. Internally, the
%          system will be converted to state-space form.
% inputid  A vector of input system IDs. The outputs of the corresponding
%          systems will be used as inputs to this system. The number of inputs
%          in this system must equal the total number of outputs in the input
%          systems. An empty vector (or zero) indicates that the system inputs
%          are unconnected.
%
% Optional arguments (assumed zero if missing/empty):
% Rc       The continuous state or input noise intensity matrix. The noise
%          vector is assumed to have the same size as x (for state-space
%          systems) or u (for transfer-function systems).
% Qc       The continuous cost function is E(Int [x(t);u(t)]'*Qc*[x(t);u(t)] dt)
%          (for state-space systems) or E(Int [y(t);u(t)]'*Qc*[y(t);u(t)] dt)
%          (for transfer-function/zpk systems).

% Sanity checks
if nargin < 4
  error('To few arguments to function: N = jtAddContSys(N,sysid,sys,inputid)');
end
origclass = class(sys);
switch origclass
 case 'ss'
 case 'tf'
 case 'zpk'
  otherwise
  error(['System class ' origclass ' not supported. SYS must be either ss,'...
    ' tf, or zpk.']);
end

if ~isct(sys)
  error('System is not continuous time.');
end

if ~isproper(sys)
  error('System is not proper.');
end

if hasdelay(sys)
	error('System has time delay, which is not supported.')
end

if ~isa(sys,'ss')
  sys = ss(sys);  % if tf or zpk, convert system to state-space form
end
[A,B,C,D] = ssdata(sys);
if norm(D) > eps
  error('The continuous system has a direct term, which is not supported.');
end

n = size(A,1); % number of states
r = size(B,2); % number of inputs
p = size(C,1); % number of outputs

if inputid == 0
  inputid = [];
end
if ismember(sysid,inputid)
  error('System cannot be connected to itself.')
end
  
if nargin < 5 || isempty(Rc)
  Rc = zeros(n);
else
  switch origclass
   case 'ss'
    if ~isequal(size(Rc),[n n])
      error(['For state-space systems, the noise Rc should be an' ...
	     ' n*n matrix, where n is the number of states.'])
    end
   case {'tf','zpk'}
     if ~isequal(size(Rc),[r r])
       error(['For transfer function systems, the noise Rc should' ...
         ' be an r*r matrix, where r is the number of inputs.'])
    else
      Rc = B*Rc*B';
    end
  end
end

if isempty(inputid)  % For unconnected systems, correct the input matrix
  B = zeros(n,0);    % and do not allow cost on the inputs
  r = 0;
end

if nargin < 6 || isempty(Qc)
  Qc = zeros(n+r);
else
  switch origclass
   case 'ss'
    if ~isequal(size(Qc),[n+r,n+r])
      error(['For state-space systems, the cost Qc should be an' ...
	     ' (n+r)*(n+r) matrix, where n is the number of states and' ...
       ' r is the number of inputs'])
    end
   case {'tf','zpk'}
    if ~isequal(size(Qc),[p+r,p+r]) 
      error(['For transfer function systems, the cost Qc should be an' ...
	     ' (p+r)*(p+r) matrix, where p is the number of outputs and' ...
       ' r is the number of inputs'])
    else
      xutoyu = blkdiag(C,eye(r));
      Qc = xutoyu'*Qc*xutoyu;
    end
  end
end

% Store all relevant data in a struct
S = struct('id',sysid,'type',1);

S.A = A;
S.B = B;
S.C = C;
S.Rc = Rc;
S.Qc = Qc;
S.inputid = inputid;
S.n = n;
S.r = r;
S.p = p;
S.origclass = origclass;
S.resetDynamics = false;

N.sys = [N.sys(:)' {S}];

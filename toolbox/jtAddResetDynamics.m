function N = jtAddResetDynamics(N,sysid,A,B,inputid)
% N = jtAddResetDynamics(N,sysid,A)
% N = jtAddResetDynamics(N,sysid,A,B,inputid)
%
% Specify reset dynamics A (and optionally B) for the continuous-time
% system sysid. When reset using jtResetSys(sysid[,ver]), the system state
% is momentarily updated as x(t+) = A * x(t) [+ B u(t)].
%
% Arguments:
% N        The JitterTime system.
% sysid    The ID number of a previously defined continuous-time system.
% A        State update matrix for the reset dynamics.
%
% Optional arguments (assumed zero if missing/empty):
% B        Input matrix for the reset dynamics.
% inputid  A vector of system IDs. The outputs of the corresponding systems
%          will be used as inputs for the reset update. The number of columns
%          in B must equal the total number of outputs in the input systems.
%
% For multiple-version dynamics, A, B and inputid can be given as cell arrays
% (which must all have the same length, or length 1, or empty).

% Sanity checks
if (nargin < 3)
  error('To few arguments to function: N = jtAddResetDynamics(N,sysid,A,B,inputid)')
end

base = 0;
for s = 1:length(N.sys)
  if (N.sys{s}.id == sysid && N.sys{s}.type == 1)
    base = s;
  end
end
if (base == 0)
  error('No continuous-time system with ID %d found. Define using jtAddContSys first.', sysid);
end
S = N.sys{base};
if ~isequal(S.origclass,'ss')
  error('System with ID %d not a state-space system', sysid)
end

if (nargin < 4)
  B = [];
end
if (nargin < 5)
  inputid = [];
end

if iscell(A)
	if ~isrow(sysd)
		error('A must be a row cell array')
  end
  if ~iscell(B) || ~isrow(B)
    error('B must be a row cell array')
  end
	versions = length(A);
  if length(B) ~= versions
    error('B array must have same length as A')
  end
else
  versions = 1;
  A = {A};
  B = {B};
end

for k = 1:versions
  if ~isequal(size(A{k}),[S.n S.n])
    error('A must have size %d x %d to match the original system',S.n,S.n)
  end
  if (~isempty(B{k}) && size(B{k},1) ~= S.n)
    error('B must have %d rows to match the original system',S.n)
  end
end

% Check inputid argument
if nargin < 5 || isequal(inputid,0)
  inputid = [];
end
if iscell(inputid)
	if (~isempty(inputid) && length(inputid) ~= versions)
		error('Length of inputid array must match length of Phi array')
  end
else % Copy inputid vector to each version
	inputid = squeeze(num2cell(repmat(inputid,[1 1 versions 1]),[1 2]))';
end
for v = 1:length(versions)
  if ismember(sysid,inputid{v})
    error('System cannot be connected to itself.')
  end
end

S.resetDynamics = true;
S.resetA = A;
S.resetB = B;
S.resetInputid = inputid;
S.resetVersions = versions;

N.sys{base} = S;

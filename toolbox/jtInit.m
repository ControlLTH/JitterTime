function N = jtInit
% N = jtInit
%
% Initialize a new JitterTime system.
%
% Return values:
% N         The JitterTime system which must be passed to all other functions.

if (nargin > 0)
  error(['Too many arguments to function N = jsInit.']);
end

N = struct('sys',0);
N.sys = cell(1,0);


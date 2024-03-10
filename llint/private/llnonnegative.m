function [odeFcn,thresholdNonNegative] = llnonnegative(ode,y0,threshold,idxNonNegative)  
%LLNONNEGATIVE  Helper function for handling nonnegative solution constraints
%   Modify the derivative function to prevent the solution from crossing zero.
%

neq = numel(y0);
thresholdNonNegative = [];
if any( (idxNonNegative < 1) | (idxNonNegative > neq) )
  error(message('llint:llnonnegative:NonNegativeIndicesInvalid'));
end
if any(y0(idxNonNegative) < 0)
  error(message('llint:llnonnegative:NonNegativeViolatedAtT0'));
end  
if length(threshold) == 1
  thresholdNonNegative = threshold(ones(size(idxNonNegative)));
else
  thresholdNonNegative = threshold(idxNonNegative);
end
thresholdNonNegative = thresholdNonNegative(:);
odeFcn = @local_odeFcn_nonnegative;   

% -----------------------------------------------------------
% Nested function: ODE with nonnegativity constraints imposed
%
  function yp = local_odeFcn_nonnegative(t,y,varargin)
    yp = feval(ode,t,y,varargin{:}); 
    ndx = idxNonNegative( find(y(idxNonNegative) <= 0) );
    yp(ndx) = max(yp(ndx),0);
  end  % local_odeFcn_nonnegative
% -----------------------------------------------------------

end  % llnonnegative


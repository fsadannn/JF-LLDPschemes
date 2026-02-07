function [o,varargout] = LLDPparams( options,neq )
%LLPARAMS Summary of this function goes here
%   Detailed explanation goes here
if isfield(options,'dKmax')   % s could be empty.
  kdmax = subsref(options, struct('type','.','subs','dKmax'));
  if isempty(kdmax)
    kdmax = min(30,neq);
  else
      kdmax = min(kdmax,neq);
  end
else
  kdmax = min(30,neq);
end

if isfield(options,'dKmin')   % s could be empty.
  kdmin = subsref(options, struct('type','.','subs','dKmin'));
  if isempty(kdmin)
    kdmin = 4;
  end
else
  kdmin = 4;
end

if isfield(options,'debug')   % s could be empty.
  debug = subsref(options, struct('type','.','subs','debug'));
  if isempty(debug)
    debug = 0;
  end
else
  debug = 0;
end

if isfield(options,'gamma')   % s could be empty.
  gamma = subsref(options, struct('type','.','subs','gamma'));
  if isempty(gamma)
    gamma = 0.001;
  end
else
  gamma = 0.001;
end

if isfield(options,'jacgap')   % s could be empty.
  jacgap = subsref(options, struct('type','.','subs','jacgap'));
  if isempty(jacgap)
    jacgap = 8;
  end
else
  jacgap = 8;
end

if isfield(options,'Jorder')   % s could be empty.
  jorder = subsref(options, struct('type','.','subs','Jorder'));
  if isempty(jorder)
    jorder = 1;
  end
else
  jorder = 1;
end

o = kdmax;
varargout{1}=kdmin;
varargout{2}=debug;
varargout{3}=gamma;
if nargout == 5
    varargout{4}=jacgap;
end
if nargout == 6
    varargout{5}=jorder;
end

end


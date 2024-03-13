function [neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, Fx, args, odeFcn, odeFxcn, ...
          options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, ...
          dataType ] =   ...
    llarguments(FcnHandlesUsed, FcnHandlesUsed2, solver, ode, odeFx, ...
                    tspan, y0, options, extras)
%ODEARGUMENTS  Helper function that processes arguments for LL solvers.
%

if FcnHandlesUsed  % function handles used
  if isempty(tspan) || isempty(y0) 
    error(message('llint:llarguments:TspanOrY0NotSupplied', solver));
  end      
  if length(tspan) < 2
    error(message('llint:llarguments:SizeTspan', solver));
  end  
  htspan = abs(tspan(2) - tspan(1));  
  tspan = tspan(:);
  ntspan = length(tspan);
  t0 = tspan(1);  
  next = 2;       % next entry in tspan
  tfinal = tspan(end);     
  args = extras;                 % use f(t,y,p1,p2...) 

else
  % Get default tspan and y0 from the function if none are specified.
  if isempty(tspan) || isempty(y0) 
    if exist(ode)==2 && ( nargout(ode)<3 && nargout(ode)~=-1 ) 
      error(message('llint:llarguments:NoDefaultParams', funstring( ode ), solver, funstring( ode )));      
    end
    [def_tspan,def_y0,def_options] = feval(ode,[],[],'init',extras{:});
    if isempty(tspan)
      tspan = def_tspan;
    end
    if isempty(y0)
      y0 = def_y0;
    end
    options = llset(def_options,options);
  end  
  tspan = tspan(:);
  ntspan = length(tspan);
  if ntspan == 1    % Integrate from 0 to tspan   
    t0 = 0;          
    next = 1;       % Next entry in tspan.
  else              
    t0 = tspan(1);  
    next = 2;       % next entry in tspan
  end
  htspan = abs(tspan(next) - t0);
  tfinal = tspan(end);   
  
  % The input arguments of f determine the args to use to evaluate f.
  if (exist(ode)==2)
    if (nargin(ode) == 2)           
      args = {};                   % f(t,y)
    else
      args = [{''} extras];        % f(t,y,'',p1,p2...)
    end
  else  % MEX-files, etc.
    try 
      args = [{''} extras];        % try f(t,y,'',p1,p2...)     
      feval(ode,tspan(1),y0(:),args{:});   
    catch
      args = {};                   % use f(t,y) only
    end
  end
end

if ~FcnHandlesUsed2  % function handles used
    % Get default tspan and y0 from the function if none are specified.
    
    if exist(odeFx)==2 && ( nargout(odeFx)<3 && nargout(odeFx)~=-1 )
        error(message('llint:llarguments:NoDefaultParams', funstring( odeFx ), solver, funstring( ode )));
    end
  
  % The input arguments of fx determine the args to use to evaluate fx.
  if (exist(odeFx)==2)
    if (nargin(odeFx) == 2)           
      args = {};                   % f(t,y)
    else
      args = [{''} extras];        % f(t,y,'',p1,p2...)
    end
  else  % MEX-files, etc.
    try 
      args = [{''} extras];        % try f(t,y,'',p1,p2...)     
      feval(odeFx,tspan(1),y0(:),args{:});   
    catch
      args = {};                   % use f(t,y) only
    end
  end
end

y0 = y0(:);
neq = length(y0);

% Test that tspan is internally consistent.
if any(isnan(tspan))
  error(message('llint:llarguments:TspanNaNValues'));
end
if t0 == tfinal
  error(message('llint:llarguments:TspanEndpointsNotDistinct'));
end
tdir = sign(tfinal - t0);
if any( tdir*diff(tspan) <= 0 )
  error(message('llint:llarguments:TspanNotMonotonic'));
end

f0 = feval(ode,t0,y0,args{:});
Fx = feval(odeFx,t0,y0,args{:});
[m,n] = size(f0);
if n > 1
  error(message('llint:llarguments:FoMustReturnCol', funstring( ode )));
elseif m ~= neq
    error(message('MATLAB:LLarguments:SizeIC', funstring( ode ), m, neq, funstring( ode )));
end

% Determine the dominant data type
classT0 = class(t0);
classY0 = class(y0);
classF0 = class(f0);
  
dataType = superiorfloat(t0,y0,f0);

if ~( strcmp(classT0,dataType) && strcmp(classY0,dataType) && ...
        strcmp(classF0,dataType))
    input1 = '''t0'', ''y0''';
    input2 = '''f(t0,y0)''';
    warning(message('llint:llarguments:InconsistentDataType',input1,input2,solver));
end


% Get the error control options, and set defaults.
rtol = odeget(options,'RelTol',1e-3,'fast');
if (length(rtol) ~= 1) || (rtol <= 0)
  error(message('llint:llarguments:RelTolNotPosScalar'));
end
if rtol < 100 * eps(dataType) 
  rtol = 100 * eps(dataType);
  warning(message('llint:llarguments:RelTolIncrease', sprintf( '%g', rtol )))
end
atol = odeget(options,'AbsTol',1e-6,'fast');
if any(atol <= 0)
  error(message('llint:llarguments:AbsTolNotPos'));
end
normcontrol = strcmp(odeget(options,'NormControl','off','fast'),'on');
if normcontrol
  if length(atol) ~= 1
    error(message('llint:llarguments:NonScalarAbsTol'));
  end
  normy = norm(y0);
else
  if (length(atol) ~= 1) && (length(atol) ~= neq)
    error(message('llint:llarguments:SizeAbsTol', funstring( ode ), neq)); 
  end
  atol = atol(:);
  normy = [];
end
threshold = atol / rtol;

% By default, hmax is 1/10 of the interval.
safehmax = 16.0*eps(dataType)*max(abs(t0),abs(tfinal));  % 'inf' for tfinal = inf
defaulthmax = max(0.1*(abs(tfinal-t0)), safehmax);
hmax = min(abs(tfinal-t0), abs(llget(options,'MaxStep',defaulthmax,'fast')));
if hmax <= 0
  error(message('llint:llarguments:MaxStepLEzero'));
end
htry = abs(odeget(options,'InitialStep',[],'fast'));
if ~isempty(htry) && (htry <= 0)
  error(message('llint:llarguments:InitialStepLEzero'));
end

odeFcn = ode;
odeFxcn = odeFx;

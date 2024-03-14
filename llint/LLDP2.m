function varargout = LLDP2(ode,tspan,y0,options,varargin)
%   LLDP2 Jacobian-free Locally Linearized Runge-Kutta 
%    method of Dormand and Prince with fixed order 1 finite difference in the Arnoldi algorithm
%
%   [TOUT,YOUT] = LLDP2(ODEFUN,TSPAN,Y0) with TSPAN = [T0 TFINAL] integrates 
%   the system of differential equations y' = f(t,y) from time T0 to TFINAL 
%   with initial conditions Y0. ODEFUN is a function handle. For a scalar T
%   and a vector Y, ODEFUN(T,Y) must return a column vector corresponding 
%   to f(t,y). Each row in the solution array YOUT corresponds to a time 
%   returned in the column vector TOUT.  To obtain solutions at specific 
%   times T0,T1,...,TFINAL (all increasing or all decreasing), use TSPAN = 
%   [T0 T1 ... TFINAL].     
%   
%   [TOUT,YOUT] = LLDP2(ODEFUN,TSPAN,Y0,OPTIONS) solves as above with default
%   integration properties replaced by values in OPTIONS, an argument created
%   with the llset function. See llset for details. Commonly used options 
%   are scalar relative error tolerance 'RelTol' (1e-3 by default) and vector
%   of absolute error tolerances 'AbsTol' (all components 1e-6 by default).
%   If certain components of the solution must be non-negative, use
%   ODESET to set the 'NonNegative' property to the indices of these
%   components.
%   
%   LLDP2 can solve problems M(t,y)*y' = f(t,y) with mass matrix M that is
%   nonsingular. Use ODESET to set the 'Mass' property to a function handle 
%   MASS if MASS(T,Y) returns the value of the mass matrix. If the mass matrix 
%   is constant, the matrix can be used as the value of the 'Mass' option. If
%   the mass matrix does not depend on the state variable Y and the function
%   MASS is to be called with one input argument T, set 'MStateDependence' to
%   'none'.
%
%   [TOUT,YOUT,TE,YE,IE] = LLDP2(ODEFUN,TSPAN,Y0,OPTIONS) with the 'Events'
%   property in OPTIONS set to a function handle EVENTS, solves as above 
%   while also finding where functions of (T,Y), called event functions, 
%   are zero. For each function you specify whether the integration is 
%   to terminate at a zero and whether the direction of the zero crossing 
%   matters. These are the three column vectors returned by EVENTS: 
%   [VALUE,ISTERMINAL,DIRECTION] = EVENTS(T,Y). For the I-th event function: 
%   VALUE(I) is the value of the function, ISTERMINAL(I)=1 if the integration 
%   is to terminate at a zero of this event function and 0 otherwise. 
%   DIRECTION(I)=0 if all zeros are to be computed (the default), +1 if only 
%   zeros where the event function is increasing, and -1 if only zeros where 
%   the event function is decreasing. Output TE is a column vector of times 
%   at which events occur. Rows of YE are the corresponding solutions, and 
%   indices in vector IE specify which event occurred.    
%
%   SOL = LLDP2(ODEFUN,[T0 TFINAL],Y0...) returns a structure that can be
%   used with DEVAL to evaluate the solution or its first derivative at 
%   any point between T0 and TFINAL. The steps chosen by LLDP2 are returned 
%   in a row vector SOL.x.  For each I, the column SOL.y(:,I) contains 
%   the solution at SOL.x(I). If events were detected, SOL.xe is a row vector 
%   of points at which events occurred. Columns of SOL.ye are the corresponding 
%   solutions, and indices in vector SOL.ie specify which event occurred. 
%   Details are to be found in The MATLAB ODE Suite, L. F. Shampine and
%   M. W. Reichelt, SIAM Journal on Scientific Computing, 18-1, 1997.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-14-94
%   Copyright 1984-2017 The MathWorks, Inc.
%   Modified version of ode45 according  to
%   Jacobian-free Locally Linearized Runge-Kutta method of Dormand and Prince
%   for large systems of differential equations by F.S. Naranjo-Noda and J.C. Jimenez
%   Copyright (c) 2022, Frank S. Naranjo-Noda
solver_name = 'LLDP2';

% Check inputs
if nargin < 4
  options = [];
  if nargin < 3
    y0 = [];
    if nargin < 2
      tspan = [];
      if nargin < 1
        error(message('llint:LLDP_Kphi1:NotEnoughInputs'));
      end
    end
  end
end

% Stats
nsteps  = 0;
nfailed = 0;
nfevals = 0;
% nJevals = 0;

% Output
FcnHandlesUsed  = isa(ode,'function_handle');
output_sol = (FcnHandlesUsed && (nargout==1));      % sol = odeXX(...)
output_ty  = (~output_sol && (nargout > 0));  % [t,y,...] = odeXX(...)
% There might be no output requested...

sol = []; f3d = [];
if output_sol
  sol.solver = solver_name;
  sol.extdata.odefun = ode;
  sol.extdata.options = options;
  sol.extdata.varargin = varargin;
end

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, odeFcn, ...
  options, threshold, rtol, normcontrol, normy, hmax, htry, htspan, dataType] = ...
  llargumentsfj(FcnHandlesUsed, solver_name, ode, tspan, ...
                y0, options, varargin);
[kdmax,kdmin,~,~,~,jorder]=LLDPparams(options,neq);
gamma=0.001;
nfevals = nfevals + 1;

% Handle the output
if nargout > 0
  outputFcn = llget(options,'OutputFcn',[],'fast');
else
  outputFcn = llget(options,'OutputFcn',@odeplot,'fast');
end
outputArgs = {};
if isempty(outputFcn)
  haveOutputFcn = false;
else
  haveOutputFcn = true;
  outputs = llget(options,'OutputSel',1:neq,'fast');
  if isa(outputFcn,'function_handle')
    % With MATLAB 6 syntax pass additional input arguments to outputFcn.
    outputArgs = varargin;
  end
end
refine = max(1,llget(options,'Refine',4,'fast'));
if ntspan > 2
  outputAt = 1;          % output only at tspan points
elseif refine <= 1
  outputAt = 2;          % computed points, no refinement
else
  outputAt = 3;          % computed points, with refinement
  S = (1:refine-1) / refine;
end
printstats = strcmp(llget(options,'Stats','off','fast'),'on');

% Handle the event function
[haveEventFcn,eventFcn,eventArgs,valt,teout,yeout,ieout] = ...
  llevents(FcnHandlesUsed,odeFcn,t0,y0,options,varargin);

% Handle the mass matrix
[Mtype, M, Mfun] =  llmass(FcnHandlesUsed,odeFcn,t0,y0,options,varargin);
if Mtype > 0  % non-trivial mass matrix
  Msingular = llget(options,'MassSingular','no','fast');
  if strcmp(Msingular,'maybe')
    warning(message('llint:LLDP_Kphi1:MassSingularAssumedNo'));
  elseif strcmp(Msingular,'yes')
    error(message('llint:LLDP_Kphi1:MassSingularYes'));
  end
  % Incorporate the mass matrix into odeFcn and odeArgs.
  [odeFcn,odeArgs] = llmassexplicit(FcnHandlesUsed,Mtype,odeFcn,odeArgs,Mfun,M);
  f0 = feval(odeFcn,t0,y0,odeArgs{:});
  nfevals = nfevals + 1;
end

% Non-negative solution components
idxNonNegative = llget(options,'NonNegative',[],'fast');
nonNegative = ~isempty(idxNonNegative);
if nonNegative  % modify the derivative function
  [odeFcn,thresholdNonNegative] = llnonnegative(odeFcn,y0,threshold,idxNonNegative);
  f0 = feval(odeFcn,t0,y0,odeArgs{:});
  nfevals = nfevals + 1;
end

t = t0;
y = y0;

% Allocate memory if we're generating output.
nout = 0;
tout = []; yout = [];
if nargout > 0
  if output_sol
    chunk = min(max(100,50*refine), refine+floor((2^11)/neq));
    tout = zeros(1,chunk,dataType);
    yout = zeros(neq,chunk,dataType);
    f3d  = zeros(neq,7,chunk,dataType);
  else
    if ntspan > 2                         % output only at tspan points
      tout = zeros(1,ntspan,dataType);
      yout = zeros(neq,ntspan,dataType);
    else                                  % alloc in chunks
      chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
      tout = zeros(1,chunk,dataType);
      yout = zeros(neq,chunk,dataType);
    end
  end
  nout = 1;
  tout(nout) = t;
  yout(:,nout) = y;
end

% Initialize method parameters.
pow = 1/5;
A = [1/5, 3/10, 4/5, 8/9, 1, 1]; % Still used by restarting criteria
% B = [
%     1/5         3/40    44/45   19372/6561      9017/3168       35/384
%     0           9/40    -56/15  -25360/2187     -355/33         0
%     0           0       32/9    64448/6561      46732/5247      500/1113
%     0           0       0       -212/729        49/176          125/192
%     0           0       0       0               -5103/18656     -2187/6784
%     0           0       0       0               0               11/84
%     0           0       0       0               0               0
%     ];
% E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];

% Same values as above extracted as scalars (1 and 0 values ommitted)
a2=cast(1/5,dataType);
a3=cast(3/10,dataType);
a4=cast(4/5,dataType);
a5=cast(8/9,dataType);

% b11=cast(1/5,dataType);
% b21=cast(3/40,dataType);
% b31=cast(44/45,dataType);
% b41=cast(19372/6561,dataType);
% b51=cast(9017/3168,dataType);
% b61=cast(35/384,dataType);
b22=cast(9/40,dataType);
b32=cast(-56/15,dataType);
b42=cast(-25360/2187,dataType);
b52=cast(-355/33,dataType);
b33=cast(32/9,dataType);
b43=cast(64448/6561,dataType);
b53=cast(46732/5247,dataType);
b63=cast(500/1113,dataType);
b44=cast(-212/729,dataType);
b54=cast(49/176,dataType);
b64=cast(125/192,dataType);
b55=cast(-5103/18656,dataType);
b65=cast(-2187/6784,dataType);
b66=cast(11/84,dataType);

% e1=cast(71/57600,dataType);
e3=cast(-71/16695,dataType);
e4=cast(71/1920,dataType);
e5=cast(-17253/339200,dataType);
e6=cast(22/525,dataType);
e7=cast(-1/40,dataType);

hmin = 16*eps(t);
if isempty(htry)
  % Compute an initial step size h using y'(t).
  absh = min(hmax, htspan);
  if normcontrol
    rh = (norm(f0) / max(normy,threshold)) / (0.8 * rtol^pow);
  else
    rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8 * rtol^pow);
  end
  if absh * rh > 1
    absh = 1 / rh;
  end
  absh = max(absh, hmin);
else
  absh = min(hmax, max(hmin, htry));
end
f1 = f0;

% Initialize the output function.
if haveOutputFcn
  feval(outputFcn,[t tfinal],y(outputs),'init',outputArgs{:});
end

% Cleanup the main ode function call
FcnUsed = isa(odeFcn,'function_handle');
odeFcn_main = odefcncleanup(FcnUsed,odeFcn,odeArgs);

% THE MAIN LOOP

done = false;
atol=options.AbsTol;
Fxphi = zeros(neq, 5);

%krilov consants and parameters
kmin = kdmax;
kmax = kdmin;
pmin = 6;
pmax = 3;

%stats
ksum = 0;
nexpo=0;
breakcont=0;
reusedsubspace=0;
korder_krylov = 100;
korder_rk = 100;
alpha2 = 100;
korder_krylov_max = 0;
korder_rk_max = 0;
alpha2_max = 0;
eta_v_2_cont = 0; 


%adaptive for krilov
kdim=kdmin;
fac=1/log(2);
%adaptive Jacobian

while ~done

  % By default, hmin is a small number such that t+hmin is only slightly
  % different than t.  It might be 0 if t is 0.
  hmin = 16*eps(t);
  absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
  h = tdir * absh;

  % Stretch the step if within 10% of tfinal-t.
  if 1.1*absh >= abs(tfinal - t)
    h = tfinal - t;
    absh = abs(h);
    done = true;
  end
  ny = norm(y);

  % LOOP FOR ADVANCING ONE STEP.
  nofailed = true;                      % no failed attempts
  while true
     
     

     % matrix exponential calculation
     if nofailed
         [phi,kerror,kdim,nexpcont,breakdown,knfeval,padepq,korder] = phi1LLDP2(odeFcn_main,f1,h,hmin,t,y,ny,kdim,...
                rtol,atol,kdmax,kdmin,gamma);
         nfevals=nfevals+knfeval;
     else
         [phi,kerror,kdim,nexpcont,breakdown,knfeval,padepq,korder] = phi1LLDP2(odeFcn_main,f1,h,hmin,t,y,ny,kdim,...
                rtol,atol,kdmax,kdmin,gamma, 1);
         nfevals=nfevals+knfeval;
        reusedsubspace = reusedsubspace+1;
     end

     % update information
     ksum = ksum + kdim;
     breakcont=breakcont+breakdown;
     nexpo = nexpo + nexpcont;
    

    ny2 = sqrt((1+ny)*eps);
    nyp1 =  ny2/(norm(phi(:,1))+eps);
    nyp2 =  ny2/(norm(phi(:,2))+eps);
    nyp3 =  ny2/(norm(phi(:,3))+eps);
    nyp4 =  ny2/(norm(phi(:,4))+eps);
    nyp5 =  ny2/(norm(phi(:,5))+eps);

        
    logh=log(h);
    alpha_2 = min([log(nyp1),log(nyp2),log(nyp3),log(nyp4),log(nyp5)])/logh;
    eta2 = 1;

    if alpha_2 <3 && (alpha_2+2)<korder
    
        Fxphi(:,1) = FreeJ_f_w(2,odeFcn_main,y,phi(:,1),t,nyp1);
        Fxphi(:,2) = FreeJ_f_w(2,odeFcn_main,y,phi(:,2),t,nyp2);
        Fxphi(:,3) = FreeJ_f_w(2,odeFcn_main,y,phi(:,3),t,nyp3);
        Fxphi(:,4) = FreeJ_f_w(2,odeFcn_main,y,phi(:,4),t,nyp4);
        Fxphi(:,5) = FreeJ_f_w(2,odeFcn_main,y,phi(:,5),t,nyp5);
        nfevals = nfevals + 10;
        eta2 = 2;

    else
        Fxphi(:,1) = (1/nyp1).*(odeFcn_main(t, y+nyp1.*phi(:,1))-f1);
        Fxphi(:,2) = (1/nyp2).*(odeFcn_main(t, y+nyp2.*phi(:,2))-f1);
        Fxphi(:,3) = (1/nyp3).*(odeFcn_main(t, y+nyp3.*phi(:,3))-f1);
        Fxphi(:,4) = (1/nyp4).*(odeFcn_main(t, y+nyp4.*phi(:,4))-f1);
        Fxphi(:,5) = (1/nyp5).*(odeFcn_main(t, y+nyp5.*phi(:,5))-f1);
        nfevals = nfevals + 5;
        eta_v_2_cont = eta_v_2_cont + 1;
        
    end
    

    yLL2 = y + phi(:,5);
    common = f1 + Fxphi(:,5);

    y2 = y + phi(:,1);
    t2 = t + h*a2;
    f2 = odeFcn_main(t2, y2)- f1 - Fxphi(:,1);

    y3 = y + (h*b22).*f2 + phi(:,2);
    t3 = t + h*a3;
    f3 = odeFcn_main(t3, y3)- f1 - Fxphi(:,2);

    y4 = y + (h*b32).*f2 + (h*b33).*f3 + phi(:,3);
    t4 = t + h*a4;
    f4 = odeFcn_main(t4, y4)- f1 - Fxphi(:,3);

    y5 = y + (h*b42).*f2 + (h*b43).*f3 + (h*b44).*f4 + phi(:,4);
    t5 = t + h*a5;
    f5 = odeFcn_main(t5, y5)- f1 - Fxphi(:,4);

    y6 = yLL2 + (h*b52).*f2 + (h*b53).*f3 + (h*b54).*f4 + (h*b55).*f5;
    t6 = t + h;
    fprev = odeFcn_main(t6, y6);
    f6 = fprev - common;

    tnew = t + h;
    if done
      tnew = tfinal;   % Hit end point exactly.
    end
    h = tnew - t;      % Purify h.

    ynew = yLL2 +  (h*b63).*f3 + (h*b64).*f4 + (h*b65).*f5 + (h*b66).*f6 ;
    F_ynew = odeFcn_main(tnew,ynew);
    f7 = F_ynew - common;

    nfevals = nfevals + 6;
    
    
    rk = (1+alpha_2)*eta2+1;
%     korderp =  max([1,min([rk,korder,4])])+1;
%     pow = 1/korderp;
    
%  disp([rk,korder,korderp])

    % Estimate the error.
    NNrejectStep = false;
    if normcontrol
      normynew = norm(ynew);
      errwt = max(max(normy,normynew),threshold);
      err = absh * (norm(e3.*f3 + e4.*f4 + e5.*f5 + e6.*f6 + e7.*f7) / errwt);
      if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
        errNN = norm( max(0,-ynew(idxNonNegative)) ) / errwt ;
        if errNN > rtol
          err = errNN;
          NNrejectStep = true;
        end
      end
    else
      err = absh * norm((e3.*f3 + e4.*f4 + e5.*f5 + e6.*f6 + e7.*f7) ./ max(max(abs(y),abs(ynew)),threshold),inf);
      if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
        errNN = norm( max(0,-ynew(idxNonNegative)) ./ thresholdNonNegative, inf);
        if errNN > rtol
          err = errNN;
          NNrejectStep = true;
        end
      end
    end

    % Accept the solution only if the weighted error is no more than the
    % tolerance rtol.  Estimate an h that will yield an error of rtol on
    % the next step or the next try at taking this step, as the case may be,
    % and use 0.8 of this value to avoid failures.
    if err > rtol                       % Failed step
      nfailed = nfailed + 1;
      if absh <= hmin
        warning(message('llint:LLDP_Kphi1:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));
        solver_output = llfinalize(solver_name, sol,...
                                    outputFcn, outputArgs,...
                                    printstats, [nsteps, nfailed, nfevals],...
                                    nout, tout, yout,...
                                    haveEventFcn, teout, yeout, ieout,...
                                    {f3d,idxNonNegative});
        if nargout > 0
          varargout = solver_output;
        end
        return;
      end

      if nofailed
        nofailed = false;
        if NNrejectStep
          absh = max(hmin, 0.5*absh);
        else
          absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
        end
      else
        absh = max(hmin, 0.5 * absh);
      end
      h = tdir * absh;
      done = false;

      rfacmax=max(1,kdim/3);
      rfacmin=1;
      kdnew = ceil(kdim + max(rfacmax,min(log(kerror/gamma)*fac,rfacmin)));
      kdim=max(kdmin,min(kdmax,kdnew));


    else                                % Successful step

      NNreset_f7 = false;
      if nonNegative && any(ynew(idxNonNegative)<0)
        ynew(idxNonNegative) = max(ynew(idxNonNegative),0);
        if normcontrol
          normynew = norm(ynew);
        end
        NNreset_f7 = true;
      end
      
      kmin = min(kdim,kmin);
      kmax = max(kmax,kdim);
      pmin = min(padepq,pmin);
      pmax = max(padepq,pmax);

      afacmin=kdim/3;
      afacmax=-kdim/4;
      kdnew = floor(kdim + max(afacmax,min(log(kerror/gamma)*fac,afacmin)));
      kdim=max(kdmin,min(kdmax,kdnew));
      
      
      korder_krylov = min(korder_krylov, korder);
      korder_krylov_max = max(korder_krylov_max, korder);
      
      korder_rk = min(rk,korder_rk);
      korder_rk_max = max(rk,korder_rk_max);
      
      alpha2 = min(alpha_2,alpha2);
      alpha2_max = max(alpha_2,alpha2_max);
      
      break;

    end
    
  end
  nsteps = nsteps + 1;

  if haveEventFcn
    f = [zeros(length(f2),1) f2 f3 f4 f5 f6 f7];
    [te,ye,ie,valt,stop] = ...
        odezero(@ntrpLLRK45_Krilov_fj,eventFcn,eventArgs,valt,t,y,t0,h,f,idxNonNegative,odeFcn_main,f1,kdmax);
    if ~isempty(te)
      if output_sol || (nargout > 2)
        teout = [teout, te];
        yeout = [yeout, ye];
        ieout = [ieout, ie];
      end
      if stop               % Stop on a terminal event.
        % Adjust the interpolation data to [t te(end)].

        % Update the derivatives using the interpolating polynomial.
        taux = t + (te(end) - t)*A;
        [~,f(:,2:7)] = ntrpLLRK45_Krilov_fj(taux,t,y,h,f,idxNonNegative,odeFcn_main,f1,kdmax);
        f2 = f(:,2); f3 = f(:,3); f4 = f(:,4); f5 = f(:,5); f6 = f(:,6); f7 = f(:,7);

        tnew = te(end);
        ynew = ye(:,end);
        h = tnew - t;
        done = true;
      end
    end
  end

  if output_sol
    nout = nout + 1;
    if nout > length(tout)
      tout = [tout, zeros(1,chunk,dataType)];  % requires chunk >= refine
      yout = [yout, zeros(neq,chunk,dataType)];
%       f3d  = cat(3,f3d,zeros(neq,7,chunk,dataType));
    end
    tout(nout) = tnew;
    yout(:,nout) = ynew;
%     f3d(:,:,nout) = [zeros(length(f2),1) f2 f3 f4 f5 f6 f7];
  end

  if output_ty || haveOutputFcn
    switch outputAt
     case 2      % computed points, no refinement
      nout_new = 1;
      tout_new = tnew;
      yout_new = ynew;
     case 3      % computed points, with refinement
      tref = t + (tnew-t)*S;
      nout_new = refine;
      tout_new = [tref, tnew];
      yntrp45 = ntrpLLRK45_Krilov_fj(tref,t,y,h,[zeros(length(f2),1) f2 f3 f4 f5 f6 f7],idxNonNegative,odeFcn_main,f1,kdmax);
      yout_new = [yntrp45, ynew];
     case 1      % output only at tspan points
      nout_new =  0;
      tout_new = [];
      yout_new = [];
      while next <= ntspan
        if tdir * (tnew - tspan(next)) < 0
          if haveEventFcn && stop     % output tstop,ystop
            nout_new = nout_new + 1;
            tout_new = [tout_new, tnew];
            yout_new = [yout_new, ynew];
          end
          break;
        end
        nout_new = nout_new + 1;
        tout_new = [tout_new, tspan(next)];
        if tspan(next) == tnew
          yout_new = [yout_new, ynew];
        else
          yntrp45 = ntrpLLRK45_Krilov_fj(tspan(next),t,y,h,[zeros(length(f2),1) f2 f3 f4 f5 f6 f7],idxNonNegative,odeFcn_main,f1,kdmax);
          yout_new = [yout_new, yntrp45];
        end
        next = next + 1;
      end
    end

    if nout_new > 0
      if output_ty
        oldnout = nout;
        nout = nout + nout_new;
        if nout > length(tout)
          tout = [tout, zeros(1,chunk,dataType)];  % requires chunk >= refine
          yout = [yout, zeros(neq,chunk,dataType)];
        end
        idx = oldnout+1:nout;
        tout(idx) = tout_new;
        yout(:,idx) = yout_new;
      end
      if haveOutputFcn
        stop = feval(outputFcn,tout_new,yout_new(outputs,:),'',outputArgs{:});
        if stop
          done = true;
        end
      end
    end
  end

  if done
    break
  end

  % If there were no failures compute a new h.
  if nofailed
    % Note that absh may shrink by 0.8, and that err may be 0.
    temp = 1.25*(err/rtol)^pow;
    if temp > 0.2
      absh = absh / temp;
    else
      absh = 5.0*absh;
    end
  end

  % Advance the integration one step.
  t = tnew;
  y = ynew;
  if normcontrol
    normy = normynew;
  end
  if NNreset_f7
    % Used f7 for unperturbed solution to interpolate.
    % Now reset f7 to move along constraint.
%     f7 = odeFcn_main(tnew,ynew);
    F_ynew= deFcn_main(tnew,ynew) - common;
    nfevals = nfevals + 1;
  end
  f1 = F_ynew;  % Already have f(tnew,ynew)
%   nJevals = nJevals + 1;
end

solver_output = llfinalize(solver_name, sol,...
                            outputFcn, outputArgs,...
                            printstats, [nsteps, nfailed, nfevals],...
                            nout, tout, yout,...
                            haveEventFcn, teout, yeout, ieout,...
                            {f3d,idxNonNegative});
if nargout == 1
    solver_output{1}.stats.nJevals=0;
    solver_output{1}.stats.nexpm=nexpo;
    solver_output{1}.stats.Kdim_sum = ksum;
    solver_output{1}.stats.Kdim_min = kmin;
    solver_output{1}.stats.Kdim_max = kmax;
    solver_output{1}.stats.pade_min = pmin;
    solver_output{1}.stats.pade_max = pmax;
    solver_output{1}.stats.breakdown = breakcont;
    solver_output{1}.stats.reusedsubspace = reusedsubspace;
    solver_output{1}.stats.korder =  min([korder_rk,korder_krylov,5]);
    solver_output{1}.stats.korder_krylov = korder_krylov;
    solver_output{1}.stats.eta = 1;
    solver_output{1}.stats.alpha = korder_krylov-1;
    solver_output{1}.stats.korder_rk = korder_rk;
    solver_output{1}.stats.alpha2 = alpha2;
    solver_output{1}.stats.korder_krylov_max = korder_krylov_max;
    solver_output{1}.stats.eta_max = 1;
    solver_output{1}.stats.alpha_max = korder_krylov_max-1;
    solver_output{1}.stats.korder_rk_max = korder_rk_max;
    solver_output{1}.stats.alpha2_max = alpha2_max;
    solver_output{1}.stats.eta_v_2_cont = eta_v_2_cont;
%     if dextinf~=0
%     extinf = struct('kdim',kdimarr);
%         solver_output{1}.extinf = extinf;
%     end
end
%debug

if nargout > 0
  varargout = solver_output;
end

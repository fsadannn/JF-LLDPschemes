function [yinterp,ypinterp] = ntrpLLRK45_Krilov(tinterp,t,y,h,f,idxNonNegative,Fx,F,kdim)
%NTRP45  Interpolation helper function for ODE45.
%   YINTERP = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F,IDX) uses data computed in ODE45
%   to approximate the solution at time TINTERP.  TINTERP may be a scalar 
%   or a row vector. 
%   The arguments TNEW and YNEW do not affect the computations. They are 
%   required for consistency of syntax with other interpolation functions. 
%   Any values entered for TNEW and YNEW are ignored.
%    
%   [YINTERP,YPINTERP] = NTRP45(TINTERP,T,Y,TNEW,YNEW,H,F,IDX) returns also the
%   derivative of the polynomial approximating the solution. 
%
%   IDX has indices of solution components that must be non-negative. Negative 
%   YINTERP(IDX) are replaced with zeros and the derivative YPINTERP(IDX) is 
%   set to zero.
%   
%   See also ODE45, DEVAL.

%   Mark W. Reichelt and Lawrence F. Shampine, 6-13-94
%   Copyright 1984-2009 The MathWorks, Inc.
%   Modified version for Local Lonearized Dormand and Prince Runge-Kutta method
%   Copyright (c) 2022, Frank S. Naranjo-Noda

BI = [
    1       -183/64      37/12       -145/128
    0          0           0            0
    0       1500/371    -1000/159    1000/371
    0       -125/32       125/12     -375/64 
    0       9477/3392   -729/106    25515/6784
    0        -11/7        11/3        -55/28
    0         3/2         -4            5/2
    ];

neq = length(y);
u = tinterp - t;
LL=zeros(neq,length(tinterp));
for i=1:length(tinterp) 
   tt=u(i);
   LL(:,i)=phi1LLDP_single(Fx,F,tt,kdim);
end
s = u/h;   
yinterp = y(:,ones(size(tinterp))) + LL + f*(h*BI)*cumprod([s;s;s;s]);

ypinterp = [];  
if nargout > 1
  ypinterp = f*BI*[ ones(size(s)); cumprod([2*s;3/2*s;4/3*s])];
end

% Non-negative solution
if ~isempty(idxNonNegative)
  idx = find(yinterp(idxNonNegative,:)<0); % vectorized
  if ~isempty(idx)
    w = yinterp(idxNonNegative,:);
    w(idx) = 0;
    yinterp(idxNonNegative,:) = w;
    if nargout > 1   % the derivative
      w = ypinterp(idxNonNegative,:);
      w(idx) = 0;
      ypinterp(idxNonNegative,:) = w;
    end      
  end
end  


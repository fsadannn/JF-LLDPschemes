function dxdt = f_dnd( t, y, theta )
%F_DND Summary of this function goes here
%   Detailed explanation goes here
N = length(y);

n1 = (N+1)*(N+1)/73;

dxdt = zeros(N,1);

dxdt(1) = ((y(2)-1)^2*0.25 + y(1)*(1-2.0*y(1)+y(2)))*n1 + y(1)*(1.0-y(1));
tt = y(2:N-1);
t2 = y(3:N);
dxdt(2:N-1) = ((t2-y(1:N-2)).^2.*0.25 + tt.*(y(1:N-2) - ...
    2.0.*tt+t2)).*n1 + tt.*(1.0-tt);

dxdt(N) = (y(N-1)*y(N-1)*0.25 + y(N)*(y(N-1)-2.0*y(N)))*n1 + y(N)*(1.0-y(N));
end

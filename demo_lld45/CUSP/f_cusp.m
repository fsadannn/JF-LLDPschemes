function ddt = f_cusp( t, y, theta )
%F_CUSP Summary of this function goes here
%   Detailed explanation goes here
N = length(y)/3;
epsc=1.0e-4;
epsc1 = 1/(epsc);
sigmac = 1/144;

% if nargin == 2
%     theta = 0;
% end;

% n1 = N;
sn2 = sigmac*N*N;
ddt = zeros(3*N,1);

ddt(1) = -epsc1*(y(1)^3+y(N+1)*y(1)+y(2*N+1)) + sn2*(y(N)-2.0*y(1)+y(2));
u = (y(1)-0.7)*(y(1)-1.3);
v = u/(u+0.1);
ddt(N+1) = y(2*N+1)+ 0.07*v + sn2*(y(2*N)-2*y(N+1)+y(N+2));
ddt(2*N+1) = (1-y(N+1)^2)*y(2*N+1)-y(N+1)-0.4*y(1) + 0.035*v + ...
            sn2*(y(3*N)-2.0*y(2*N+1)+y(2*N+2));

index = 2:N-1;    
ddt(index) = -epsc1.*(y(index).^3+y(N+2:2*N-1).*y(index)+y(2*N+2:3*N-1)) + ...
        sn2.*(y(1:N-2)-2.0.*y(index)+y(3:N));
u = (y(index)-0.7).*(y(index)-1.3);
v = u./(u+0.1);
ddt(N+2:2*N-1) = y(2*N+2:3*N-1)+ 0.07.*v + ...
            sn2.*(y(N+1:2*N-2)-2.*y(N+2:2*N-1)+y(N+3:2*N));
ddt(2*N+2:3*N-1) = (1-y(N+2:2*N-1).^2).*y(2*N+2:3*N-1)-y(N+2:2*N-1)-0.4.*y(2:N-1) + ...
            0.035.*v +  sn2.*(y(2*N+1:3*N-2)-2.0.*y(2*N+2:3*N-1)+y(2*N+3:3*N));

ddt(N) = -epsc1*(y(N)^3+y(2*N)*y(N)+y(3*N)) + sn2*(y(N-1)-2*y(N)+y(1));
u = (y(N)-0.7)*(y(N)-1.3);
v = u/(u+0.1);
ddt(2*N) = y(3*N)+ 0.07*v + sn2*(y(2*N-1)-2*y(2*N)+y(N+1));
ddt(3*N) = (1-y(2*N)^2)*y(3*N)-y(2*N)-0.4*y(N) + 0.035*v + ...
            sn2*(y(3*N-1)-2*y(3*N)+y(2*N+1));
        
end


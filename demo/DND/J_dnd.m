function dfdt = J_dnd( t, y, theta )
%J_DND Summary of this function goes here
%   Detailed explanation goes here

N = length(y);
persistent dxdt;

% if nargin == 2
%     theta = 0;
% end

n1 = (N+1)*(N+1)/73;

if length(dxdt)~= N
    dxdt = zeros(N,1);
end

[t1,t2] = size(y);

if t1<t2
    y = y(:);
end

yy=n1.*y;
dxdt(1) = 1.0-4.0*yy(1) + yy(2);
dxdt(2:N-1) = yy(1:N-2)-4.0.*yy(2:N-1)+yy(3:N);
dxdt(N) = yy(N-1)-4.0*yy(N);
dfdt = spdiags([0.5.*yy dxdt+1.0-2.0.*y 1.5.*yy],-1:1, N, N);
end


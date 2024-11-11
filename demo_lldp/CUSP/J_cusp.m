function dfdt = J_cusp( t, y, theta )
%J_CUSP Summary of this function goes here
%   Detailed explanation goes here
N = length(y)/3;
epsc=1.0e-4;
epsc1 = 1/(epsc);
sigmac = 1/144;
if nargin == 2
    theta = 0;
end
persistent ee;
%persistent K;
%persistent D1;
%persistent D2;
persistent D3;
%persistent D4;
persistent D5;
persistent D6;
%persistent D7;
%persistent D8;
%persistent D9;
%persistent SP;
%persistent cKKK;

[t1,t2] = size(y);
if t1<t2
    y = y(:);
end

sn2 = sigmac*N*N;

if length(ee)~=N
    ee = ones(N,1);
    D3 = spdiags(-epsc1.*ee, 0, N, N);
    D6 = speye(N,N);
    ee = sn2.*ee;
    D5 = spdiags([ee -2.0.*ee ee], -1:1, N, N);
    D5(1,N-2)=sn2;
    D5(N,3)=sn2;
end

tt = (y(1:N)-1.0)./((y(1:N).*2-2.0.*y(1:N)+1.01).^2);

D1 = spdiags([ee -2.0.*ee-epsc1.*(3.0.*y(1:N).*2+y(N+1:2*N)) ee], -1:1, N, N);
D2 = spdiags(-epsc1.*y(1:N), 0, N, N);
D4 = spdiags(0.014.*tt, 0, N, N);
D7 = spdiags(-0.4+0.007.*tt, 0, N, N);
D8 = spdiags(-1.0-2.0.*y(2*N+1:3*N).*y(N+1:2*N), 0, N, N);
D9 = spdiags([ee -2.0.*ee+1.0-y(N+1:2*N).*2 ee], -1:1, N, N);

D1(1,N-2)=sn2;
D1(N,3)=sn2;
D9(1,N-2)=sn2;
D9(N,3)=sn2;
dfdt = [D1,D2,D3;D4,D5,D6;D7,D8,D9];

end


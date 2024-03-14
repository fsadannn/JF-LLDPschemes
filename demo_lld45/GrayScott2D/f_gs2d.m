function [ dy ] = f_gs2d( t, y, theta )
%F_GS2D Summary of this function goes here
%   Detailed explanation goes here
Nd2 = length(y)/2;
N = sqrt(Nd2);
gsdu=0.2;
gsdv=0.1;
gsa=0.04;
gsb=0.06;
% persistent GS_I1;
% persistent GS_I2;
persistent GS_JAC;
persistent N2;

% if length(GS_I1)~= Nd2
if isempty(N2) || N2~= Nd2
    N2 = N*N;
%     GS_I1 = 1:N2;
%     GS_I2 = N2+1:2*N2;
    e=-ones(N,1);
    T = spdiags([e,-2.0.*e,e],-1:1,N,N);
    T(1,2) = -2.0;  T(N,N-1) = -2.0;
    c = (N-1)^2;
    GS_JAC = c*(kron(T,speye(N,N)) + kron(speye(N,N),T));
end

% tmp = (y(GS_I2).^2).*y(GS_I1);

tmp = (y(N2+1:2*N2).^2).*y(1:N2);

dy = [
    (-gsdu).*(GS_JAC*y(1:N2))+gsa.*(1-y(1:N2))-tmp;
        (-gsdv).*(GS_JAC*y(N2+1:2*N2))+ tmp - (gsa+gsb).*y(N2+1:2*N2)
    ];



% dy = [
%     (-gsdu).*(GS_JAC*y(GS_I1))+gsa.*(1-y(GS_I1))-tmp;
%         (-gsdv).*(GS_JAC*y(GS_I2))+ tmp - (gsa+gsb).*y(GS_I2)
%     ];


% dy = zeros(2*N2,1);  
% dy(GS_I1) = -gsdu.*(GS_JAC*y(GS_I1));
% dy(GS_I2) = -gsdv.*(GS_JAC*y(GS_I2));
% tmp = (y(GS_I2).^2).*y(GS_I1);
% dy(GS_I1) = dy(GS_I1) + gsa.*(1-y(GS_I1)) - tmp;
% dy(GS_I2) = dy(GS_I2) + tmp - (gsa+gsb).*y(GS_I2);

end
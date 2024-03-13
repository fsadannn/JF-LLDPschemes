%  phi1LLDP computes an approximation of phi(A)*u, phi(1/5*A)*u
%  phi(3/10*A)*u, phi(4/5*A)*u, phi(8/9*A)*u
%  for a general matrix A using Krylov subspace projection techniques.
%  Here, phi(z) = (exp(z)-I)/z and this phis are use in LLDP scheme.
%  It does not compute the matrix functions in isolation but instead,
%  it computes directly the action of these functions on the
%  operand vectors. This way of doing so allows for addressing large
%  sparse problems. The matrix under consideration interacts only
%  via matrix-vector products (matrix-free method).
%  It only compute one matrix exponential and use the relation betwin
%  the coefficients and the exponential properties to compute the
%  others phi.
function phi = phi1LLDP_single(A, u, h, m)

n=length(u);

btol  = 2*eps;
mb    = m;
beta = norm(u);
k1 = 3;
% begin Arnoldi
V = zeros(n,m+1);
H = zeros(m+3,m+3);
V(:,1) = (1/beta).*u;
if any(isnan(V(:,1)))
   phi = sparse(length(u),1);
   return;
end
A=h.*A;
for j = 1:m
    p = A*V(:,j);
    s = V(:,1:j);
    H(1:j,j) = s.'*p;
    p = p - s*H(1:j,j);
    s = norm(p);
    if s < btol
        k1 = 0;
        mb = j;
        break;
    end
    H(j+1,j) = s;
    V(:,j+1) = (1/s).*p;
end
H = (1/h).*H;
% build \hat{H}
if k1 == 0
    if mb>1
        mb=mb-1;
    end
    m=mb;
    hk = H(m+1,m);
    H=[H(1:m,1:m),zeros(m,3);zeros(3,m+3)];
    V=V(:,1:m+1);
else
    hk = H(m+1,m);
    H(m+1,m) = 0;
end
H(1,m+1) = 1;
H(m+1,m+2) = 1; H(m+2,m+3) = 1;

nhC= h*norm(H,'inf');
% scaling calculation
[~,e] = log2(nhC);
s = max(0,e+1);
% exponential calculation
pd=6;

M1 = expm64v41(h.*H,pd,s);

% calating \hat{E}
M1(m+1,m+1) = hk*M1(m,m+2);
M1(m+2,m+1) = hk*M1(m,m+3);

mx = m + 1;
% Matix projection

phi = V*(beta*M1(1:mx,mx));

end
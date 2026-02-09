function [phi, err, m, nexpo, breakdown, padepq] = phi1LLDP_hJ_f_gamma(A, b, h, hmin, y, hnormA, m, ...
                                                 rtol, atol,kdmax,kdmin, ...
                                                 gamma, reuse)
%PHI1LLDP_hJ_f Summary of this function goes here
%   Detailed explanation goes here
persistent V
persistent H
persistent mold
persistent hm1m

n=length(b);
m=min(m,n);
if isempty(mold)
    mold = m;
end
if nargin==12
    reuse=0;
end

% the coeficient is in rfacs( m )
% but is expensive call a function
% and coeficient are hardcode
rfacmin=1;
fac=1/log(2);

cte = 1./90;

nexpo=0;
padepq = 0;

y=abs(y);
tolr=[1.0e-9 1.0e-12 1.0e-15];
Table=[3 4 4;4 5 6];

breakdown=0;
btol = 2*eps;
mb = m;
rndoff= hnormA*eps;
k1 = 3;

% begin Arnoldi
if size(V,1)~=n
    hm1m=0;
    V = zeros(n,m+1);
    H = zeros(m+3,m+3);
else
    if mold>m
        H = [H(1:m,1:m),zeros(m,3);zeros(3,m+3)];
    elseif  m>mold
        H = [H(1:mold,1:mold),zeros(mold,m-mold+3);zeros(m-mold+3,m+3)]; 
    end
    if size(V,2)<m+1
        V = [V,zeros(n,m-size(V,2)+1)];
    end
end

reused = 0;
A=h.*A;
if reuse~=0 && mold>0
    if m==mold
        reused = 1; 
    elseif m>mold && size(V,1)==n
        reused = 2;
    end
    H(mold+1,mold)=hm1m;
    H = h.*H;
end

normb = norm(b);
if reused==0
   % build base of hA and b
   % begin Arnoldi
   V(:,1)=(1/normb).*b;
   % test if norm(b) is 0
   if any(isnan(V(:,1)))
        phi = sparse(n,5);
        err = 0;
        mold = 0;
        hm1m = 0;
        return;
   end
   for j = 1:m
        p = A*V(:,j);
        s = V(:,1:j);
        H(1:j,j) = s.'*p;
        p = p - s*H(1:j,j);
        s = norm(p);
        if s < btol
            k1 = 0;
            mb = j;
            breakdown=1;
            break
        end
        H(j+1,j) = s;
        V(:,j+1) = (1/s).*p;
   end
   % end Arnoldi
elseif reused==2
    % expand Krilov base of hA and b
    % begin Arnoldi
    mtemp=mold;
    mb=m;
    k1 = 3;
    for j = mtemp+1:m
        p = A*V(:,j);
        s = V(:,1:j);
        H(1:j,j) = s.'*p;
        p = p - s*H(1:j,j);
        s = norm(p);
        if s < btol
            k1 = 0;
            mb = j;
            breakdown=1;
            break;
        end
        H(j+1,j) = s;
        V(:,j+1) = (1/s).*p;
    end
    % end Arnoldi
end
% using scaling invariance property of Arnoldi
% rescaling H to have H for the subspace of A and b
H = (1/h).*H;
% build \hat{H}
if k1 == 0
    if mb>1
        mb=mb-1;
    end
    m=mb;
    hk = H(m+1,m);
    H=[H(1:m,1:m),zeros(m,3);zeros(3,m+3)];
else
    hk = H(m+1,m);
    H(m+1,m) = 0;
end
H(1,m+1) = 1;
H(m+1,m+2) = 1; H(m+2,m+3) = 1;
avm1dot=A*V(:,m+1);

work=1;

while work
    % select p-p of Pade
    nhC = h*cte*norm(H,'inf');
    col = find(tolr>=rtol,1,'last');
    if isempty(col)
        pd = 3;
    else
        fil = (nhC>=1)+1;
        pd = Table(fil,col);
    end
    padepq = max(padepq,pd);

    % scaling calculation
    [~,e] = log2(nhC);
    s = max(0,e+1);

    % exponential calculation
    M190 = expm64v41((h*cte).*H,pd,s);
    nexpo = nexpo + 1;

    M1 = M190*M190;     % M290
    M15 = M1*M1;        % m490
    M310 = M15*M15;     % M890
    M45 = M310*M310;    % M1690
    M1 = M45*M45;       % M3290
    M8090 = M1*M45*M1;  % M8090
    M45 = M310*M190;    % M110
    M15 = M45*M45;      % M15
    M310 = M45*M15;     % M310
    M1 = M15*M15;       % M25
    M45 = M1*M1;        % M45
    M1 = M45*M15;       % M1

    % calating \hat{E}
    beta = normb;
    %error relative
    % the divsion by h is because Av_{m+1} is in reality
    % hAv{m+1}(avm1dot) and need rescaling
%     err_rel=norm(((hk*M1(m,m+3)*beta/h).*avm1dot)./(atol+rtol.*y));
    err_rel=sqrt((1/n)*sum((((hk*M1(m,m+3)*beta/h).*avm1dot)./(atol+rtol.*y)).^2));
    if err_rel/gamma<1
       break
    else
        if breakdown==1 || m==kdmax
            hnew = max(hmin,0.9*h);
            if abs(h-hnew)<2*eps
                break
            end
            h = hnew;
        else
            % the coeficient is in rfacs( m )
            % but is expensive call a function
            % and coeficient are hardcode
            %  [ rfacmin,rfacmax ] = rfacs( m );
            rfacmax=max(1,m/3);
            knew =  log(err_rel/gamma)*fac;
            knew = ceil(m + min(rfacmax,max(knew,rfacmin)));
            knew = max(kdmin, min(kdmax,knew));

            H = [h.*H(1:m,1:m),zeros(m,knew-m+3);zeros(knew-m+3,knew+3)];
            H(m+1,m) = hk.*h;
            if size(V,2)<knew+1
                V = [V,zeros(n,knew-size(V,2)+1)];
            end
            mtemp=m+1;
            m=knew;
            mb=m;
            k1 = 3;
            
            j=mtemp;
            s = V(:,1:j);
            H(1:j,j) = s.'*avm1dot;
            avm1dot = avm1dot - s*H(1:j,j);
            s = norm(avm1dot);
            if s < btol
                k1 = 0;
                mb = j;
                breakdown=1;
                mtemp=m;
            end
            H(j+1,j) = s;
            V(:,j+1) = (1/s)*avm1dot;

            for j = mtemp+1:m
                p = A*V(:,j);
                s = V(:,1:j);
                H(1:j,j) = s.'*p;
                p = p - s*H(1:j,j);
                s = norm(p);
                if s < btol
                    k1 = 0;
                    mb = j;
                    breakdown=1;
                    break;
                end
                H(j+1,j) = s;
                V(:,j+1) = (1/s)*p;
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
            else
                hk = H(m+1,m);
                H(m+1,m) = 0;
            end
            H(1,m+1) = 1;
            H(m+1,m+2) = 1; H(m+2,m+3) = 1;
            avm1dot=A*V(:,m+1);

        end
    end

end

M1(m+1,m+1) = hk*M1(m,m+2);
M1(m+2,m+1) = hk*M1(m,m+3);
M15(m+1,m+1) = hk*M15(m,m+2);
M15(m+2,m+1) = hk*M15(m,m+3);
M310(m+1,m+1) = hk*M310(m,m+2);
M310(m+2,m+1) = hk*M310(m,m+3);
M45(m+1,m+1) = hk*M45(m,m+2);
M45(m+2,m+1) = hk*M45(m,m+3);
M8090(m+1,m+1) = hk*M8090(m,m+2);
M8090(m+2,m+1) = hk*M8090(m,m+3);

phi = zeros(n,5);
mx = m + 1;

avm1dot = V(:,1:mx);
% Matix projection
phi(:,1) = avm1dot*(beta*M15(1:mx,mx));
phi(:,2) = avm1dot*(beta*M310(1:mx,mx));
phi(:,3) = avm1dot*(beta*M45(1:mx,mx));
phi(:,4) = avm1dot*(beta*M8090(1:mx,mx));
phi(:,5) = avm1dot*(beta*M1(1:mx,mx));

err = max(err_rel,rndoff);
mold = m;
hm1m = hk;

end


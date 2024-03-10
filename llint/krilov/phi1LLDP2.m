function [phi, err, m, nexpo,breakdown,nfeval,padepq, korder] = phi1LLDP2(func, b, h, hmin, t, y, normy, m, ...
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
if nargin==13
    reuse=0;
end

jorder = 1;

% the coeficient is in rfacs( m )
% but is expensive call a function
% and coeficient are hardcode
rfacmin=1;
fac=1/log(2);

cte = 1./90;

nexpo=0;
nfeval=0;
padepq = 0;

%y=abs(y);
tolr=[1.0e-9 1.0e-12 1.0e-15];
Table=[3 4 4;4 5 6];

breakdown=0;
btol = 2*eps;
mb = m;
rndoff= eps;
k1 = 3;
ww = atol+rtol.*abs(y);

% hsubspace = 1;
delta = sqrt((1+normy)*eps);
idelta = 1/delta;
korder = abs(log(delta)/log(h))+1;

% if m == kdmin
%     kdmin = max([floor(korder)-1,1]);
%     m = kdmin;
% end

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
if reuse~=0 && mold>0
    if m==mold
        reused = 1;
    elseif m>mold && size(V,1)==n
        reused = 2;
    end
    H(mold+1,mold)=hm1m;
%     H = hsubspace.*H;
end

% if jorder==1
%     % order 1
%     hdelta = hsubspace/delta;
%     freejf = @(vector) FreeJ_f_w(1,func,y,vector,t,delta,hdelta,b);
%     % order 1
%     deltah = delta;
%     idelta = hsubspace/delta;
%     freejf = @(vector) idelta.*(func(t,y+deltah.*(vector))-b);
%     freejf = @(vector) (func(t,y+(delta*hsubspace).*(vector))-func(t,y))/(delta);
% else
%     % order 2
%     h2delta = h/(2*delta);
%     freejf = @(vector)FreeJ_f_w(2,func,y,vector,t,delta,h2delta);
% end

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
        %p = A*V(:,j);
%         if jorder==1
%             % order 1
%             p = FreeJ_f_w(1,func,y,V(:,j),t,delta,hdelta,b);
%         else
%             % order 2
%             p = FreeJ_f_w(2,func,y,V(:,j),t,delta,h2delta);
%         end
%         p = freejf(V(:,j));
        p = idelta.*(func(t,y+delta.*V(:,j))-b);
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
    nfeval=jorder*mb;
    % end Arnoldi
elseif reused==2
    % expand Krilov base of hA and b
    % begin Arnoldi
    mtemp=mold;
    mb=m;
    k1 = 3;
    for j = mtemp+1:m
        %p = A*V(:,j);
%         if jorder==1
%             % order 1
%             p = FreeJ_f_w(1,func,y,V(:,j),t,delta,hdelta,b);
%         else
%             % order 2
%             p = FreeJ_f_w(2,func,y,V(:,j),t,delta,h2delta);
%         end
%         p = freejf(V(:,j));
        p = idelta.*(func(t,y+delta.*V(:,j))-b);
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
    nfeval=jorder*abs(mb-mtemp+1);
    % end Arnoldi
end
% using scaling invariance property of Arnoldi
% rescaling H to have H for the subspace of A and b
% H = (1/hsubspace).*H;
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
% avm1dot=A*V(:,m+1);
% if jorder==1
%     % order 1
%     avm1dot = FreeJ_f_w(1,func,y,V(:,m+1),t,delta,hdelta,b);
% else
%     % order 2
%     avm1dot = FreeJ_f_w(2,func,y,V(:,m+1),t,delta,h2delta);
% end
% avm1dot = freejf(V(:,m+1));
avm1dot = idelta.*(func(t,y+delta.*V(:,m+1))-b);
nfeval=nfeval+jorder;

work=1;

while work
     % select p-p of Pade
    normH = norm(H,'inf');
    nhC = h*normH;
    
    if nhC > 600
        hnew = max(600/normH,hmin);
        h= hnew;
        nhC = h*normH;
    end
    
    nhC = cte*nhC;
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
%     err_rel=norm(((hk*M1(m,m+3)*beta/h).*avm1dot)./(atol+rtol.*abs(y)));
%     err_rel=sqrt((1/n)*sum((((hk*M1(m,m+3)*beta).*avm1dot)./ww).^2));
    phi_n=3;
    C=sqrt((1/n)*sum((((hk*beta).*avm1dot)./ww).^2));
%     C = max(((hk*beta).*avm1dot)./ww);
    err_rel=abs(C*M1(m,m+3));
    
    p1 = abs(M1(m,m+2));
    p2 = abs(M1(m,m+3))*norm(avm1dot);

    
    if p1<p2
%         disp(p1/p2);
%         err_rel=sqrt((1/n)*sum((((hk*M1(m,m+2)*beta).*V(:,m+1))./ww).^2));
         phi_n=2;
        C = sqrt((1/n)*sum((((hk*beta).*V(:,m+1))./ww).^2));
        err_rel=abs(C*M1(m,m+2));
    end
    
    
    if err_rel/gamma<1
        break;
    else
        if breakdown==1 || m==kdmax
%              hnew = max(hmin,0.9*h);
            HH = (H(1:m,1:m))^(m-1);
            hhh = HH(m,1); % e^T_m*H_m^(m-1)*e1
            lhh = log(hhh);
            ff = sum(log(2:(phi_n+m-1))); % log((r+n-1)!)
            
            ht = abs((log(gamma)-log(C)-lhh+ff)/(phi_n+m-1));
            h_new = exp(ht);
            
            hnew = max(hmin,min(0.9*h,h_new));  
            
%             disp([t,h,hnew,h_new]);
            
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
            
            H = [H(1:m,1:m),zeros(m,knew-m+3);zeros(knew-m+3,knew+3)];
            H(m+1,m) = hk;
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
                %p = A*V(:,j);
%                 if jorder==1
%                     % order 1
%                     p = FreeJ_f_w(1,func,y,V(:,j),t,delta,hdelta,b);
%                 else
%                     % order 2
%                     p = FreeJ_f_w(2,func,y,V(:,j),t,delta,h2delta);
%                 end
%                 p = freejf(V(:,j));
                p = idelta.*(func(t,y+delta.*V(:,j))-b);
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
            nfeval=nfeval+jorder*abs(mb-mtemp+1);
%             H = (1/hsubspace).*H;
            % build \hat{H}
            if k1 == 0
                if mb>1
                    mb=mb-1;
                end
                m=mb;
                hk = H(m+1,m);
                H=[H(1:m,1:m),zeros(m,3);zeros(3,m+3)];
                % V=V(:,1:m+1);
            else
                hk = H(m+1,m);
                H(m+1,m) = 0;
            end
            H(1,m+1) = 1;
            H(m+1,m+2) = 1; H(m+2,m+3) = 1;
            %avm1dot=A*V(:,m+1);
%             if jorder==1
%                 % order 1
%                 avm1dot = FreeJ_f_w(1,func,y,V(:,m+1),t,delta,hdelta,b);
%             else
%                 % order 2
%                 avm1dot = FreeJ_f_w(2,func,y,V(:,m+1),t,delta,h2delta);
%             end
%             avm1dot = freejf(V(:,m+1));
            avm1dot = idelta.*(func(t,y+delta.*V(:,m+1))-b);
            nfeval=nfeval+jorder;
        end
    end
end

phi = zeros(n,5);
mx = m + 1;

M1(mx,mx) = hk*M1(m,m+2);
M1(m+2,mx) = hk*M1(m,m+3);
M15(mx,mx) = hk*M15(m,m+2);
M15(m+2,mx) = hk*M15(m,m+3);
M310(mx,mx) = hk*M310(m,m+2);
M310(m+2,mx) = hk*M310(m,m+3);
M45(mx,mx) = hk*M45(m,m+2);
M45(m+2,mx) = hk*M45(m,m+3);
M8090(mx,mx) = hk*M8090(m,m+2);
M8090(m+2,mx) = hk*M8090(m,m+3);

if p1<p2
    avm1dot2 = V(:,1:m);
    % Matix projection
    phi(:,1) = avm1dot2*(beta*M15(1:m,mx)) + avm1dot.*(beta*M15(mx+1,mx));
    phi(:,2) = avm1dot2*(beta*M310(1:m,mx))+ avm1dot.*(beta*M310(mx+1,mx));
    phi(:,3) = avm1dot2*(beta*M45(1:m,mx))+ avm1dot.*(beta*M45(mx+1,mx));
    phi(:,4) = avm1dot2*(beta*M8090(1:m,mx))+ avm1dot.*(beta*M8090(mx+1,mx));
    phi(:,5) = avm1dot2*(beta*M1(1:m,mx))+ avm1dot.*(beta*M1(mx+1,mx));
else
    avm1dot = V(:,1:mx);
    % Matix projection
    phi(:,1) = avm1dot*(beta*M15(1:mx,mx));
    phi(:,2) = avm1dot*(beta*M310(1:mx,mx));
    phi(:,3) = avm1dot*(beta*M45(1:mx,mx));
    phi(:,4) = avm1dot*(beta*M8090(1:mx,mx));
    phi(:,5) = avm1dot*(beta*M1(1:mx,mx));
end

err = max(err_rel,rndoff);
mold = m;
hm1m = hk;
end


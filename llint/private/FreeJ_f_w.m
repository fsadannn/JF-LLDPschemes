function [ hJv ]  = FreeJ_f_w( order, f, u, v, t, delta, hdelta, fu )
% FREEJ compute h*J(u)*v where J is the Jacobian evaluate at u
% this is the same as the directional derivate and this product
% is compute by difference quotients
% order is 1 or 2
% if order is 1 the formula h*(f(u+delta*v)-f(u))/delta is used
% if order is 2 the formula h*(f(u+delta*v)-f(u-delta*v))/delta is used
% f is the function handle
% u is the vector where the Jacobian is evaluate
% v is the vector that multiplied by Jacobian
% fu is the evaluation of the function in the vector u and this arg is
% optional
% delta is the factaror (f(u+delta*v)-f(u))/delta
% hdelta is the cocienta betwin h and delta
efu = 1;
if order == 1
    if nargin == 7
        efu = 0;
    end
    if efu==0
            f1 = feval(f,t,u + delta.*v);
            f2 = feval(f,t,u);
        hJv = hdelta.*(f1-f2);
    else
        hJv = hdelta.*(feval(f,t,u + delta.*v) - fu);
    end
    
elseif order == 2
    deltav=delta.*v;
    f1 = f(t,u + deltav);
    f2 = f(t,u - deltav);
    if nargin == 6
        hJv = (1/(2*delta)).*(f1-f2);
    else
        hJv = hdelta.*(f1-f2);
    end
     
end

end


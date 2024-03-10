fname = @f_dnd;
Jname = @J_dnd;
N=500;
ptsx=50.*(1:N)/(N+1)+(1-(1:N)/(N+1)).*(-23);
sgx=sign(ptsx);
x0=((sgx+1)/2).*exp(-0.8284.*ptsx)+(1-sgx)/2;
clear ptsx;
clear sgx;
IT=[0 0.1];

options15sExact=odeset('RelTol',1.0e-12,'AbsTol',1.0e-14,'Jacobian',Jname);

Tol={"Refined";"Refined";"Mild";"Mild";"Crude";"Crude"};
Name = {"LLDP1";"LLDP2";"LLDP1";"LLDP2";"LLDP1";"LLDP2"};
RError = [0;0;0;0;0;0];
ASteps = [0;0;0;0;0;0];
RSteps = [0;0;0;0;0;0];
fEval = [0;0;0;0;0;0];
KSubspace = [0;0;0;0;0;0];
ME = [0;0;0;0;0;0];
mtotal = [0;0;0;0;0;0];
mmin = [0;0;0;0;0;0];
mmax = [0;0;0;0;0;0];
pmin = [0;0;0;0;0;0];
pmax = [0;0;0;0;0;0];

% Refined
ATol = 1.0e-12;
RTol = 1.0e-9;
optionsLL  = llset('RelTol',RTol,'AbsTol',ATol,'dKmin',4,'dKmax',50,'debug',0);

SolLL = LLDP1(fname,IT,x0,optionsLL);
TLL = SolLL.x;
YLL = real(SolLL.y)';
[~,Y] = ode15s(fname,TLL,x0,options15sExact);
LLRE = RelError(Y,YLL);
RError(1)=LLRE;
ASteps(1)=SolLL.stats.nsteps;
RSteps(1)=SolLL.stats.nfailed;
fEval(1)=SolLL.stats.nfevals;
KSubspace(1)=SolLL.stats.nsteps;
ME(1)=SolLL.stats.nexpm;
mtotal(1)=SolLL.stats.Kdim_sum;
mmin(1)=SolLL.stats.Kdim_min;
mmax(1)=SolLL.stats.Kdim_max;
pmin(1)=SolLL.stats.pade_min;
pmax(1)=SolLL.stats.pade_max;

SolLL = LLDP2(fname,IT,x0,optionsLL);
TLL = SolLL.x;
YLL = real(SolLL.y)';
[~,Y] = ode15s(fname,TLL,x0,options15sExact);
LLRE = RelError(Y,YLL);
RError(2)=LLRE;
ASteps(2)=SolLL.stats.nsteps;
RSteps(2)=SolLL.stats.nfailed;
fEval(2)=SolLL.stats.nfevals;
KSubspace(2)=SolLL.stats.nsteps;
ME(2)=SolLL.stats.nexpm;
mtotal(2)=SolLL.stats.Kdim_sum;
mmin(2)=SolLL.stats.Kdim_min;
mmax(2)=SolLL.stats.Kdim_max;
pmin(2)=SolLL.stats.pade_min;
pmax(2)=SolLL.stats.pade_max;

% Mild
ATol = 1.0e-9;
RTol = 1.0e-6;
optionsLL  = llset('RelTol',RTol,'AbsTol',ATol,'dKmin',4,'dKmax',50,'debug',0);

SolLL = LLDP1(fname,IT,x0,optionsLL);
TLL = SolLL.x;
YLL = real(SolLL.y)';
[~,Y] = ode15s(fname,TLL,x0,options15sExact);
LLRE = RelError(Y,YLL);
RError(3)=LLRE;
ASteps(3)=SolLL.stats.nsteps;
RSteps(3)=SolLL.stats.nfailed;
fEval(3)=SolLL.stats.nfevals;
KSubspace(3)=SolLL.stats.nsteps;
ME(3)=SolLL.stats.nexpm;
mtotal(3)=SolLL.stats.Kdim_sum;
mmin(3)=SolLL.stats.Kdim_min;
mmax(3)=SolLL.stats.Kdim_max;
pmin(3)=SolLL.stats.pade_min;
pmax(3)=SolLL.stats.pade_max;

SolLL = LLDP2(fname,IT,x0,optionsLL);
TLL = SolLL.x;
YLL = real(SolLL.y)';
[~,Y] = ode15s(fname,TLL,x0,options15sExact);
LLRE = RelError(Y,YLL);
RError(4)=LLRE;
ASteps(4)=SolLL.stats.nsteps;
RSteps(4)=SolLL.stats.nfailed;
fEval(4)=SolLL.stats.nfevals;
KSubspace(4)=SolLL.stats.nsteps;
ME(4)=SolLL.stats.nexpm;
mtotal(4)=SolLL.stats.Kdim_sum;
mmin(4)=SolLL.stats.Kdim_min;
mmax(4)=SolLL.stats.Kdim_max;
pmin(4)=SolLL.stats.pade_min;
pmax(4)=SolLL.stats.pade_max;


% Crude
ATol = 1.0e-6;
RTol = 1.0e-3;
optionsLL  = llset('RelTol',RTol,'AbsTol',ATol,'dKmin',4,'dKmax',50,'debug',0);

SolLL = LLDP1(fname,IT,x0,optionsLL);
TLL = SolLL.x;
YLL = real(SolLL.y)';
[~,Y] = ode15s(fname,TLL,x0,options15sExact);
LLRE = RelError(Y,YLL);
RError(5)=LLRE;
ASteps(5)=SolLL.stats.nsteps;
RSteps(5)=SolLL.stats.nfailed;
fEval(5)=SolLL.stats.nfevals;
KSubspace(5)=SolLL.stats.nsteps;
ME(5)=SolLL.stats.nexpm;
mtotal(5)=SolLL.stats.Kdim_sum;
mmin(5)=SolLL.stats.Kdim_min;
mmax(5)=SolLL.stats.Kdim_max;
pmin(5)=SolLL.stats.pade_min;
pmax(5)=SolLL.stats.pade_max;

SolLL = LLDP2(fname,IT,x0,optionsLL);
TLL = SolLL.x;
YLL = real(SolLL.y)';
[~,Y] = ode15s(fname,TLL,x0,options15sExact);
LLRE = RelError(Y,YLL);
RError(6)=LLRE;
ASteps(6)=SolLL.stats.nsteps;
RSteps(6)=SolLL.stats.nfailed;
fEval(6)=SolLL.stats.nfevals;
KSubspace(6)=SolLL.stats.nsteps;
ME(5)=SolLL.stats.nexpm;
mtotal(6)=SolLL.stats.Kdim_sum;
mmin(6)=SolLL.stats.Kdim_min;
mmax(6)=SolLL.stats.Kdim_max;
pmin(6)=SolLL.stats.pade_min;
pmax(6)=SolLL.stats.pade_max;

% Table
disp("DND: M=500, d=500");
disp(table(Tol,Name,RError,ASteps,RSteps,fEval,KSubspace,ME,mtotal,mmin,mmax,pmin,pmax));
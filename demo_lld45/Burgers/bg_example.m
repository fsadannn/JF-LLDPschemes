fname=@f_burgers;
Jname=@J_burgers;
n=500;
N=n;
x0=((sin((3.0*pi/(N+1)).*(1:N))).^2).*((1.0-1/(N+1).*(1:N)).^(3/2));
IT=[0 0.5];

options15sExact=odeset('RelTol',1.0e-12,'AbsTol',1.0e-14,'Jacobian',Jname);

Tol={"Refined";"Mild";"Crude"};
Name = {"LLDP45";"LLDP45";"LLDP45"};
RError = [0;0;0];
ASteps = [0;0;0];
RSteps = [0;0;0];
fEval = [0;0;0];
JEval = [0;0;0];
KSubspace = [0;0;0];
ME = [0;0;0];
mtotal = [0;0;0];
mmin = [0;0;0];
mmax = [0;0;0];
pmin = [0;0;0];
pmax = [0;0;0];

% Refined
ATol = 1.0e-12;
RTol = 1.0e-9;
optionsLL  = llset('RelTol',RTol,'AbsTol',ATol,'dKmin',4,'dKmax',50);

SolLL = LLDP45(fname,Jname,IT,x0,optionsLL);
TLL = SolLL.x;
YLL = real(SolLL.y)';
[~,Y] = ode15s(fname,TLL,x0,options15sExact);
LLRE = RelError(Y,YLL);
RError(1)=LLRE;
ASteps(1)=SolLL.stats.nsteps;
RSteps(1)=SolLL.stats.nfailed;
fEval(1)=SolLL.stats.nfevals;
JEval(1)=SolLL.stats.nJevals;
KSubspace(1)=SolLL.stats.nsteps;
ME(1)=SolLL.stats.nexpm;
mtotal(1)=SolLL.stats.Kdim_sum;
mmin(1)=SolLL.stats.Kdim_min;
mmax(1)=SolLL.stats.Kdim_max;
pmin(1)=SolLL.stats.pade_min;
pmax(1)=SolLL.stats.pade_max;

% Mild
ATol = 1.0e-9;
RTol = 1.0e-6;
optionsLL  = llset('RelTol',RTol,'AbsTol',ATol,'dKmin',4,'dKmax',50,'debug',0);

SolLL = LLDP45(fname,Jname,IT,x0,optionsLL);
TLL = SolLL.x;
YLL = real(SolLL.y)';
[~,Y] = ode15s(fname,TLL,x0,options15sExact);
LLRE = RelError(Y,YLL);
RError(2)=LLRE;
ASteps(2)=SolLL.stats.nsteps;
RSteps(2)=SolLL.stats.nfailed;
fEval(2)=SolLL.stats.nfevals;
JEval(2)=SolLL.stats.nJevals;
KSubspace(2)=SolLL.stats.nsteps;
ME(2)=SolLL.stats.nexpm;
mtotal(2)=SolLL.stats.Kdim_sum;
mmin(2)=SolLL.stats.Kdim_min;
mmax(2)=SolLL.stats.Kdim_max;
pmin(2)=SolLL.stats.pade_min;
pmax(2)=SolLL.stats.pade_max;

% Crude
ATol = 1.0e-6;
RTol = 1.0e-3;
optionsLL  = llset('RelTol',RTol,'AbsTol',ATol,'dKmin',4,'dKmax',50,'debug',0);

SolLL = LLDP45(fname,Jname,IT,x0,optionsLL);
TLL = SolLL.x;
YLL = real(SolLL.y)';
[~,Y] = ode15s(fname,TLL,x0,options15sExact);
LLRE = RelError(Y,YLL);
RError(3)=LLRE;
ASteps(3)=SolLL.stats.nsteps;
RSteps(3)=SolLL.stats.nfailed;
fEval(3)=SolLL.stats.nfevals;
JEval(3)=SolLL.stats.nJevals;
KSubspace(3)=SolLL.stats.nsteps;
ME(3)=SolLL.stats.nexpm;
mtotal(3)=SolLL.stats.Kdim_sum;
mmin(3)=SolLL.stats.Kdim_min;
mmax(3)=SolLL.stats.Kdim_max;
pmin(3)=SolLL.stats.pade_min;
pmax(3)=SolLL.stats.pade_max;


% Table
disp("Burgers: M=500, d=500");
disp(table(Tol,Name,RError,ASteps,RSteps,fEval,JEval,KSubspace,ME,mtotal,mmin,mmax,pmin,pmax));
fname = @f_cusp;
Jname = @J_cusp;
N=500;
x0= [zeros(1,N),-2*cos((2*pi/(N))*(1:N)),2*sin((2*pi/(N))*(1:N))];
IT=[0 1.0e-4];

optionsLLexact = llset('RelTol',1e-12,'AbsTol',1.0e-14,'dKmin',100,...
        'dKmax',100,'jacgap', inf);

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
[~,Y] = LLDP45(fname,Jname,TLL,x0,optionsLLexact);
clear phi1LLDP_hJ_f_gamma;
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
[~,Y] = LLDP45(fname,Jname,TLL,x0,optionsLLexact);
clear phi1LLDP_hJ_f_gamma;
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
[~,Y] = LLDP45(fname,Jname,TLL,x0,optionsLLexact);
clear phi1LLDP_hJ_f_gamma;
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
disp("Cusp: M=500, d=1500");
disp(table(Tol,Name,RError,ASteps,RSteps,fEval,JEval,KSubspace,ME,mtotal,mmin,mmax,pmin,pmax));
clear phi1LLDP_hJ_f_gamma;
%This file solves and simulates the model developed in Eichenbaum,
%Rebelo and Trabandt (2020),'The Macroeconomics of Epidemics'.

%Matlab 2016b used for calculations.

clear all; clc; close all; tic;


%Scenario: optimal containment for 12 weeks, then 
%unanticipated cut of muc back to zero

%Implementation: take initial conditions from week 12 and solve and simulate  
%model forward with those initial conditions and zero containment policy  
load results_vaccine_treatment_medprep_optimal.mat

exitdate=45; %week when optimal containment will be stopped unexpectedly

%State variables in week of exit (initial conditions for 
%simulation with muc zero going forward)
global pop_ini_break i_ini_break s_ini_break d_ini_break r_ini_break

pop_ini_break=Pop(exitdate); 
i_ini_break=I(exitdate); %Note, I(1) is a state variable, i.e. I_0 for Cons. in period 1. That is, at exitdate, the entry I(exitdate) is the state variable at time exitdate
s_ini_break=S(exitdate);
d_ini_break=D(exitdate);
r_ini_break=R(exitdate); 


%save allocations with optimal containment policy 
muc_opt=muc;
I_opt=I;
S_opt=S;
R_opt=R;
D_opt=D;
T_opt=T;
Pop_opt=Pop;
cs_opt=cs;
ns_opt=ns;
Us_opt=Us;
RnotSIRmacro_opt=RnotSIRmacro;
aggC_opt=aggC;
aggH_opt=aggH;
ci_opt=ci;
cr_opt=cr;
ni_opt=ni;
nr_opt=nr;
Ui_opt=Ui;
Ur_opt=Ur;
U_opt=U;
pid_endo_opt=pid_endo;


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters, calibration targets and other settings%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betta=0.96^(1/52);  %Weekly household discount factor
pid=7*0.005/18;     %Weekly probability of dying
pir=7*1/18-pid;     %Weekly probability of recovering
phii=0.8;           %Productivity of infected people

deltav=1/52;        %Weekly probability of discovering a vaccine
deltac=1/52;        %Weekly probability of discovering a treatment
kappa=0.9;          %Slope of pid-function in endog. pid scenario 
                    %(medical preparedness scenario)

%Calibration targets for hours and income
n_target=28;         %Weekly hours
inc_target=58000/52; %weekly income

%Calibation targets for shares of pis-terms (pi terms in latest draft) in T-function in SIR model
pis3_shr_target=2/3;                   %share of T_0 jump due general infections
pis1_shr_target=(1-pis3_shr_target)/2; %share of T_0 jump due to consumption-based infections
pis2_shr_target=(1-pis3_shr_target)/2; %share of T_0 jump due to work-based infections
RplusD_target=0.60;                    %total share of people infected and then either recovered or dead after epidemic

pop_ini=1;          %Initial population
i_ini=0.001;        %Initial infected

HH=250;             %Number of periods to solve and simulate the model

%containment policy
muc=zeros(HH,1);    %exogenous path for muc (mu in latest draft) over time.
                    %if you want e.g. a containment policy of
                    %10 percent for weeks 1,...52, then set muc(1:52)=0.1;
                    %
                    %Make sure that muc=0 at the end of the solution and
                    %simulation horizon. Steady state assumes muc=0;
                    %
                    %With optimal policy, muc path will be chosen to
                    %maximize PV utility (switch below).                                           
                    
do_opt_policy=0;    %switch: if = 0, model is solved and simulated with
                    %given path for containment policy muc
                    %
                    %if =1, model is solved and simulated with optimal
                    %containment path muc. Path for muc set above will be
                    %overwritten with optimal path. COMPUTATIONS TAKE A
                    %WHILE IF YOU SET do_opt_taxes=1
                    
use_parallel=0;     %when optimal policy is computed, use_parallel=1 uses 
                    %parallel computing to maximize PV utility using fmincon.

%nonlinear solver and minimizer settings
opts_fsolve=optimoptions('fsolve','Display','iter','TolFun',1e-9); %options for fsolve
opts_fsolve_fmincon=optimoptions('fsolve','Display','off','TolFun',1e-9); %options for fsolve used opt. policy calcs. (fmincon)
if use_parallel==0
    opts_fmincon=optimoptions('fmincon','Display','iter','TolFun',1e-7,'MaxFunctionEvaluations',2000,'FiniteDifferenceStepSize',1e-2); %options for fmincon w/o parallel comp.
elseif use_parallel==1
    opts_fmincon=optimoptions('fmincon','Display','iter','TolFun',1e-7,'MaxFunctionEvaluations',10000,'UseParallel',true,'FiniteDifferenceStepSize',1e-3);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Steady State Calculations%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta=1/n_target^2;     %calculate disutility of labor parameter theta
                        %so that pre-infection steady state labor is equal to
                        %desired target (using n=(1/theta)^(1/2) pre-infection 
                        %steady state equation)
A=inc_target/n_target;  %calculate parameter A such that income is equal to
                        %desired target, c=inc_target=A*n_target

%steady states
nrss=(1/theta)^(1/2);           %labor recovered (same as post-infection steady state)
crss=A*nrss;                    %consumption recovered
urss=log(crss)-theta/2*nrss^2;  %utility recovered
Urss=1/(1-betta)*urss;          %PV utility recovered
UrssConsUnits=Urss*crss;        %PV utility in cons. units (Urss*Marg.Util.Cons); value of life
niss=(1/theta)^(1/2);           %labor infected
ciss=phii*A*niss;               %consumption infected
uiss=log(ciss)-theta/2*niss^2;  %utility infected
Uiss=1/(1-(1-deltac)*betta*(1-pir-pid))*(uiss...
+(1-deltac)*betta*pir*Urss+deltac*betta*Urss);  %PV utility infected


%Check level of present value utility
if Uiss-Urss>0, error(['Error: parameterization implies Uiss>Urss: ',num2str(Uiss-Urss)]);end

%calibrate the pis's in T-function
go_calibrate_pis;

%initial guess for optimal muc and load last solution of optimal policy 
%allocations if you dont want to start maximization of PV utility from scratch.
if do_opt_policy==1      
    muc_guess=zeros(HH,1);                  %initial guess for opt. cont. policy
    
    %load last_solution_opt_policy    %uncomment if you want to use last
                                      %solution as initial guess
    %muc_guess=muc_sol;
    %if numel(muc_guess)~=HH
    %    error('Initial guess for optimal policy loaded from disk has different dimension than HH.');
    %end
end

%initial guess of vectors of ns, ni and nr to solve nonlinear
%equilibrium equations
n_vec_guess=nrss*ones(3*HH,1); %guess of vectors for ns,ni,nr

%If optimal policy is desired, find optimal path for muc to maximize PV utility
if do_opt_policy==1 %optimal policy
    
    %minimize negative PV utility (i.e. max PV utility) to find opt. path for muc; nonlinear
    %model equations are solved inside getU.m using function get_err.m
    LB=muc_guess*0-2;UB=muc_guess*0+2;%lower and upper bounds
    muc_sol = fmincon(@getU,muc_guess,[],[],[],[],LB,UB,[],opts_fmincon,n_vec_guess,opts_fsolve_fmincon,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uiss,HH,crss,nrss,Urss,phii,deltav,deltac,kappa);
    
    muc=zeros(HH,1);muc(1:numel(muc_sol))=muc_sol;
    save last_solution_opt_policy muc_sol;%save solution for possible use in subsequent maximization   
end

%Given either optimal path for muc or exogenous path for muc,
%solve nonlinear equilibrium model equations (i.e. adjust guesses ns,nr,ni)
[n_vec,fval,exitflag]=fsolve(@get_err,n_vec_guess,opts_fsolve,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uiss,HH,crss,nrss,Urss,muc,phii,deltav,deltac,kappa);
if exitflag~=1
    error('Fsolve could not solve the model');
end    
    
%get allocations given either exogenous or optimal path for muc at ns,ni,nr
%solution
[err_exit,I_exit,S_exit,R_exit,D_exit,T_exit,Pop_exit,cs_exit,ns_exit,Us_exit,RnotSIRmacro_exit,aggC_exit,aggH_exit,ci_exit,cr_exit,ni_exit,nr_exit,Ui_exit,Ur_exit,U_exit,pid_endo_exit] = get_err(n_vec,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,Uiss,HH,crss,nrss,Urss,muc,phii,deltav,deltac,kappa);
disp(['Max. abs. error in equilib. equations:',num2str(max(abs(err)))]);
disp(' ');
RnotSIRmacro_exit
muc_exit=muc;

%Put together data from allocations of optimal containment for 12 weeks followed by 
%allocations with zero containment thereafter (exit)
I=[I_opt(1:exitdate-1,1);I_exit(1:HH-exitdate+1,1)];
S=[S_opt(1:exitdate-1,1);S_exit(1:HH-exitdate+1,1)];
R=[R_opt(1:exitdate-1,1);R_exit(1:HH-exitdate+1,1)];
D=[D_opt(1:exitdate-1,1);D_exit(1:HH-exitdate+1,1)];

muc=[muc_opt(1:exitdate-1,1);muc_exit(1:HH-exitdate+1,1)];
T=[T_opt(1:exitdate-1,1);T_exit(1:HH-exitdate+1,1)];
Pop=[Pop_opt(1:exitdate-1,1);Pop_exit(1:HH-exitdate+1,1)];
cs=[cs_opt(1:exitdate-1,1);cs_exit(1:HH-exitdate+1,1)];
ns=[ns_opt(1:exitdate-1,1);ns_exit(1:HH-exitdate+1,1)];
Us=[Us_opt(1:exitdate-1,1);Us_exit(1:HH-exitdate+1,1)];
aggC=[aggC_opt(1:exitdate-1,1);aggC_exit(1:HH-exitdate+1,1)];
aggH=[aggH_opt(1:exitdate-1,1);aggH_exit(1:HH-exitdate+1,1)];
ci=[ci_opt(1:exitdate-1,1);ci_exit(1:HH-exitdate+1,1)];
cr=[cr_opt(1:exitdate-1,1);cr_exit(1:HH-exitdate+1,1)];
ni=[ni_opt(1:exitdate-1,1);ni_exit(1:HH-exitdate+1,1)];
nr=[nr_opt(1:exitdate-1,1);nr_exit(1:HH-exitdate+1,1)];
Ui=[Ui_opt(1:exitdate-1,1);Ui_exit(1:HH-exitdate+1,1)];
Ur=[Ur_opt(1:exitdate-1,1);Ur_exit(1:HH-exitdate+1,1)];
U=[U_opt(1:exitdate-1,1);U_exit(1:HH-exitdate+1,1)];
pid_endo=[pid_endo_opt(1:exitdate-1,1);pid_endo_exit(1:HH-exitdate+1,1)];


%plotting
ia=2;ib=2;fsize=12;
horz=HH;
time=0:1:horz-1;

figure;
subplot(ia,ib,1)
plot(time,I(1:horz),'b-','LineWidth',2);hold on;
plot(time(exitdate),I(exitdate),'ro','MarkerSize',5,'MarkerFaceColor','r');
box off;
title('Infected, I','FontSize',fsize);
ylabel('Share of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,2)
plot(time,S(1:horz),'b-','LineWidth',2);hold on;
plot(time(exitdate),S(exitdate),'ro','MarkerSize',5,'MarkerFaceColor','r');
box off;
title('Susceptibles, S','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,3)
plot(time,R(1:horz),'b-','LineWidth',2);hold on;
plot(time(exitdate),R(exitdate),'ro','MarkerSize',5,'MarkerFaceColor','r');
box off;
title('Recovered, R','FontSize',fsize);
ylabel('Share of Initial Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,4)
plot(time,D(1:horz),'b-','LineWidth',2);hold on;
plot(time(exitdate),D(exitdate),'ro','MarkerSize',5,'MarkerFaceColor','r');
box off;
title('Deaths, D','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);

suptitle('The Evolution of an Epidemic');
orient landscape
print -dpdf -fillpage SIRmacro_epidemic_simulation_fig1


figure;
subplot(ia,ib,3)
plot(time,100*(cs(1:horz)-crss)/crss,'b-','LineWidth',2); hold on
plot(time,100*(ci(1:horz)-crss)/crss,'r--','LineWidth',2); hold on
plot(time,100*(cr(1:horz)-crss)/crss,'k-.','LineWidth',2);hold on
plot(time(exitdate),100*(cs(exitdate)-crss)/crss,'ro','MarkerSize',5,'MarkerFaceColor','r');hold on
plot(time(exitdate),100*(ci(exitdate)-crss)/crss,'ro','MarkerSize',5,'MarkerFaceColor','r');hold on
plot(time(exitdate),100*(cr(exitdate)-crss)/crss,'ro','MarkerSize',5,'MarkerFaceColor','r');
box off;
title('Consumption by Type','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Cons. Susceptibles','Cons. Infected','Cons. Recovered','Location','best');
legend boxoff

subplot(ia,ib,4)
plot(time,100*(ns(1:horz)-nrss)/nrss,'b-','LineWidth',2); hold on
plot(time,100*(ni(1:horz)-nrss)/nrss,'r--','LineWidth',2); hold on
plot(time,100*(nr(1:horz)-nrss)/nrss,'k-.','LineWidth',2); hold on
plot(time(exitdate),100*(ns(exitdate)-nrss)/nrss,'ro','MarkerSize',5,'MarkerFaceColor','r');hold on
plot(time(exitdate),100*(ni(exitdate)-nrss)/nrss,'ro','MarkerSize',5,'MarkerFaceColor','r');hold on
plot(time(exitdate),100*(nr(exitdate)-nrss)/nrss,'ro','MarkerSize',5,'MarkerFaceColor','r');
box off;
title('Hours by Type','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Hours Susceptibles','Hours Infected','Hours Recovered','Location','best');
legend boxoff

subplot(ia,ib,1)
plot(time,0*100*(aggC(1:horz)-crss)/crss,'m:','LineWidth',1.5);hold on
plot(time,100*(aggC(1:horz)-crss)/crss,'b-','LineWidth',2); hold on
plot(time(exitdate),100*(aggC(exitdate)-crss)/crss,'ro','MarkerSize',5,'MarkerFaceColor','r');
box off;
title('Aggregate Consumption, C','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,2)
plot(time,0*100*(aggC(1:horz)-crss)/crss,'m:','LineWidth',1.5);hold on
plot(time,100*(aggH(1:horz)-nrss)/nrss,'b-','LineWidth',2);hold on;
plot(time(exitdate),100*(aggH(exitdate)-nrss)/nrss,'ro','MarkerSize',5,'MarkerFaceColor','r');
box off;
title('Aggregate Hours, H','FontSize',fsize);
set(gca,'FontSize',fsize);

suptitle('The Evolution of an Epidemic');
orient landscape
print -dpdf -fillpage SIRmacro_epidemic_simulation_fig2



figure;
subplot(3,1,1)
plot(time,Us(1:HH),'b-','LineWidth',2); hold on;
plot(time,Ui(1:HH),'r--','LineWidth',2);hold on;
plot(time,Ur(1:HH),'k-.','LineWidth',2);hold on;
plot(time,U(1:HH),'m:','LineWidth',2);hold on;
plot(time(exitdate),Us(exitdate),'ro','MarkerSize',5,'MarkerFaceColor','r');
hold on;
plot(time(exitdate),Ui(exitdate),'ro','MarkerSize',5,'MarkerFaceColor','r');
hold on;
plot(time(exitdate),Ur(exitdate),'ro','MarkerSize',5,'MarkerFaceColor','r');
hold on;
plot(time(exitdate),U(exitdate),'ro','MarkerSize',5,'MarkerFaceColor','r');

set(gca,'FontSize',fsize);
legend('Susceptibles','Infected','Recovered','Total','Location','best');
set(gca,'FontSize',fsize);
title('Present Value Utility','FontSize',fsize);
legend boxoff;

subplot(3,1,2)
plot(time,0*muc,'k:','LineWidth',1.5); hold on
plot(time,muc,'b-','LineWidth',2);  hold on;
plot(time(exitdate),muc(exitdate),'ro','MarkerSize',5,'MarkerFaceColor','r');
set(gca,'FontSize',fsize);
title('Containment Policy, \mu_c','FontSize',fsize);
ylabel('%','FontSize',fsize)

subplot(3,1,3)
plot(time,100*pid*18/7,'k:','LineWidth',1.5); hold on
plot(time,100*pid_endo*18/7,'b-','LineWidth',2); hold on;
plot(time(exitdate),100*pid_endo(exitdate)*18/7,'ro','MarkerSize',5,'MarkerFaceColor','r');
set(gca,'FontSize',fsize);
title('Mortality Rate, \pi_d','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
ylabel('%','FontSize',fsize)

suptitle('The Evolution of an Epidemic');
orient portrait
print -dpdf -fillpage SIRmacro_epidemic_simulation_fig3


%Save workspace
save results_vaccine_treatment_medprep_optimal_earlyexit


%output some data used in paper and robustness table
aggCons_trough_percent=min((100*(aggC-crss)/crss))
aggCons_avg_first_year_percent=mean((100*(aggC(1:52)-crss)/crss))
terminal_one_minus_susceptibles_percent=100*(1-S(end))
peak_infection_percent=max(100*I)
terminal_death_share_percent=100*D(end)
terminal_number_deaths_US_millions=terminal_death_share_percent/100*330

toc;
%This file solves and simulates the model developed in Eichenbaum,
%Rebelo and Trabandt (2020),'The Macroeconomics of Epidemics'.

%Matlab 2016b used for calculations.

clear all; clc; close all; tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters, calibration targets and other settings%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
betta=0.96^(1/52);  %Weekly household discount factor
pid=7*0.005/18;     %Weekly probability of dying
pir=7*1/18-pid;     %Weekly probability of recovering
phii=0.8;           %Productivity of infected people

deltav=1/52;        %Weekly probability of discovering a vaccine
deltac=1/52;        %Weekly probability of discovering a treatment
kappa=0.9;          %Slope of pid-function in endog. pid scenario (medical preparedness scenario)

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

HH=500;             %Number of periods to solve and simulate the model


use_parallel=1;     %when optimal policy is computed, use_parallel=1 uses
                    %parallel computing to maximize PV utility using fmincon.

%nonlinear solver and minimizer settings
opts_fsolve=optimoptions('fsolve','Display','iter','TolFun',1e-9); %options for fsolve
opts_fsolve_fmincon=optimoptions('fsolve','Display','off','TolFun',1e-9); %options for fsolve used opt. policy calcs. (fmincon)
if use_parallel==0
    opts_fmincon=optimoptions('fmincon','Display','iter','TolFun',1e-7,'MaxFunctionEvaluations',10000,'FiniteDifferenceStepSize',1e-3); %options for fmincon w/o parallel comp.
elseif use_parallel==1
    opts_fmincon=optimoptions('fmincon','Display','iter','TolFun',1e-7,'MaxFunctionEvaluations',50000,'UseParallel',true,'FiniteDifferenceStepSize',1e-3);
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

%calibrate the pis's in T-function
go_calibrate_pis;

%initial guess for smart containment
cr_guess=crss*ones(HH,1)+crss*ones(HH,1).*randn(HH,1)*0.1-200;
%ci_guess=crss*ones(HH,1);

guess=[cr_guess];

%minimize negative PV utility (i.e. max PV utility) to find opt. paths for n's and c's.
%model equations are solved inside getU.m
LB=guess*0+10; %lower bound for consumption
UB=guess*3;%upper bound for consumption
[smart_sol,fval,exitflag] = fmincon(@getU,guess,[],[],[],[],LB,UB,@nonlinconstraint,opts_fmincon,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,HH,crss,nrss,Urss,phii,deltav,deltac,kappa,Uiss);

%save solution for possible use in subsequent maximization
save last_solution_smart_containment smart_sol;

%check solution
if exitflag~=1
    disp('Fmincon has not (yet) computed maximum PV utility, U0.');
end

%Get all endogenous variabels given optimal paths for n's and c's
[U0,c,ceq,I,S,R,D,T,Pop,cs,ns,Us,RnotSIRmacro,aggC,aggH,ci,cr,ni,nr,Ui,Ur,U,pid_endo]=getU(smart_sol,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,HH,crss,nrss,Urss,phii,deltav,deltac,kappa,Uiss);
RnotSIRmacro
U0

%plotting
ia=2;ib=3;fsize=12;
horz=HH;
time=0:1:horz-1;

figure;
subplot(ia,ib,1)
plot(time,100*I(1:horz),'b-','LineWidth',2);
box off;
title('Infected, I','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,2)
plot(time,100*S(1:horz),'b-','LineWidth',2);
box off;
title('Susceptibles, S','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,3)
plot(time,100*D(1:horz),'b-','LineWidth',2);
box off;
title('Deaths, D','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,5)
plot(time,100*(cs(1:horz)-crss)/crss,'b-','LineWidth',2); hold on
plot(time,100*(ci(1:horz)-crss)/crss,'r--','LineWidth',2); hold on
plot(time,100*(cr(1:horz)-crss)/crss,'k-.','LineWidth',2); hold off
box off;
title('Consumption by Type','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
leg=legend('Susceptibles','Infected','Recovered','Location','best');
set(leg,'FontSize',fsize-3);
legend boxoff

subplot(ia,ib,6)
plot(time,100*(ns(1:horz)-nrss)/nrss,'b-','LineWidth',2); hold on
plot(time,100*(ni(1:horz)-nrss)/nrss,'r--','LineWidth',2); hold on
plot(time,100*(nr(1:horz)-nrss)/nrss,'k-.','LineWidth',2); hold off
box off;
title('Hours by Type','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
leg=legend('Susceptibles','Infected','Recovered','Location','best');
set(leg,'FontSize',fsize-3);
legend boxoff

subplot(ia,ib,4)
plot(time,0*100*(aggC(1:horz)-crss)/crss,'k:','LineWidth',1.5);hold on
plot(time,100*(aggC(1:horz)-crss)/crss,'b-','LineWidth',2);
box off;
title('Agg. Consumption, C','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
set(gca,'FontSize',fsize);
xlabel('Weeks','FontSize',fsize);


suptitle('Figure 10: Smart Containment in the Benchmark SIR-Macro Model');
orient landscape
print -dpdf -fillpage epidemic_simulation_fig10


ia=2;ib=2;fsize=12;
horz=HH;
time=0:1:horz-1;

figure;
subplot(ia,ib,1)
plot(time,100*I(1:horz),'b-','LineWidth',2);
box off;
title('Infected, I','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,2)
plot(time,100*S(1:horz),'b-','LineWidth',2);
box off;
title('Susceptibles, S','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,3)
plot(time,100*D(1:horz),'b-','LineWidth',2);
box off;
title('Deaths, D','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);
xlabel('Weeks','FontSize',fsize);

subplot(ia,ib,4)
plot(time,0*100*(aggC(1:horz)-crss)/crss,'k:','LineWidth',1.5);hold on
plot(time,100*(aggC(1:horz)-crss)/crss,'b-','LineWidth',2);
box off;
title('Agg. Consumption, C','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
set(gca,'FontSize',fsize);
xlabel('Weeks','FontSize',fsize);


suptitle('Figure 10: Smart Containment in the Benchmark SIR-Macro Model');
orient landscape
print -dpdf -fillpage epidemic_simulation_fig10

%Save workspace
save results_smart_containment



%output some data used in paper and robustness table
aggCons_trough_percent=min((100*(aggC-crss)/crss))
aggCons_avg_first_year_percent=mean((100*(aggC(1:52)-crss)/crss))
terminal_one_minus_susceptibles_percent=100*(1-S(end))
peak_infection_percent=max(100*I)
terminal_death_share_percent=100*D(end)
terminal_number_deaths_US_millions=terminal_death_share_percent/100*330

toc;

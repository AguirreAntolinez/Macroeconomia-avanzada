function [f,c,ceq,I,S,R,D,T,Pop,cs,ns,Us,RnotSIRmacro,aggC,aggH,ci,cr,ni,nr,Ui,Ur,U,pid_endo]=getU(smart_guess,A,theta,i_ini,pop_ini,pis1,pis2,pis3,pir,pid,betta,HH,crss,nrss,Urss,phii,deltav,deltac,kappa,Uiss)

%back out guesses for cr 
cr=smart_guess(1:HH);

%calculate allocations based on planner FOCs:
cr(1)=crss; %normalize initial consumption of recovered to steady state (since planner does not affect that since R(0)=0)
nr=A/theta./cr; %Planner optimality condition for cr
ni=nr*0; %hours infected
ci=0*cr+crss*0.0001; %set consumption of infected close to zero
ns=nr; %hours susceptibles
cs=cr; %consumption susceptibles

%equilibrium equations of planner setup

%pre-allocate
I=NaN*ones(HH+1,1);
S=NaN*ones(HH+1,1);
D=NaN*ones(HH+1,1);
R=NaN*ones(HH+1,1);
Pop=NaN*ones(HH+1,1);
T=NaN*ones(HH,1);

%initial conditions
Pop(1)=pop_ini;
I(1)=i_ini;
S(1)=Pop(1)-I(1);
D(1)=0;
R(1)=0;

%%Endogenous death probability
pid_endo=NaN*ones(HH,1);

%iterate on SIR equations
for j=1:1:HH    
    %epi and population dynamics
    T(j,1)=pis1*S(j)*cs(j)*I(j)*ci(j)+pis2*S(j)*ns(j)*I(j)*ni(j)+pis3*S(j)*I(j);
    pid_endo(j,1)=pid+kappa*I(j)^2;
    S(j+1,1)=S(j)-T(j);
    I(j+1,1)=I(j)+T(j)-(pir+pid_endo(j,1))*I(j);
    R(j+1,1)=R(j)+pir*I(j);
    D(j+1,1)=D(j)+pid_endo(j,1)*I(j);
    Pop(j+1,1)=Pop(j,1)-pid_endo(j,1)*I(j);
end

%Aggregate consumption and hours
aggC=S(1:HH).*cs+I(1:HH).*ci+R(1:HH).*cr;
aggH=S(1:HH).*ns+I(1:HH).*ni*phii+R(1:HH).*nr;

%resource constraint error
err=aggC-A*aggH;


%utility recovered
ur=log(cr)-theta/2*nr.^2;
Ur=NaN*ones(HH+1,1);
Ur(HH+1)=Urss;
for tt=HH:-1:1
    Ur(tt,1)=ur(tt)+betta*Ur(tt+1,1);
end

%utility infected
ui=log(ci)-theta/2*ni.^2;
Ui=NaN*ones(HH+1,1);
Ui(HH+1)=(1-deltac)^HH*Uiss+(1-(1-deltac)^HH)*Urss;%terminal condition
for tt=HH:-1:1
    Ui(tt,1)=ui(tt)+(1-deltac)*betta*((1-pir-pid_endo(tt))*Ui(tt+1,1)+pir*Ur(tt+1,1))+deltac*betta*Ur(tt+1,1);
end

%utility susceptibles
us=log(cs)-theta/2*ns.^2;
Us=NaN*ones(HH+1,1);
Usss=Urss; %PV utility of susceptibles same as recovered in steady state
Us(HH+1)=(1-deltav)^HH*Usss+(1-(1-deltav)^HH)*Urss;%terminal condition
for tt=HH:-1:1
    Us(tt,1)=us(tt)+(1-deltav)*betta*(1-T(tt)).*Us(tt+1,1)+(1-deltav)*betta*T(tt).*Ui(tt+1,1)+deltav*betta*Ur(tt+1,1);
end

%Present value of society utility
U=S(1:HH).*Us(1:HH)+I(1:HH).*Ui(1:HH)+R(1:HH).*Ur(1:HH);

RnotSIRmacro=T(1)/I(1)/(pir+pid);

%PV Utility at time, t=0
%Use minimizer (fmincon) to find optimal policy muc such that
%U1 is maximal -- hence flip sign
U0=-U(1);

%output of function needed for fmincon
f=U0; %objective
c=[];
ceq=err; %nonlinear equality constraints




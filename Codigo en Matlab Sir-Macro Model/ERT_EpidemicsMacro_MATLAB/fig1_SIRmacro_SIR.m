%Baseline
cd baseline;
load results_baseline;
cd ..

Ibase=I;
Sbase=S;
Rbase=R;
Dbase=D;
aggCbase=aggC;
aggHbase=aggH;
crbase=crss;
nrbase=nrss; 


%Perturbation 1
cd canonicalSIR;
load results_canonicalSIR;
cd ..

Iperturb1=I;
Sperturb1=S;
Rperturb1=R;
Dperturb1=D;
aggCperturb1=aggC;
aggHperturb1=aggH;
crperturb1=crss;
nrperturb1=nrss; 


%Perturbation 2
cd baseline;
load results_baseline;
cd ..
Iperturb2=I*NaN;
Sperturb2=S*NaN;
Rperturb2=R*NaN;
Dperturb2=D*NaN;
aggCperturb2=aggC*NaN;
aggHperturb2=aggH*NaN;
crperturb2=crss;
nrperturb2=nrss; 


%plotting
ia=2;
ib=3;
fsize=12;
horz=HH;

horz=150;

time=0:1:horz-1;

figure;
subplot(ia,ib,2)
plot(time,100*Sbase(1:horz),'b-','LineWidth',2);hold on
plot(time,100*Sperturb1(1:horz),'m--','LineWidth',2);hold on
plot(time,100*Sperturb2(1:horz),'r-.','LineWidth',2);hold off
box off;
title('Susceptibles, S','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,3)
plot(time,100*Rbase(1:horz),'b-','LineWidth',2);hold on
plot(time,100*Rperturb1(1:horz),'m--','LineWidth',2);hold on
plot(time,100*Rperturb2(1:horz),'r-.','LineWidth',2);hold off
box off;
title('Recovered, R','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,4)
plot(time,100*Dbase(1:horz),'b-','LineWidth',2);hold on
plot(time,100*Dperturb1(1:horz),'m--','LineWidth',2);hold on
plot(time,100*Dperturb2(1:horz),'r-.','LineWidth',2);hold off
box off;
title('Deaths, D','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,5)
plot(time,0*100*(aggCbase(1:horz)-crbase)/crbase,'m:','LineWidth',1.5);hold on
plot(time,100*(aggCbase(1:horz)-crbase)/crbase,'b-','LineWidth',2);hold on
plot(time,100*(aggCperturb1(1:horz)-crperturb1)/crperturb1,'m--','LineWidth',2);hold on
plot(time,100*(aggCperturb2(1:horz)-crperturb2)/crperturb2,'r-.','LineWidth',2);hold off
box off;
title('Aggregate Consumption, C','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,6)
plot(time,0*100*(aggCbase(1:horz)-crbase)/crbase,'m:','LineWidth',1.5);hold on
plot(time,100*(aggHbase(1:horz)-nrbase)/nrbase,'b-','LineWidth',2);hold on
plot(time,100*(aggHperturb1(1:horz)-nrperturb1)/nrperturb1,'m--','LineWidth',2);hold on
plot(time,100*(aggHperturb2(1:horz)-nrperturb2)/nrperturb2,'r-.','LineWidth',2);hold off
box off;
title('Aggregate Hours, N','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,1)
plot(time,100*Ibase(1:horz),'b-','LineWidth',2);hold on
plot(time,100*Iperturb1(1:horz),'m--','LineWidth',2);hold on
plot(time,100*Iperturb2(1:horz),'r-.','LineWidth',2);hold off
box off;
title('Infected, I','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);


leg1=legend('Basic SIR-Macro Model','SIR Model (\pi_{1}=\pi_{2}=0, model recalibrated)');
set(leg1,...
    'Position',[-0.0704335567178815 0.894078940309976 1.19107142857143 0.044047619047619],...
    'Orientation','horizontal');
legend boxoff;

suptitle2('Figure 1: Basic SIR-Macro Model vs. SIR Model');
orient landscape
print -dpdf -fillpage epidemic_simulation_fig1


 

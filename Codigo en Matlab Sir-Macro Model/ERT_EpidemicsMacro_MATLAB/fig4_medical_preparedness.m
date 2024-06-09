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
mucbase=muc;
pidbase=pid*18/7;


%Perturbation 1
cd medprep;
load results_medprep;
cd ..

Iperturb1=I;
Sperturb1=S;
Rperturb1=R;
Dperturb1=D;
aggCperturb1=aggC;
aggHperturb1=aggH;
crperturb1=crss;
nrperturb1=nrss; 
mucperturb1=muc;
pid_endo_perturb1=pid_endo*18/7; %convert to daily


%Perturbation 2
cd medprep_optimal;
load results_medprep_optimal;
cd ..
Iperturb2=I;
Sperturb2=S;
Rperturb2=R;
Dperturb2=D;
aggCperturb2=aggC;
aggHperturb2=aggH;
crperturb2=crss;
nrperturb2=nrss; 
mucperturb2=muc;
pid_endo_perturb2=pid_endo*18/7;



%plotting
ia=3;
ib=3;
fsize=12;
horz=HH;

horz=150;

time=0:1:horz-1;

figure;
subplot(ia,ib,2)
plot(time,100*Sbase(1:horz),'b-','LineWidth',2);hold on
plot(time,100*Sperturb1(1:horz),'r-.','LineWidth',2);hold on
plot(time,100*Sperturb2(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Sbase(1:horz),'b-','LineWidth',2);hold on
box off;
title('Susceptibles, S','FontSize',fsize);
ylabel('% of Ini. Pop.','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,3)
plot(time,100*Rbase(1:horz),'b-','LineWidth',2);hold on
plot(time,100*Rperturb1(1:horz),'r-.','LineWidth',2);hold on
plot(time,100*Rperturb2(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Rbase(1:horz),'b-','LineWidth',2);hold on
box off;
title('Recovered, R','FontSize',fsize);
ylabel('% of Ini. Pop.','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,4)
plot(time,100*Dbase(1:horz),'b-','LineWidth',2);hold on
plot(time,100*Dperturb1(1:horz),'r-.','LineWidth',2);hold on
plot(time,100*Dperturb2(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Dbase(1:horz),'b-','LineWidth',2);hold on
box off;
title('Deaths, D','FontSize',fsize);
ylabel('% of Ini. Pop.','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,5)
plot(time,0*100*(aggCbase(1:horz)-crbase)/crbase,'m:','LineWidth',1.5);hold on
plot(time,100*(aggCbase(1:horz)-crbase)/crbase,'b-','LineWidth',2);hold on
plot(time,100*(aggCperturb1(1:horz)-crperturb1)/crperturb1,'r-.','LineWidth',2);hold on
plot(time,100*(aggCperturb2(1:horz)-crperturb2)/crperturb2,'k--','LineWidth',2);hold on
plot(time,100*(aggCbase(1:horz)-crbase)/crbase,'b-','LineWidth',2);hold on
axis([0 horz -39 5]);
box off;
title('Aggregate Consumption, C','FontSize',fsize);
ylabel('% of Ini. St.St.','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,6)
plot(time,0*100*(aggCbase(1:horz)-crbase)/crbase,'m:','LineWidth',1.5);hold on
plot(time,100*(aggHbase(1:horz)-nrbase)/nrbase,'b-','LineWidth',2);hold on
plot(time,100*(aggHperturb1(1:horz)-nrperturb1)/nrperturb1,'r-.','LineWidth',2);hold on
plot(time,100*(aggHperturb2(1:horz)-nrperturb2)/nrperturb2,'k--','LineWidth',2);hold on
plot(time,100*(aggHbase(1:horz)-nrbase)/nrbase,'b-','LineWidth',2);hold on
axis([0 horz -39 5]);
box off;
title('Aggregate Hours, N','FontSize',fsize);
ylabel('% of Ini. St.St.','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);


subplot(ia,ib,7)
plot(time,100*pidbase*ones(horz,1),'b-','LineWidth',2);hold on
plot(time,100*pid_endo_perturb1(1:horz),'r-.','LineWidth',2);hold on
plot(time,100*pid_endo_perturb2(1:horz),'k--','LineWidth',2);hold off
axis([0 horz 0.4 1.1]);
box off;
title('Case Fatality Rate, \pi_d','FontSize',fsize);
ylabel('%','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);



subplot(ia,ib,8)
plot(time,100*mucbase(1:horz),'b-','LineWidth',2);hold on
plot(time,100*mucperturb1(1:horz),'r-.','LineWidth',2);hold on
plot(time,100*mucperturb2(1:horz),'k--','LineWidth',2);hold on
plot(time,100*mucbase(1:horz),'b-','LineWidth',2);hold on
axis([0 horz 0 115]);
box off;
title('Best Containment Policy, \mu','FontSize',fsize);
ylabel('%','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);


subplot(ia,ib,1)
plot(time,100*Ibase(1:horz),'b-','LineWidth',2);hold on
plot(time,100*Iperturb1(1:horz),'r-.','LineWidth',2);hold on
plot(time,100*Iperturb2(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Ibase(1:horz),'b-','LineWidth',2);hold on
box off;
title('Infected, I','FontSize',fsize);
ylabel('% of Ini. Pop.','FontSize',fsize);
set(gca,'FontSize',fsize);

leg1=legend('Basic SIR-Macro Model (\pi_d constant)','Endog. Case Fatality Rate (\pi_{d} = f ( Infected )  )','Best Simple Containment Policy');
set(leg1,...
    'Position',[-0.0704335567178815 0.894078940309976 1.19107142857143 0.044047619047619],...
    'Orientation','horizontal');
legend boxoff;

suptitle2('Figure 4: Medical Preparedness');
orient landscape
print -dpdf -fillpage epidemic_simulation_fig4


 

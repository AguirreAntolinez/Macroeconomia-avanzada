%Baseline
cd vaccine_treatment_medprep_optimal;
load results_vaccine_treatment_medprep_optimal;
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


%Perturbation 1
cd vaccine_treatment_medprep_optimal_earlyexit_after12w
load results_vaccine_treatment_medprep_optimal_earlyexit;
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


%Perturbation 2
cd vaccine_treatment_medprep_optimal_earlyexit_after44w
load results_vaccine_treatment_medprep_optimal_earlyexit;
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


%plotting
ia=2;
ib=4;
fsize=11;
horz=HH;

horz=120;

time=0:1:horz-1;

figure;

%Early Exit after 12 Weeks
subplot(ia,ib,2)
plot(time,100*Dbase(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Dperturb1(1:horz),'r-','LineWidth',2);hold on
plot(time,100*Dbase(1:horz),'k--','LineWidth',2);hold on
box off;
title('Deaths, D','FontSize',fsize);
%ylabel('% of Initial Population','FontSize',fsize-2);
set(gca,'FontSize',fsize);
axis([0 horz 0 0.4]);

subplot(ia,ib,3)
plot(time,0*100*(aggCbase(1:horz)-crbase)/crbase,'k:','LineWidth',1);hold on
plot(time,100*(aggCbase(1:horz)-crbase)/crbase,'k--','LineWidth',2);hold on
plot(time,100*(aggCperturb1(1:horz)-crperturb1)/crperturb1,'r-','LineWidth',2);hold on
plot(time,100*(aggCbase(1:horz)-crbase)/crbase,'k--','LineWidth',2);hold on
box off;
title('Agg. Consumption, C','FontSize',fsize);
%ylabel('% Dev. from Initial Steady State','FontSize',fsize-1);
set(gca,'FontSize',fsize);
axis([0 horz -30 2]);

subplot(ia,ib,4)
plot(time,0*100*mucbase(1:horz),'k:','LineWidth',1);hold on
plot(time,100*mucbase(1:horz),'k--','LineWidth',2);hold on
plot(time,100*mucperturb1(1:horz),'r-','LineWidth',2);hold on
plot(time,0*100*mucbase(1:horz),'k:','LineWidth',1);hold on
box off;
title('Best Containment Policy, \mu','FontSize',fsize);
%ylabel('%','FontSize',fsize-1);
set(gca,'FontSize',fsize);
axis([0 horz -10 80]);


subplot(ia,ib,1)
plot(time,100*Ibase(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Iperturb1(1:horz),'r-','LineWidth',2);hold on
plot(time,100*Ibase(1:horz),'k--','LineWidth',2);hold on
box off;
title('Infected, I','FontSize',fsize);
%ylabel('% of Initial Population','FontSize',fsize-1);
set(gca,'FontSize',fsize);
axis([0 horz 0 5]);

xlim=get(gca,'XLim');
ylim=get(gca,'YLim');
ht = text(xlim(1)-70,ylim(1),'Panel A: Exit after 12 Weeks');
set(ht,'Rotation',90)
set(ht,'FontSize',fsize+2,'FontWeight','bold');



 







%Early Exit after 44 Weeks

subplot(ia,ib,6)
plot(time,100*Dbase(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Dperturb2(1:horz),'r-','LineWidth',2);hold on
plot(time,100*Dbase(1:horz),'k--','LineWidth',2);hold on
box off;
title('Deaths, D','FontSize',fsize);
%ylabel('% of Initial Population','FontSize',fsize-1);
%xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
axis([0 horz 0 0.4]);

subplot(ia,ib,7)
plot(time,0*100*(aggCbase(1:horz)-crbase)/crbase,'k:','LineWidth',1);hold on
plot(time,100*(aggCbase(1:horz)-crbase)/crbase,'k--','LineWidth',2);hold on
plot(time,100*(aggCperturb2(1:horz)-crperturb2)/crperturb2,'r-','LineWidth',2);hold on
plot(time,100*(aggCbase(1:horz)-crbase)/crbase,'k--','LineWidth',2);hold on
box off;
title('Agg. Consumption, C','FontSize',fsize);
%ylabel('% Dev. from Initial Steady State','FontSize',fsize-1);
%xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
axis([0 horz -30 2]);

subplot(ia,ib,8)
plot(time,0*100*mucbase(1:horz),'k:','LineWidth',1);hold on
plot(time,100*mucbase(1:horz),'k--','LineWidth',2);hold on
plot(time,100*mucperturb2(1:horz),'r-','LineWidth',2);hold on
plot(time,100*mucbase(1:horz),'k--','LineWidth',2);hold on
box off;
title('Best Containment Policy, \mu','FontSize',fsize);
%ylabel('%','FontSize',fsize-1);
%xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
axis([0 horz -10 80]);


subplot(ia,ib,5)
plot(time,100*Ibase(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Iperturb2(1:horz),'r-','LineWidth',2);hold on
plot(time,100*Ibase(1:horz),'k--','LineWidth',2);hold on
box off;
title('Infected, I','FontSize',fsize);
%ylabel('% of Initial Population','FontSize',fsize-1);
%xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
axis([0 horz 0 5]);

xlim=get(gca,'XLim');
ylim=get(gca,'YLim');
ht = text(xlim(1)-70,ylim(1),'Panel B: Exit after 44 Weeks');
set(ht,'Rotation',90)
set(ht,'FontSize',fsize+2,'FontWeight','bold');

xlim=get(gca,'XLim');
ylim=get(gca,'YLim');
ht = text(xlim(1)-50,ylim(1)-1.5,'Notes: x-axis in weeks; infected and deaths in % of ini. population; consumption in % dev. from ini. steady state; best containment policy in %.');
set(ht,'Rotation',0)
set(ht,'FontSize',fsize,'FontWeight','normal');



leg1=legend('Best Simple Containment Policy in Benchmark Model','Early Exit from Best Simple Containment Policy');
set(leg1,...
    'Position',[-0.0704335567178815 0.894078940309976 1.19107142857143 0.044047619047619],...
    'Orientation','horizontal','FontSize',fsize+1);
legend boxoff;


suptitle4('Figure 8: Benchmark SIR-Macro Model (Vaccines, Treatment, Med. Preparedness)');
orient landscape
print -dpdf -fillpage epidemic_simulation_fig8


 

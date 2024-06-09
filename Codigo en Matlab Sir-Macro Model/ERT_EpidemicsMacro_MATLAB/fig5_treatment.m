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


%Perturbation 1
cd treatment;
load results_treatment;
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
cd treatment_optimal;
load results_treatment_optimal;
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
ib=3;
fsize=12;
horz=HH;

horz=150;

time=0:1:horz-1;

figure;
subplot(ia,ib,2)
plot(time,100*Sbase(1:horz),'b-','LineWidth',2);hold on
plot(time,NaN*100*Sperturb1(1:horz),'r-.o','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*Sperturb2(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Sperturb1(1:horz),'r-.','LineWidth',2,'MarkerSize',4);hold on
plot(time([1, 20, 30, 45, 55, 90, 140]),100*Sperturb1([1, 20, 30, 45, 55, 90, 140]),'ro','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*Sbase(1:horz),'b-','LineWidth',2);hold on
box off;
title('Susceptibles, S','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,3)
plot(time,100*Rbase(1:horz),'b-','LineWidth',2);hold on
plot(time,NaN*100*Rperturb1(1:horz),'r-.o','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*Rperturb2(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Rperturb1(1:horz),'r-.','LineWidth',2,'MarkerSize',4);hold on
plot(time([1, 20, 30, 45, 55, 90, 140]),100*Rperturb1([1, 20, 30, 45, 55, 90, 140]),'ro','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*Rbase(1:horz),'b-','LineWidth',2);hold on
box off;
title('Recovered, R','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,4)
plot(time,100*Dbase(1:horz),'b-','LineWidth',2);hold on
plot(time,NaN*100*Dperturb1(1:horz),'r-.o','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*Dperturb2(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Dperturb1(1:horz),'r-.','LineWidth',2,'MarkerSize',4);hold on
plot(time([1, 20, 30, 45, 55, 90, 140]),100*Dperturb1([1, 20, 30, 45, 55, 90, 140]),'ro','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*Dbase(1:horz),'b-','LineWidth',2);hold on
box off;
title('Deaths, D','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);

subplot(ia,ib,5)
plot(time,0*100*(aggCbase(1:horz)-crbase)/crbase,'m:','LineWidth',1.5);hold on
plot(time,100*(aggCbase(1:horz)-crbase)/crbase,'b-','LineWidth',2);hold on
plot(time,NaN*100*(aggCperturb1(1:horz)-crperturb1)/crperturb1,'r-.o','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*(aggCperturb2(1:horz)-crperturb2)/crperturb2,'k--','LineWidth',2);hold on
plot(time,100*(aggCperturb1(1:horz)-crperturb1)/crperturb1,'r-.','LineWidth',2,'MarkerSize',4);hold on
plot(time([1, 20, 30, 45, 55, 90, 140]),100*(aggCperturb1([1, 20, 30, 45, 55, 90, 140])-crperturb1)/crperturb1,'ro','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*(aggCbase(1:horz)-crbase)/crbase,'b-','LineWidth',2);hold on
box off;
title('Aggregate Consumption, C','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);


subplot(ia,ib,6)
plot(time,100*mucbase(1:horz),'b-','LineWidth',2);hold on
plot(time,NaN*100*mucperturb1(1:horz),'r-.o','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*mucperturb2(1:horz),'k--','LineWidth',2);hold on
plot(time,100*mucperturb1(1:horz),'r-.','LineWidth',2,'MarkerSize',4);hold on
plot(time([1, 20, 30, 45, 55, 90, 140]),100*mucperturb1([1, 20, 30, 45, 55, 90, 140]),'ro','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*mucbase(1:horz),'b-','LineWidth',2);hold on
box off;
title('Best Containment Policy, \mu','FontSize',fsize);
ylabel('%','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);


subplot(ia,ib,1)
plot(time,100*Ibase(1:horz),'b-','LineWidth',2);hold on
plot(time,NaN*100*Iperturb1(1:horz),'r-.o','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*Iperturb2(1:horz),'k--','LineWidth',2);hold on
plot(time,100*Iperturb1(1:horz),'r-.','LineWidth',2,'MarkerSize',4);hold on
plot(time([1, 20, 30, 45, 55, 90, 140]),100*Iperturb1([1, 20, 30, 45, 55, 90, 140]),'ro','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','r');hold on
plot(time,100*Ibase(1:horz),'b-','LineWidth',2);hold on
box off;
title('Infected, I','FontSize',fsize);
ylabel('% of Initial Population','FontSize',fsize);
set(gca,'FontSize',fsize);


leg1=legend('Basic SIR-Macro Model','Model with Treatment','Best Simple Containment Policy with Treatment');
set(leg1,...
    'Position',[-0.0704335567178815 0.894078940309976 1.19107142857143 0.044047619047619],...
    'Orientation','horizontal');
legend boxoff;

suptitle2('Figure 5: SIR-Macro Model With Treatments');
orient landscape
print -dpdf -fillpage epidemic_simulation_fig5


 

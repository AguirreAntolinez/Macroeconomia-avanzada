cd vaccine_treatment_medprep_smart_containment

load results_smart_containment

cd .. 


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
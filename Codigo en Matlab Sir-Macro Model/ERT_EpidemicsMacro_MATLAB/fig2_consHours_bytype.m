%Plot Baseline

cd baseline;
load results_baseline;
cd ..

%plotting
ia=1;
ib=2;
fsize=12;
horz=HH;

horz=150;

time=0:1:horz-1;

figure;
subplot(ia,ib,1)
plot(time,100*(cs(1:horz)-crss)/crss,'b-','LineWidth',2); hold on
plot(time,100*(ci(1:horz)-crss)/crss,'r--','LineWidth',2); hold on
plot(time,100*(cr(1:horz)-crss)/crss,'k-.','LineWidth',2); hold off
box off;
axis([0 horz-1 -22 1.99]);
title('Consumption by Type','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
set(gca,'FontSize',fsize);
legend('Consumption Susceptibles','Consumption Infected','Consumption Recovered','Location','best');
legend boxoff

subplot(ia,ib,2)
plot(time,100*(ns(1:horz)-nrss)/nrss,'b-','LineWidth',2); hold on
plot(time,100*(ni(1:horz)-nrss)/nrss,'r--','LineWidth',2); hold on
plot(time,100*(nr(1:horz)-nrss)/nrss,'k-.','LineWidth',2); hold on
plot(time,100*(ni(1:horz)-nrss)/nrss,'r--','LineWidth',2); hold off
axis([0 horz-1 -22 1.99]);
box off;
title('Hours by Type','FontSize',fsize);
xlabel('Weeks','FontSize',fsize);
ylabel('% Dev. from Initial Steady State','FontSize',fsize);
set(gca,'FontSize',fsize); 
legend('Hours Susceptibles','Hours Infected','Hours Recovered','Location','best');
legend boxoff


suptitle('Figure 2: Consumption and Hours by Type in Basic SIR-Macro Model');
orient landscape
print -dpdf -fillpage epidemic_simulation_fig2

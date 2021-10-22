close all
flags.cost_control = 0; 

mod1 = load('blah1.mat');
mod2 = load('blah2.mat');
mod3 = load('blah3.mat');

figure
hold on
for k=1:15
    plot(mod1.t{k},sum(mod1.y{k}(:,[10;14;18;38;42;46;64;68;72]),2),'b-','linewidth',2);
    plot(mod2.t{k},sum(mod2.y{k}(:,[10;14;18;38;42;46;64;68;72]),2),'r--','linewidth',2);
    plot(mod3.t{k},sum(mod3.y{k}(:,[10;14;18;38;42;46;64;68;72]),2),'k-.','linewidth',2);
    
end
hold off
legend({'Normal','Placebo','Immune'},'location','best');

figure
hold on
for k=1:15
    plot(mod1.t{k},mod1.y{k}(:,27),'b-','linewidth',2);
    plot(mod2.t{k},mod2.y{k}(:,27),'r--','linewidth',2);
    plot(mod3.t{k},mod3.y{k}(:,27),'k-.','linewidth',2);
end
hold off
legend({'Normal','Placebo','Immune'},'location','best');

figure
hold on
for k=1:15
    plot(mod1.t{k},params.N_crit/params.N*(sum(mod1.y{k}(:,7:26),2) + sum(mod1.y{k}(:,35:54),2) + sum(mod1.y{k}(:,61:80),2)),'b-','linewidth',2);
    plot(mod2.t{k},params.N_crit/params.N*(sum(mod2.y{k}(:,7:26),2) + sum(mod2.y{k}(:,35:54),2) + sum(mod2.y{k}(:,61:80),2)),'r--','linewidth',2);
    plot(mod3.t{k},params.N_crit/params.N*(sum(mod3.y{k}(:,7:26),2) + sum(mod3.y{k}(:,35:54),2) + sum(mod3.y{k}(:,61:80),2)),'k-.','linewidth',2);
end
plot(316*ones(1000,1),linspace(0,1,1000),'k--');
hold off
legend({'Normal','Placebo','Immune'},'location','best');

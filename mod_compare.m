close all
flags.cost_control = 0; 

win_len = 15;

mod1 = load('blah1.mat');
mod2 = load('blah2.mat');
mod3 = load('blah3.mat');

figure
hold on
for k=1:win_len
    plot(mod1.t{k},params.N_crit*sum(mod1.y{k}(:,[10;14;18;38;42;46;64;68;72]),2),'b-','linewidth',2);
    plot(mod2.t{k},params.N_crit*sum(mod2.y{k}(:,[10;14;18;38;42;46;64;68;72]),2),'r--','linewidth',2);
    plot(mod3.t{k},params.N_crit*sum(mod3.y{k}(:,[10;14;18;38;42;46;64;68;72]),2),'k-.','linewidth',2);
    
end
plot(DATA_T,DATA_pos,'go');
hold off
legend({'Normal','Placebo','Immune'},'location','best');

figure
hold on
for k=1:win_len
    plot(mod1.t{k},params.N_crit*mod1.y{k}(:,27),'b-','linewidth',2);
    plot(mod2.t{k},params.N_crit*mod2.y{k}(:,27),'r--','linewidth',2);
    plot(mod3.t{k},params.N_crit*mod3.y{k}(:,27),'k-.','linewidth',2);
end
plot(DATA_T,DATA_tot,'go');
hold off
legend({'Normal','Placebo','Immune'},'location','best');

figure
hold on
for k=1:win_len
    plot(mod1.t{k},params.N_crit/params.N*(sum(mod1.y{k}(:,7:26),2) + sum(mod1.y{k}(:,35:54),2) + sum(mod1.y{k}(:,61:80),2)),'b-','linewidth',2);
    plot(mod2.t{k},params.N_crit/params.N*(sum(mod2.y{k}(:,7:26),2) + sum(mod2.y{k}(:,35:54),2) + sum(mod2.y{k}(:,61:80),2)),'r--','linewidth',2);
    plot(mod3.t{k},params.N_crit/params.N*(sum(mod3.y{k}(:,7:26),2) + sum(mod3.y{k}(:,35:54),2) + sum(mod3.y{k}(:,61:80),2)),'k-.','linewidth',2);
end
plot(316*ones(1000,1),linspace(0,1,1000),'k--');
hold off
legend({'Normal','Placebo','Immune'},'location','best');

figure
hold on
for k=1:win_len
    plot(mod1.t{k},(sum(mod1.y{k}(:,7:26),2) + sum(mod1.y{k}(:,35:54),2) + sum(mod1.y{k}(:,61:80),2))./mod1.y{k}(:,27),'b-','linewidth',2);
    plot(mod2.t{k},(sum(mod2.y{k}(:,7:26),2) + sum(mod2.y{k}(:,35:54),2) + sum(mod2.y{k}(:,61:80),2))./mod2.y{k}(:,27),'r--','linewidth',2);
    plot(mod3.t{k},(sum(mod3.y{k}(:,7:26),2) + sum(mod3.y{k}(:,35:54),2) + sum(mod3.y{k}(:,61:80),2))./mod3.y{k}(:,27),'k-.','linewidth',2);
end
plot(316*ones(1000,1),linspace(0,1,1000),'k--');
hold off
legend({'Normal','Placebo','Immune'},'location','best');

figure
hold on
for k=1:win_len
    plot(mod1.t{k},params.N_crit*(mod1.Params{k}.rho0*(sum(mod1.y{k}(:,1:10),2)+sum(mod1.y{k}(:,15:26),2))+mod1.Params{k}.rhoI*(sum(mod1.y{k}(:,11:14),2))+mod1.Params{k}.rhoV0*(sum(mod1.y{k}(:,29:38),2)+sum(mod1.y{k}(:,43:64),2)+sum(mod1.y{k}(:,69:80),2))+mod1.Params{k}.rhoVI*(sum(mod1.y{k}(:,39:42),2)+sum(mod1.y{k}(:,65:68),2))),'b-','linewidth',2);
    plot(mod2.t{k},params.N_crit*(mod2.Params{k}.rho0*(sum(mod2.y{k}(:,1:10),2)+sum(mod2.y{k}(:,15:26),2))+mod2.Params{k}.rhoI*(sum(mod2.y{k}(:,11:14),2))+mod2.Params{k}.rhoV0*(sum(mod2.y{k}(:,29:38),2)+sum(mod2.y{k}(:,43:64),2)+sum(mod2.y{k}(:,69:80),2))+mod2.Params{k}.rhoVI*(sum(mod2.y{k}(:,39:42),2)+sum(mod2.y{k}(:,65:68),2))),'r--','linewidth',2);
    plot(mod3.t{k},params.N_crit*(mod3.Params{k}.rho0*(sum(mod3.y{k}(:,1:10),2)+sum(mod3.y{k}(:,15:26),2))+mod3.Params{k}.rhoI*(sum(mod3.y{k}(:,11:14),2))+mod3.Params{k}.rhoV0*(sum(mod3.y{k}(:,29:38),2)+sum(mod3.y{k}(:,43:64),2)+sum(mod3.y{k}(:,69:80),2))+mod3.Params{k}.rhoVI*(sum(mod3.y{k}(:,39:42),2)+sum(mod3.y{k}(:,65:68),2))),'k-.','linewidth',2);
end
plot(DATA_T,DATA_test,'go');
hold off
legend({'Normal','Placebo','Immune'},'location','best');
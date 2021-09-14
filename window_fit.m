close all
clear variables

date_flag = 1;

if date_flag
    base_date = datenum('10-Mar-2020','dd-mmm-yy');
else
    base_date = 0;
end

run_params;
flags.cost_control=0;
flags.model = 3;
flags.lsq_refine = 1; %turn on if you want to refine the parameter fitting

window(1) = 150;
window(2) = 308;
window(3) = 384;%384;%350;
window(4) = 407;%407;%497;
window(5) = 497;

%Use the Ontario data to fit Kc, K0, Mc, M0, and initial sick
year = '2021';
mday = 'Jul20';

load(['DATA-',year,'-',mday,'.mat']);

DATA_pos = DATA_pos - DATA_pos(1);
DATA_tot = DATA_tot - DATA_tot(1);

spot = zeros(size(window));

for k=1:length(window)
    spot(k) = find(DATA_T==window(k));
end

%Initialize parameters
v0 = [0.0607 0.0263 0.0095];

LB = [0;0;0];

UB =[inf;inf;inf];

flags.diffy0 = 0;

var_names = {'S';'S_1';'S_2';'E';'E_1';'E_2';'P';'P_1';'P_2';'P_M';'I_S';'I_{S1}';'I_{S2}';'I_{SM}';'I_A';'I_{A1}';'I_{A2}';'I_{AM}';'R_S';'R_{S1}';'R_{S2}';'R_{SM}';'R_A';'R_{A1}';'R_{A2}';'R_{AM}';'M';'C'};

options = optimoptions('lsqcurvefit','UseParallel',true);
odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);

flags.cases = 0;
for k=1:length(window)
    fprintf('Fitting window %i\n',k);
    if k==1
        current_times = DATA_T(1:window(1));
        current_pos = DATA_pos(1:window(1));
        current_tot = DATA_tot(1:window(1));
    else
        current_times = DATA_T(window(k-1):window(k));
        current_pos = DATA_pos(window(k-1):window(k));
        current_tot = DATA_tot(window(k-1):window(k));
        v0 = V{k-1};
        y0 = y{k-1}(end,:);
    end
    if flags.lsq_refine
        %Update parameters to convergence
        converge_done = 0;
        converge_tol = 1E-8;
        while ~converge_done
            [V{k},RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@covid_lsq,v0,current_times,[current_pos/params.N_crit,current_tot/params.N_crit],LB,UB,options,flags,y0);
%             norm(V{k}-v0,2)
            if norm(V{k}-v0,2)<=converge_tol
                converge_done = 1;
            end
            v0 = V{k};

        end
    else
        [V{k},RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@covid_lsq,v0,current_times,[current_pos/params.N_crit,current_tot/params.N_crit],LB,UB,options,flags,y0);
    end
    
    switch flags.cases
        case 0
            params.Kc = V{k}(1);
            params.Mc = V{k}(2);
            params.rho0 = V{k}(3);

            params.rhoI = 4*V{k}(3);
            params.K0 = 4*params.Kc;
            params.M0 = 2*params.Mc;
    end

    params.beta = params.R0*params.phi*params.gamma/(params.alpha*(params.gamma+params.phi)+(1-params.alpha)*params.q*params.phi);
    
    [t{k},y{k}] = ode23t(@ODEf,current_times,y0,odeopts,params,flags);
    
    PM{k} = y{k}(:,10);
    IsM{k} = y{k}(:,14);
    IaM{k} = y{k}(:,18);
    M{k} = y{k}(:,27);

    AM{k} = PM{k} + IsM{k} + IaM{k};
    AT{k} = y{k}(:,7)+y{k}(:,8)+y{k}(:,9)+y{k}(:,10)+y{k}(:,11)+y{k}(:,12)+y{k}(:,13)+y{k}(:,14)+y{k}(:,15)+y{k}(:,16)+y{k}(:,17)+y{k}(:,18);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
for k=1:length(window)
    plot(t{k}+base_date,y{k}(:,1)+y{k}(:,2)+y{k}(:,3),'linewidth',2);
end
hold off

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end

figure
hold on
for k=1:length(window)
    plot(t{k}+base_date,y{k}(:,2)+y{k}(:,3),'linewidth',2);
end
hold off

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end

figure
hold on
for k=1:length(window)
    plot(t{k}+base_date,y{k}(:,28),'linewidth',2);
end
hold off

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end

figure
plot(t{1}+base_date,AM{1},'linewidth',2);
hold on
for k=2:length(window)
    plot(t{k}+base_date,AM{k},'linewidth',2);
end
plot(DATA_T+base_date,(DATA_pos)/params.N0,'o');
hold off

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end

figure
plot(t{1}+base_date,M{1},'linewidth',2);
hold on
for k=2:length(window)
    plot(t{k}+base_date,M{k},'linewidth',2);
end
plot(DATA_T+base_date,(DATA_tot)/params.N0,'o');
hold off

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end
close all
clear variables

date_flag = 0;

if date_flag
    base_date = datenum('10-Mar-2020','dd-mmm-yy');
else
    base_date = 0;
end

flags.cost_control=0;
flags.model = 4;
flags.lsq_refine = 0; %turn on if you want to refine the parameter fitting
flags.fit_tot = 0; %turn on if you want to fit to total data as well as active data
flags.fit_test = 0; %turn on if you want to fit to test data

run_params;

% window = [150;308;384;407;552];

% window=[89;163;290;309;342;393;457;476;492;587];
window=[89;163;290;315;324;342;367;393;425;457;462;476;492;540;587];
% window=[35;89;163;290;315;324;342;367;393;425;457;462;476;492;540;587];
%vaccine change at 309,393,476
param_flag = [0;0;0;0;2;2;2;2;2;2;2;2;2;2;2];
% param_flag = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];

%for param flag 0 means behaviour(b) 1 means  vaccine(v), 2 means b and v

% window=[93;163;291;310;343;394;458;477;493;552];


%Use the Ontario data to fit Kc, K0, Mc, M0, and initial sick
year = '2021';
mday = 'Oct18';

load(['DATA-',year,'-',mday,'.mat']);
Vdat = load(['vaccine_DATA-',year,'-',mday,'.mat']);

DATA_pos = DATA_pos - DATA_pos(1);
DATA_tot = DATA_tot - DATA_tot(1);

DATA_test = [DATA_test(2:end);0];
%Assume that first days of COVID had about an 8% testing rate as well
test_spot = find(DATA_test~=0,1);
for k=2:test_spot-1
    DATA_test(k)=(DATA_tot(k)-DATA_tot(k-1))/0.08;
end

DATA_vac = zeros(length(DATA_pos),1);
DATA_vac(Vdat.DATA_T(2+14):end) = Vdat.DATA_tot1(1:end-14);

spot = zeros(size(window));

for k=1:length(window)
    spot(k) = find(DATA_T==window(k));
end

%Initialize parameters
v0 = [0.0607 0.0263 0.0095 0];

LB = [0;0;0;0];

UB =[inf;inf;inf;inf];

flags.diffy0 = 0;

switch flags.model
    case 3
        var_names = {'S';'S_1';'S_2';'E';'E_1';'E_2';'P';'P_1';'P_2';'P_M';'I_S';'I_{S1}';'I_{S2}';'I_{SM}';'I_A';'I_{A1}';'I_{A2}';'I_{AM}';'R_S';'R_{S1}';'R_{S2}';'R_{SM}';'R_A';'R_{A1}';'R_{A2}';'R_{AM}';'M';'C'};
    case 4
        var_names = {'S';'S_1';'S_2';'E';'E_1';'E_2';'P';'P_1';'P_2';'P_M';'I_S';'I_{S1}';'I_{S2}';'I_{SM}';'I_A';'I_{A1}';'I_{A2}';'I_{AM}';'R_S';'R_{S1}';'R_{S2}';'R_{SM}';'R_A';'R_{A1}';'R_{A2}';'R_{AM}';'M';'C';...
            'VS';'VS_1';'VS_2';'VE';'VE_1';'VE_2';'VP';'VP_1';'VP_2';'VP_M';'VI_S';'VI_{S1}';'VI_{S2}';'VI_{SM}';'VI_A';'VI_{A1}';'VI_{A2}';'VI_{AM}';'VR_S';'VR_{S1}';'VR_{S2}';'VR_{SM}';'VR_A';'VR_{A1}';'VR_{A2}';'VR_{AM}';...
            'WS';'WS_1';'WS_2';'WE';'WE_1';'WE_2';'WP';'WP_1';'WP_2';'WP_M';'WI_S';'WI_{S1}';'WI_{S2}';'WI_{SM}';'WI_A';'WI_{A1}';'WI_{A2}';'WI_{AM}';'WR_S';'WR_{S1}';'WR_{S2}';'WR_{SM}';'WR_A';'WR_{A1}';'WR_{A2}';'WR_{AM}'};
end

options = optimoptions('lsqcurvefit','UseParallel',true);
odeopts = odeset('NonNegative',1,'RelTol',1e-8,'AbsTol',1e-9);

% flags.cases = 0;
for k=1:length(window)
    fprintf('Fitting window %i\n',k);
    flags.cases = param_flag(k);
    if k==1
        current_times = DATA_T(1:window(1));
        current_pos = DATA_pos(1:window(1));
        current_tot = DATA_tot(1:window(1));
        current_test = DATA_test(1:window(1));
        current_vac = DATA_vac(1:window(1));
    else
        current_times = DATA_T(window(k-1):window(k));
        current_pos = DATA_pos(window(k-1):window(k));
        current_tot = DATA_tot(window(k-1):window(k));
        current_test = DATA_test(window(k-1):window(k));
        current_vac = DATA_vac(window(k-1):window(k));
        v0 = V{k-1};
        y0 = y{k-1}(end,:);
    end
    %Update parameters to convergence
    converge_done = 0;
    converge_tol = 1E-8;
    case_dat = [current_pos/params.N_crit];
    if flags.fit_tot
        case_dat = [case_dat;current_tot/params.N_crit];
    end
    if flags.fit_test
        case_dat = [case_dat;current_test/params.N_crit];
    end
    while ~converge_done
        switch flags.cases
            case 0
                [V{k},RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@covid_lsq,v0(1:3),current_times,case_dat,LB(1:3),UB(1:3),options,flags,y0,v0);
                V{k}(4) = v0(4);
            case 1
                [temp,RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@covid_lsq,v0(4),current_times,[current_vac/params.N_crit],LB(4),UB(4),options,flags,y0,v0);
                V{k}(1:3) = v0(1:3);
                V{k}(4) = temp;
            case 2
                [V{k},RESNORM,RESIDUAL,EXITFLAG] = lsqcurvefit(@covid_lsq,v0,current_times,[case_dat;current_vac/params.N_crit],LB,UB,options,flags,y0);
        end
        if flags.lsq_refine
            if norm(V{k}-v0,2)<=converge_tol
                converge_done = 1;
            end
            v0 = V{k};
        else
            converge_done = 1;
        end
    end
    
    params.Kc = V{k}(1);
    params.Mc = V{k}(2);
    params.rho0 = V{k}(3);

    params.rhoI = 4*params.rho0;
    params.K0 = 4*params.Kc;
    params.M0 = 2*params.Mc;

    params.p = V{k}(4);
    params.rhoV0 = params.rho0/10;
    params.rhoVI = 4*params.rhoV0;

    params.beta = params.R0*params.phi*params.gamma/(params.alpha*(params.gamma+params.phi)+(1-params.alpha)*params.q*params.phi);
    
    Params{k} = params;
    
    [t{k},y{k}] = ode23t(@ODEf,current_times,y0,odeopts,params,flags);
    
    PM{k} = y{k}(:,10);
    IsM{k} = y{k}(:,14);
    IaM{k} = y{k}(:,18);
    VPM{k} = y{k}(:,38);
    VIsM{k} = y{k}(:,42);
    VIaM{k} = y{k}(:,46);
    WPM{k} = y{k}(:,64);
    WIsM{k} = y{k}(:,68);
    WIaM{k} = y{k}(:,72);
    M{k} = y{k}(:,27);

    AM{k} = PM{k} + IsM{k} + IaM{k} + VPM{k} + VIsM{k} + VIaM{k} + WPM{k} + WIsM{k} +WIaM{k};
    AT{k} = y{k}(:,7)+y{k}(:,8)+y{k}(:,9)+y{k}(:,10)+y{k}(:,11)+y{k}(:,12)+y{k}(:,13)+y{k}(:,14)+y{k}(:,15)+y{k}(:,16)+y{k}(:,17)+y{k}(:,18);
    VT{k} = sum(y{k}(:,29:80),2);
    sick{k} = sum(y{k}(:,7:26),2) + sum(y{k}(:,35:54),2) + sum(y{k}(:,61:80),2);
    tests{k} = params.rho0*(sum(y{k}(:,1:10),2)+sum(y{k}(:,15:26),2))+params.rhoI*(sum(y{k}(:,11:14),2))+params.rhoV0*(sum(y{k}(:,29:38),2)+sum(y{k}(:,43:64),2)+sum(y{k}(:,69:80),2))+params.rhoVI*(sum(y{k}(:,39:42),2)+sum(y{k}(:,65:68),2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
for k=1:length(window)
    plot(t{k}+base_date,params.N_crit/params.N*(y{k}(:,1)+y{k}(:,2)+y{k}(:,3)+y{k}(:,55)+y{k}(:,56)+y{k}(:,57)),'linewidth',2);
    plot(t{k}+base_date,params.N_crit/params.N*(y{k}(:,29)+y{k}(:,30)+y{k}(:,31)),'--','linewidth',2);
end
hold off
title('Susceptible');

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end

figure
hold on
for k=1:length(window)
    plot(t{k}+base_date,sick{k}./M{k},'linewidth',2);
end
hold off
title('Ascertainment');

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end
ylim([0,30]);

figure
hold on
for k=1:length(window)
    plot(t{k}(2:end)+base_date,(M{k}(2:end)-M{k}(1:end-1))./tests{k}(2:end),'linewidth',2);
end
hold off
title('Positivity Rate');

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end

figure
hold on
for k=1:length(window)
    plot(t{k}+base_date,y{k}(:,2)+y{k}(:,3)+y{k}(:,5)+y{k}(:,6)+y{k}(:,8)+y{k}(:,9)+y{k}(:,12)+y{k}(:,13)+y{k}(:,16)+y{k}(:,17),'linewidth',2);
end
hold off
title('Distancing');

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
title('C');

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end

figure
plot(t{1}+base_date,AM{1}*params.N0,'linewidth',2);
hold on
for k=2:length(window)
    plot(t{k}+base_date,AM{k}*params.N0,'linewidth',2);
end
plot(DATA_T+base_date,(DATA_pos),'o');
hold off

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end

figure
plot(t{1}+base_date,M{1}*params.N0,'linewidth',2);
hold on
for k=2:length(window)
    plot(t{k}+base_date,M{k}*params.N0,'linewidth',2);
end
plot(DATA_T+base_date,(DATA_tot),'o');
hold off

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end

figure
hold on
for k=1:length(window)
    plot(t{k}+base_date,tests{k}*params.N0,'linewidth',2);
end
plot(DATA_T+base_date,(DATA_test),'o');
hold off

xlim([base_date,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end

figure
plot(t{1}+base_date,VT{1}*params.N0,'linewidth',2);
hold on
for k=2:length(window)
    plot(t{k}+base_date,VT{k}*params.N0,'linewidth',2);
end
plot(DATA_T+base_date,(DATA_vac),'o');
hold off

xlim([base_date+302,base_date+DATA_T(end)]);

if date_flag
    datetick('x','dd-mmm-yy','keepticks','keeplimits');
    xtickangle(45);
end
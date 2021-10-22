close all
clear variables

date_flag = 0;

loop_flag = 0;

%% BASIC SETUP-DO NOT NEED TO EDIT GENERALLY
if date_flag
    base_date = datenum('10-Mar-2020','dd-mmm-yy');
else
    base_date = 0;
end

pcols = {1/255*[0,196,255],1/255*[225,12,62],1/255*[55,176,116],1/255*[244,201,107],1/255*[149,52,235],1/255*[176,176,176]};
symbs = {'o','x','.','+','s','v'};

curpath = pwd;
newpath = strcat(curpath(1:length(curpath)-6),'\data_figs\');
plotpath = strcat(curpath(1:length(curpath)-6),'\plots\vac_plots\');
runpath = strcat(curpath,'\runs\');

pos = [2,4.5,6,4.5];

exp_len = 1;
sub_len = [0];
%% 

%0 - plot variants and vaccines

done = 0;
experiment = 0;
sub_experiment = 0;

if loop_flag
    experiment =0;
    sub_experiment =0;
end

while ~done
    close all
    clear fig
    clr_flag = 0;
    if ~loop_flag
        done = 1;
    else
        fprintf('Currently running experiment %i sub experiment %i...\n',experiment,sub_experiment);
    end

    switch experiment
        case 0
            DATA = load('variant2_DATA-2021-Sept13.mat');
            DATA2 = load('Re_DATA-2021-Sept13.mat');
            DATA3 = load('vaccine_DATA-2021-Oct18.mat');
            [fig(1),figprop(1)]=new_figure('variants_plot',pos);
            
            hold on
%             for k=1:DATA.DATA_notes.varsize
%                 plot(base_date+DATA.DATA_T(14:end),DATA.var_avg14(:,k)./DATA.new_avg14,'Marker',symbs{k},'LineStyle','none')
%             end
            for k=1:DATA.DATA_notes.varsize
                plot(base_date+DATA.DATA_T,DATA.DATA_var(:,k)./max(sum(DATA.DATA_var'),1)','Marker',symbs{k},'LineStyle','none');
            end
            plot(base_date+DATA2.DATA_T,DATA2.DATA_Re,'p');
            hold off
            xlim([base_date+DATA.DATA_notes.time_var,base_date+length(DATA.DATA_T)]);
            if date_flag
                datetick('x','dd-mmm-yy','keepticks','keeplimits');
                xtickangle(45);
            else
                figprop(1).xlab=xlabel('$t$ (days since March 10, 2020)');
            end
            figprop(1).ylab=ylabel('Percent of cases');
            figprop(1).ax = gca;
%             figprop(1).leg = legend({'$\alpha$','$\beta$','$\gamma$','$\delta$','Re'},'location','northwest');
            figprop(1).leg = legend({'$\alpha$','$\beta$','$\gamma$','$\delta$','Re'},'location','northwest');
            figprop(1).fnt_size = 14;
            
            [fig(2),figprop(2)]=new_figure('vaccine_plot',pos);
            
            hold on
            plot(base_date+DATA3.DATA_T(1:end),DATA3.DATA_tot1,'o');
            plot(base_date+DATA3.DATA_T(1:end),DATA3.DATA_tot2,'x');
            hold off
            
            xlim([base_date+DATA.DATA_notes.time_var,base_date+length(DATA.DATA_T)]);
            if date_flag
                datetick('x','dd-mmm-yy','keepticks','keeplimits');
                xtickangle(45);
            else
                figprop(2).xlab=xlabel('$t$ (days since March 10, 2020)');
            end
            figprop(2).ylab=ylabel('Vaccine Doses');
            figprop(2).ax = gca;
            figprop(2).leg = legend({'1 dose','2 doses'},'location','northwest');
            figprop(2).fnt_size = 14;
            
            
            
            
            
    end
    
    if loop_flag
        if sub_len(experiment+1)>0
            sub_experiment = sub_experiment +1;
                if sub_experiment>sub_len(experiment+1)
                    sub_experiment = 0;
                    experiment = experiment +1;
                    clr_flag = 1;
                end
        else
            experiment = experiment +1;
            clr_flag = 1;
        end
        if clr_flag
            clearvars -except date_flag base_date loop_flag pcols curpath newpath plotpath runpath pos exp_len sub_len done experiment sub_experiment fig figprop
        end
        if experiment > exp_len
            done =1;
        end
    end
    figure_revise(pos,plotpath,fig,figprop);
end
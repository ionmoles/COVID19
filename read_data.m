close all
clear variables

year = '2021';
mday = 'Sept13';
data_flag = 2;
% 0 - covid cases
% 1 - vaccine data
% 2 - variant cases

switch data_flag
    case {0,2}
        basetitle='covidtesting';
    case 1
        basetitle='vaccine_doses';
end

curpath = pwd;
datpath = strcat(curpath(1:length(curpath)-6),'\Data\');
filetit = [datpath,year,'\',basetitle,' - ',mday,'.csv'];

T = readtable(filetit,'TreatAsEmpty','');

switch data_flag
    case 0
        DATA_notes.time={'first tested case Jan 28, 2020','first time point Mar 10, 2020'};
        start_date='2020-03-10';
        date_spot = find(T{:,1}==start_date);
        Tlen = size(T,1)-date_spot+1;
        DATA_T = linspace(0,Tlen-1,Tlen)';
        DATA_pos = T{date_spot:end,5};
        DATA_tot = T{date_spot:end,8};
    case 1
        DATA_notes.time={'first vaccine data Jan 7, 2021','first point is 302 days from March 10'};
        start_date='2021-01-06';
        date_spot = find(T{:,1}==start_date);
        Tlen = size(T,1)-date_spot+1;
        DATA_T = linspace(302,302+Tlen-1,Tlen)';
        DATA_tot1 = T{date_spot:end,6};
        DATA_totfull = T{date_spot:end,8};
    case 2
        T = fillmissing(T,'previous');
        DATA_notes.time={'variants first recorded Jan 29, 2021','first time point Mar 10, 2020','variants first measured 326 days after March 10'};
        DATA_notes.time_var = 326; %corresponds to variants starting 326 days after Mar 10
        DATA_notes.varsize=4;
        start_date='2020-03-10';
        date_spot = find(T{:,1}==start_date);
        variants_nan = isnan(T{:,24:24+DATA_notes.varsize-1});
        for k=1:DATA_notes.varsize
            T{variants_nan(:,k),24+(k-1)}=0;
        end
        Tlen = size(T,1)-date_spot+1;
        DATA_T = linspace(0,Tlen-1,Tlen)';
        DATA_tot = T{date_spot:end,8};
        DATA_new = T{date_spot:end,8}-T{date_spot-1:end-1,8};
        DATA_var = T{date_spot:end,24:24+DATA_notes.varsize-1}-T{date_spot-1:end-1,24:24+DATA_notes.varsize-1};
        new_avg = (DATA_new(1:end-6)+DATA_new(2:end-5)+DATA_new(3:end-4)+DATA_new(4:end-3)+DATA_new(5:end-2)+DATA_new(6:end-1)+DATA_new(7:end))/7;
        var_avg = (DATA_var(1:end-6,:)+DATA_var(2:end-5,:)+DATA_var(3:end-4,:)+DATA_var(4:end-3,:)+DATA_var(5:end-2,:)+DATA_var(6:end-1,:)+DATA_var(7:end,:))/7;
end

switch data_flag
    case 0
        save(['DATA-',year,'-',mday,'.mat'],'DATA_notes','DATA_T','DATA_pos','DATA_tot');
    case 1
        save(['vDATA-',year,'-',mday,'.mat'],'DATA_notes','DATA_T','DATA_tot1','DATA_totfull');
    case 2
        save(['variant_DATA-',year,'-',mday,'.mat'],'DATA_notes','DATA_T','DATA_tot','DATA_new','DATA_var','new_avg','var_avg');
end
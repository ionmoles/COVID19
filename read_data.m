close all
clear variables

year = '2021';
mday = 'Jul20';
data_flag = 1;
% 0 - covid cases
% 1 - vaccine data

switch data_flag
    case 0
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
end

switch data_flag
    case 0
        save(['DATA-',year,'-',mday,'.mat'],'DATA_notes','DATA_T','DATA_pos','DATA_tot');
    case 1
        save(['vDATA-',year,'-',mday,'.mat'],'DATA_notes','DATA_T','DATA_tot1','DATA_totfull');
end
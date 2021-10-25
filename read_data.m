close all
clear variables

year = '2021';
mday = 'Oct18';
data_flag = 0;
% 0 - covid cases
% 1 - vaccine data
% 2 - variant cases (from covidtesting)
% 3 - Re
% 4 - variant cases (from other PHO source)

switch data_flag
    case {0,2}
        basetitle='covidtesting';
    case 1
        basetitle='vaccine_doses';
    case 3
        basetitle='Re';
    case 4
        basetitle = 'variant';
end

curpath = pwd;
cd ..
uppath = cd;
cd (curpath);
datpath = strcat('..\Data\');
filetit = [datpath,year,'\',basetitle,' - ',mday,'.csv'];

try
    T = readtable(filetit,'TreatAsEmpty',{'','-'});
catch
    error('Matlab is trying to read %s, but file and/or directory do not exist.\n Make sure you have created the directory %s and have stored the correct CSV file.\n',strcat(uppath,'\Data\',year,'\',basetitle,' - ',mday,'.csv'),strcat(uppath,'\Data\',year,'\'));
end

switch data_flag
    case 0
        DATA_notes.time={'first tested case Jan 28, 2020','first time point Mar 10, 2020'};
        start_date='2020-03-10';
        date_spot = find(T{:,1}==start_date);
        Tlen = size(T,1)-date_spot+1;
        DATA_T = linspace(0,Tlen-1,Tlen)';
        DATA_pos = T{date_spot:end,5};
        DATA_tot = T{date_spot:end,8};
        DATA_test = T{date_spot:end,10};
        DATA_test(isnan(DATA_test)) = 0;
    case 1
        DATA_notes.time={'first vaccine data Jan 7, 2021','first point is 302 days from March 10'};
        start_date='2021-01-06';
        date_spot = find(T{:,1}==start_date);
        Tlen = size(T,1)-date_spot+1;
        DATA_T = linspace(302,302+Tlen-1,Tlen)';
        DATA_tot1 = T{date_spot:end,6};
        DATA_tot2 = T{date_spot:end,8};
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
        new_avg7 = 0;
        var_avg7 = 0;
        
        for k=1:7
            new_avg7 = new_avg7 + DATA_new(k:end-(7-k));
            var_avg7 = var_avg7 + DATA_var(k:end-(7-k),:);
        end
        
        new_avg7 = new_avg7/7;
        var_avg7 = var_avg7/7;
        
        new_avg14 = 0;
        var_avg14 = 0;
        
        for k=1:14
            new_avg14 = new_avg14 + DATA_new(k:end-(14-k));
            var_avg14 = var_avg14 + DATA_var(k:end-(14-k),:);
        end
        
        new_avg14 = new_avg14/14;
        var_avg14 = var_avg14/14;
        
        new_avg21 = 0;
        var_avg21 = 0;
        
        for k=1:21
            new_avg21 = new_avg21 + DATA_new(k:end-(21-k));
            var_avg21 = var_avg21 + DATA_var(k:end-(21-k),:);
        end
        
        new_avg21 = new_avg21/21;
        var_avg21 = var_avg21/21;
        
%         new_avg = (DATA_new(1:end-6)+DATA_new(2:end-5)+DATA_new(3:end-4)+DATA_new(4:end-3)+DATA_new(5:end-2)+DATA_new(6:end-1)+DATA_new(7:end))/7;
%         var_avg = (DATA_var(1:end-6,:)+DATA_var(2:end-5,:)+DATA_var(3:end-4,:)+DATA_var(4:end-3,:)+DATA_var(5:end-2,:)+DATA_var(6:end-1,:)+DATA_var(7:end,:))/7;
    case 3
        DATA_notes.time={'Re first recorded March 19, 2021','first point is 9 days after Marh 10'};
        Tlen = size(T,1);
        DATA_T = linspace(9,9+Tlen-1,Tlen)';
        DATA_Re = T{:,2};
    case 4
        DATA_notes.time={'variants first recorded Dec 20, 2020','variants first measured 285 days after March 10'};
        DATA_notes.time_var = 285; %corresponds to variants starting 285 days after Mar 10
        DATA_notes.varsize=4;
        start_date='10-Mar-0020';
        date_spot = find(T{:,1}==start_date);
        variants_nan = isnan(T{:,16:2:16+2*(DATA_notes.varsize+2-1)});
        for k=1:DATA_notes.varsize+2
            T{variants_nan(:,k),16+2*(k-1)}=0;
        end
        Tlen = size(T,1)-date_spot+1;
        DATA_T = linspace(0,Tlen-1,Tlen)';
        DATA_var = T{date_spot:end,16:2:16+2*(DATA_notes.varsize+2-1)};
end

switch data_flag
    case 0
        save(['DATA-',year,'-',mday,'.mat'],'DATA_notes','DATA_T','DATA_pos','DATA_tot','DATA_test');
    case 1
        save(['vaccine_DATA-',year,'-',mday,'.mat'],'DATA_notes','DATA_T','DATA_tot1','DATA_tot2');
    case 2
        save(['variant_DATA-',year,'-',mday,'.mat'],'DATA_notes','DATA_T','DATA_tot','DATA_new','DATA_var','new_avg7','var_avg7','new_avg14','var_avg14','new_avg21','var_avg21');
    case 3
        save(['Re_DATA-',year,'-',mday,'.mat'],'DATA_notes','DATA_T','DATA_Re');
    case 4
        save(['variant2_DATA-',year,'-',mday,'.mat'],'DATA_notes','DATA_T','DATA_var');
end
%% APPLY FILTERS TO RAW DATA

% STEPS:
% - load raw stats, wake steering information, and list of valid time (qctime)
% - calculate turbulence intensity and append to channels
% - remove analog signals
% - statsfilt1: filtered by qctime, wind direction filter, turbine status, min power requirement
% - statsfilt2 on/off: filtered statsfilt1 with wake steering on/off
% - create binned database of each stats database

clear; close all; dbstop if error;

%% PROCESSING

% load raw stats data
load('results\T3_stats_raw.mat')
stats = t3_stats;

% load qc wake information
qcwaketable = readtable('qc_wake.csv');
%qcwake = table2array(qcwaketable(:,2:end));
qcwake = table2array(qcwaketable(4394:end,2:10)); %% HARDCODED!!!!
% load qc time information
qctimetable = readtable('qc_time.csv');
qctime = table2array(qctimetable(:,1));

fns = fieldnames(stats);
statslist = ["mean10min","max10min","min10min","stdev10min"];

% calculate turbulence intensity & add as field in struct
TI = stats.T3_Wind_Speed.stdev10min ./ stats.T3_Wind_Speed.mean10min;
stats.T3_TI.mean10min = TI;
stats.T3_TI.max10min = TI;
stats.T3_TI.min10min = TI;
stats.T3_TI.stdev10min = TI;

% remove analog fields
stats = rmfield(stats,fns(end-3:end));
stats = rmfield(stats,["AnalogIN_28","AnalogIN_29","AnalogIN_30","AnalogIN_31","AnalogIN_19"]);

fns = fieldnames(stats);

% check again if all variable stats are the same length
slength = length(t3_stats.(fns{1}).mean10min);
for k = 1:length(fns)
    for x = 1:length(statslist)
        if ~isequal(slength,length(stats.(fns{k}).(statslist{x})))
            fprintf('WARNING: variable length differs: %s-%s \n',fns{k},statslist{x});
        end
    end
end

% remove time period where blades were not installed (hardcoded!!!!)
for k = 1:length(fns)
    for x = 1:length(statslist)
        stats.(fns{k}).(statslist{x}) = stats.(fns{k}).(statslist{x})(4394:end);
    end
end

%% APPLY FILTERS

% filter for data qc process
rawtime = datetime(stats.MS_Excel_Timestamp.mean10min,'convertfrom','excel');
rawtime.Second = round(rawtime.Second);
indx_qc = find(ismember(rawtime,qctime));

% filter for wind direction
filter_WD = [320;350]; % in line with T2 is about 325 deg
%indx_WD = find(stats.T3_WindDir.mean10min > filter_WD(1) & stats.T3_WindDir.mean10min < filter_WD(2));
indx_WD2 = find(qcwake(:,7) > filter_WD(1) & qcwake(:,7) < filter_WD(2) & qcwake(:,4)==2); % use this for a more proper comparison
%indx_WD_all = intersect(indx_WD,indx_WD2);

% status signal = 2
filter_status = [2.05;2.15];
indx_status = find(stats.T3_Turbine_Status.mean10min > filter_status(1) & stats.T3_Turbine_Status.mean10min < filter_status(2));

% apply min power filter here
filter_power = 0;
indx_power = find(stats.T3_Active_Power.min10min > filter_power);

% get final index for filter1
indx_filt1 = intersect(indx_status,intersect(indx_WD2,intersect(indx_power,indx_qc)));

% add extra filter for MS torque above zero
% remove any outlier files

% apply filt1
for x = 1:length(fns)
    for y = 1:length(statslist)
        stats_filt1.(fns{x}).(statslist{y}) = stats.(fns{x}).(statslist{y})(indx_filt1);
    end
end

% get wake steering index
indx_wkqcOn = find(qcwake(:,1)==1 & qcwake(:,4)==2 & qcwake(:,5)==2);
indx_wkqcOff = find(qcwake(:,1)==0 & qcwake(:,4)==2 & qcwake(:,5)==2);

indx_wkOn = intersect(indx_wkqcOn,indx_filt1);
indx_wkOff = intersect(indx_wkqcOff,indx_filt1);

% apply wake steering filters
statslist = ["mean10min","max10min","min10min","stdev10min"];
for x = 1:length(fns)
    for y = 1:length(statslist)
        stats_filt2on.(fns{x}).(statslist{y}) = stats.(fns{x}).(statslist{y})(indx_wkOn);
        stats_filt2off.(fns{x}).(statslist{y}) = stats.(fns{x}).(statslist{y})(indx_wkOff);
    end
end


% get wake steering actual vs commanded yaw offsets
wake_actual = 180-qcwake(:,2);
wake_control = qcwake(:,3);

% find indices for each yaw offset bin and create new database
yaw_angles = 1.25:2.5:21.25;
for x = 1:length(fns)
    for i = 1:length(yaw_angles)-1
        indx = find(wake_actual >= yaw_angles(i) & wake_actual < yaw_angles(i+1));
        indx2 = intersect(indx_wkOn,indx);
        for y = 1:length(statslist)
            stats_filt3.(fns{x}).(strcat('angle',num2str(i))).(statslist{y}) = stats.(fns{x}).(statslist{y})(indx2);
        end
    end
end


%% BIN DATA

% set bin edges for wind speed
edges = 0:1:35;

% loop through all the variables and store binned data
stats_binraw = binStats(fns,edges,stats);

% bin stats_filter1 data
stats_binfilt1 = binStats(fns,edges,stats_filt1);

% bin stats_filter2 data
stats_binfilt2on = binStats(fns,edges,stats_filt2on);
stats_binfilt2off = binStats(fns,edges,stats_filt2off);

% bin stats_filter3 data
stats_binfilt3 = binStatsAngles(fns,edges,stats_filt3);

% create new database where binned statsfilt2 is normalized by baseline wind speed
% bin @ 14.5 m/s 
for x = 3:length(fns) % start after timestamp variables
    for y = 1:length(statslist)-1
        stats_binfiltnormOn.(fns{x}).(statslist{y})(:,1) = stats_binfilt2on.(fns{x}).(statslist{y})(:,1);
        stats_binfiltnormOn.(fns{x}).(statslist{y})(:,2) = stats_binfilt2on.(fns{x}).(statslist{y})(:,2);
        stats_binfiltnormOn.(fns{x}).(statslist{y})(:,3) = stats_binfilt2on.(fns{x}).(statslist{y})(:,3)./max(abs(stats_binfilt2off.(fns{x}).mean10min(1:16,3)));
        stats_binfiltnormOn.(fns{x}).(statslist{y})(:,4) = stats_binfilt2on.(fns{x}).(statslist{y})(:,4)./max(abs(stats_binfilt2off.(fns{x}).mean10min(1:16,3)));
        stats_binfiltnormOff.(fns{x}).(statslist{y})(:,1) = stats_binfilt2off.(fns{x}).(statslist{y})(:,1);
        stats_binfiltnormOff.(fns{x}).(statslist{y})(:,2) = stats_binfilt2off.(fns{x}).(statslist{y})(:,2);
        stats_binfiltnormOff.(fns{x}).(statslist{y})(:,3) = stats_binfilt2off.(fns{x}).(statslist{y})(:,3)./max(abs(stats_binfilt2off.(fns{x}).mean10min(1:16,3)));
        stats_binfiltnormOff.(fns{x}).(statslist{y})(:,4) = stats_binfilt2off.(fns{x}).(statslist{y})(:,4)./max(abs(stats_binfilt2off.(fns{x}).mean10min(1:16,3)));
    end
end

% normalize binned DELs by baseline @ 14.5 m/s
load('Mlife\DELs\binDELoff.mat')
load('Mlife\DELs\binDELon.mat')

DELfns = fieldnames(binDELoff); % channels from mlife
for x = 1:length(DELfns)
    % store data in order [wind speed,mean,std,count]
    stats_DELnormOn.(DELfns{x})(:,1) = binDELon.(DELfns{x})(:,1);
    stats_DELnormOn.(DELfns{x})(:,2) = binDELon.(DELfns{x})(:,2)./max(binDELoff.(DELfns{x})(1:15,2));
    stats_DELnormOn.(DELfns{x})(:,3) = binDELon.(DELfns{x})(:,3)./max(binDELoff.(DELfns{x})(1:15,2));
    stats_DELnormOn.(DELfns{x})(:,4) = binDELon.(DELfns{x})(:,4);
    stats_DELnormOff.(DELfns{x})(:,1) = binDELoff.(DELfns{x})(:,1);
    stats_DELnormOff.(DELfns{x})(:,2) = binDELoff.(DELfns{x})(:,2)./max(binDELoff.(DELfns{x})(1:15,2));
    stats_DELnormOff.(DELfns{x})(:,3) = binDELoff.(DELfns{x})(:,3)./max(binDELoff.(DELfns{x})(1:15,2));
    stats_DELnormOff.(DELfns{x})(:,4) = binDELoff.(DELfns{x})(:,4);
end


%% produce mlife files
% test1 = dir('..\FastDataMATLAB\T3_FT_formatted');
% test2 = struct2cell(test1);
% test2 = test2';
% test3 = test2(4394+2:end,1);
% listfiles = test3(indx_filt1);
%save('Mlife\listfiles_filt2on.mat','listfiles')

%% save data?
% save('results\withAllLoads\t3_stats_filt1.mat','stats_filt1');
% save('results\withAllLoads\t3_stats_filt2.mat','stats_filt2on','stats_filt2off');
% save('results\withAllLoads\t3_stats_filt3.mat','stats_filt3');
% save('results\withAllLoads\t3_binned_raw.mat','stats_binraw');
% save('results\withAllLoads\t3_binned_filt1.mat','stats_binfilt1');
% save('results\withAllLoads\t3_binned_filt2.mat','stats_binfilt2on','stats_binfilt2off');
% save('results\withAllLoads\t3_binned_filt3.mat','stats_binfilt3');
% save('results\withAllLoads\t3_binned_filtnormOn.mat','stats_binfiltnormOn');
% save('results\withAllLoads\t3_binned_filtnormOff.mat','stats_binfiltnormOff');
save('Mlife\DELs\binDELnormOn.mat','stats_DELnormOn');
save('Mlife\DELs\binDELnormOff.mat','stats_DELnormOff');

%% functions
function [binned] = binStats(fns,edges,raw)
for z = 3:length(fns)
    [a,b,c,d] = binning(raw.T3_Wind_Speed.mean10min,edges,raw.(fns{z}).mean10min); 
    binned.(fns{z}).mean10min = [a,b,c,d];
    [a,b,c,d] = binning(raw.T3_Wind_Speed.mean10min,edges,raw.(fns{z}).max10min);
    binned.(fns{z}).max10min = [a,b,c,d];
    [a,b,c,d] = binning(raw.T3_Wind_Speed.mean10min,edges,raw.(fns{z}).min10min);
    binned.(fns{z}).min10min = [a,b,c,d];
end
end

% bin yaw angles
function [binned] = binStatsAngles(fns,edges,raw)
anglesF = fieldnames(raw.(fns{3}));
for z = 3:length(fns)
    for x = 1:length(anglesF)
        [a,b,c,d] = binning(raw.T3_Wind_Speed.(anglesF{x}).mean10min,edges,raw.(fns{z}).(anglesF{x}).mean10min); 
        binned.(fns{z}).(anglesF{x}).mean10min = [a,b,c,d];
        [a,b,c,d] = binning(raw.T3_Wind_Speed.(anglesF{x}).mean10min,edges,raw.(fns{z}).(anglesF{x}).max10min);
        binned.(fns{z}).(anglesF{x}).max10min = [a,b,c,d];
        [a,b,c,d] = binning(raw.T3_Wind_Speed.(anglesF{x}).mean10min,edges,raw.(fns{z}).(anglesF{x}).min10min);
        binned.(fns{z}).(anglesF{x}).min10min = [a,b,c,d];
    end
end
end



% % create new database where binned statsfilt2on is normalized by statsfilt2off
% for x = 3:length(fns) % start after timestamp variables
%     for y = 1:length(statslist)-1
%         stats_binfiltnorm.(fns{x}).(statslist{y})(:,1) = stats_binfilt2on.(fns{x}).(statslist{y})(:,1);
%         stats_binfiltnorm.(fns{x}).(statslist{y})(:,2) = stats_binfilt2on.(fns{x}).(statslist{y})(:,2)./stats_binfilt2off.(fns{x}).(statslist{y})(:,2);
%         stats_binfiltnorm.(fns{x}).(statslist{y})(:,3) = stats_binfilt2on.(fns{x}).(statslist{y})(:,3)./stats_binfilt2off.(fns{x}).(statslist{y})(:,3);
%         stats_binfiltnorm.(fns{x}).(statslist{y})(:,4) = stats_binfilt2on.(fns{x}).(statslist{y})(:,4)./stats_binfilt2off.(fns{x}).(statslist{y})(:,4);
%     end
% end

%% LOAD DATA

% NOTES:

% PROCESSING STEPS:
% - Nacelle yaw and Nacelle_WD signals are unwrapped
% - Rotor azimuth data points above 360 are forced to 360 deg
% - Rotor azimuth data points below 0 are forced to 0 deg
% - Apply slopes and offsets to strain gage signals
% - Derive wind direction from nacelle vane and yaw position
% - Apply blade matrix coefficients
% - Apply coordinate transformation
% - Test if calcstats produce nans and skip file if that is true
% - Store scaled data and stats into each file struct
% - Test if length of all stat variables is equal
% - Store raw stats and produce csv of timestamps

clear; close all; clc; dbstop if error;
%% INPUTS

% start diary
diary 'T3_run_ReRun2.txt'

% path to raw files
path_fast = '..\FastDataMATLAB\T3\';

% start and end date of analysis
% date_start = "2019-11-08";
date_start = "2019-11-08";
date_end = "2020-02-16";

% list of statistic fields
statslist = ["mean10min","max10min","min10min","stdev10min"];

% coordinate transformation info
yawname = 'T3_Nacelle_Yaw_Position'; % name of yaw channel
headingTB = 285.6; % positive heading for TB bridge 1
headingTT = 277.4; % positive heading for TT bridge 1

% wind speed edges used for binning
edges = 0:1:35;

%% PROCESSING

% IMPORT slopes and offsets
opts = spreadsheetImportOptions("NumVariables", 3);
% Specify sheet and range
opts.Sheet = "ImportPage";
opts.DataRange = "A2:C22";
% Specify column names and types
opts.VariableNames = ["Variable", "Slope", "Offset"];
opts.VariableTypes = ["string", "double", "double"];
opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
opts = setvaropts(opts, 1, "EmptyFieldRule", "auto");
% Import the data
cal_data = readtable("T3_Slopes_and_Offsets_Final.xlsx", opts, "UseExcel", false);
clear opts

% get list of days between start and end date of analysis
list_day = dir(path_fast);
list_day = struct2cell(list_day);
list_day = list_day(1,3:end);
indx = find(datetime(list_day)>=datetime(date_start) & datetime(list_day)<=datetime(date_end));
list_day = list_day(indx);

% start looping through list_days
for i = 1:length(list_day)
    list_files = dir(fullfile(strcat(path_fast,char(list_day(i))), '*.mat')); % get list of files in folder   
    % loop through all the files in the day folder
    for j = 1:length({list_files.name})
        load(strcat(path_fast,char(list_day(i)),'\',list_files(j).name));
        % loop through variables to process data
        temp = fieldnames(data_out);
        data_out.chanlist = temp(3:end);
        fns = data_out.chanlist;
        for k = 1:length(fns)
            if isequal(fns{k},'T3_Nacelle_Yaw_Position') | isequal(fns{k},'T3_Nacelle_WD')
                data_out.(fns{k}).scaled_data = vectorunwrapd(data_out.(fns{k}).data); % unwrap vector
                data_out.(fns{k}).avtype = 'vector';
            elseif isequal(fns{k},'T3_Rotor_Azimuth') % correct erroneous data points
                data_out.(fns{k}).scaled_data = data_out.(fns{k}).data;
                i_azA = find(data_out.(fns{k}).scaled_data>360);
                data_out.(fns{k}).scaled_data(i_azA) = 360;
                i_azB = find(data_out.(fns{k}).scaled_data<0);
                data_out.(fns{k}).scaled_data(i_azB) = 0;
                data_out.(fns{k}).avtype = 'vector';
            elseif isequal(fns{k},'MS_Excel_Timestamp') | isequal(fns{k},'LabVIEW_Timestamp')
                data_out.(fns{k}).scaled_data = data_out.(fns{k}).data;
                data_out.(fns{k}).avtype = 'time';
            elseif isequal(fns{k},'T3_TB_Bending_1')
                data_out.(fns{k}).scaled_data = (data_out.(fns{k}).data-cal_data.Offset(5))*cal_data.Slope(5);
                data_out.(fns{k}).avtype = 'scalar';
            elseif isequal(fns{k},'T3_TB_Bending_2')
                data_out.(fns{k}).scaled_data = (data_out.(fns{k}).data-cal_data.Offset(6))*cal_data.Slope(6);
                data_out.(fns{k}).avtype = 'scalar';
            elseif isequal(fns{k},'T3_TT_Bending_1')
                data_out.(fns{k}).scaled_data = (data_out.(fns{k}).data-cal_data.Offset(7))*cal_data.Slope(7);
                data_out.(fns{k}).avtype = 'scalar';
            elseif isequal(fns{k},'T3_TT_Bending_2')
                data_out.(fns{k}).scaled_data = (data_out.(fns{k}).data-cal_data.Offset(8))*cal_data.Slope(8);
                data_out.(fns{k}).avtype = 'scalar';
            elseif isequal(fns{k},'T3_TT_Torque')
                data_out.(fns{k}).scaled_data = (data_out.(fns{k}).data-cal_data.Offset(9))*cal_data.Slope(9);
                data_out.(fns{k}).avtype = 'scalar';
            elseif isequal(fns{k},'T3_MS_Bending_1')
                data_out.(fns{k}).scaled_data = (data_out.(fns{k}).data-cal_data.Offset(1))*cal_data.Slope(1);
                data_out.(fns{k}).avtype = 'scalar';
            elseif isequal(fns{k},'T3_MS_Bending_2')
                data_out.(fns{k}).scaled_data = (data_out.(fns{k}).data-cal_data.Offset(2))*cal_data.Slope(2);
                data_out.(fns{k}).avtype = 'scalar';
            elseif isequal(fns{k},'T3_MS_Torque')
                data_out.(fns{k}).scaled_data = (data_out.(fns{k}).data-cal_data.Offset(3))*cal_data.Slope(3);
                data_out.(fns{k}).avtype = 'scalar';
            elseif isequal(fns{k},'T3_B1_Flap')
                data_out.(fns{k}).offset_data = data_out.(fns{k}).data - cal_data.Offset(10);
                data_out.(fns{k}).avtype = 'scalar';
            elseif isequal(fns{k},'T3_B1_Edge')
                data_out.(fns{k}).offset_data = data_out.(fns{k}).data - cal_data.Offset(11);
                data_out.(fns{k}).avtype = 'scalar';
            elseif isequal(fns{k},'T3_B2_Flap')
                data_out.(fns{k}).offset_data = data_out.(fns{k}).data - cal_data.Offset(16);
                data_out.(fns{k}).avtype = 'scalar';
            elseif isequal(fns{k},'T3_B2_Edge')
                data_out.(fns{k}).offset_data = data_out.(fns{k}).data - cal_data.Offset(17);
                data_out.(fns{k}).avtype = 'scalar';
            else
                data_out.(fns{k}).scaled_data = data_out.(fns{k}).data;
                data_out.(fns{k}).avtype = 'scalar';
            end                      
        end
        %% derive wind direction
        data_out.T3_WindDir.scaled_data = data_out.T3_Nacelle_WD.scaled_data - 280 + data_out.T3_Nacelle_Yaw_Position.scaled_data;
        data_out.T3_WindDir.avtype = 'vector';
        
        %% BLADE MATRIX APPLICATION
        if isfield(data_out,'T3_B1_Flap')
            % blade 1
            data_out.T3_B1_Flap.scaled_data = cal_data.Slope(12)*data_out.T3_B1_Flap.offset_data + cal_data.Slope(13)*data_out.T3_B1_Edge.offset_data;
            data_out.T3_B1_Edge.scaled_data = cal_data.Slope(14)*data_out.T3_B1_Flap.offset_data + cal_data.Slope(15)*data_out.T3_B1_Edge.offset_data;
            % blade 2
            data_out.T3_B2_Flap.scaled_data = cal_data.Slope(18)*data_out.T3_B2_Flap.offset_data + cal_data.Slope(19)*data_out.T3_B2_Edge.offset_data;
            data_out.T3_B2_Edge.scaled_data = cal_data.Slope(20)*data_out.T3_B2_Flap.offset_data + cal_data.Slope(21)*data_out.T3_B2_Edge.offset_data;
        end
        %% COORDINATE TRANSFORMATION
        % tower base
        yaw = data_out.(yawname).scaled_data;
        TB_angles = (yaw-headingTB).*pi/180;
        data_out.T3_TB_momentFA.scaled_data = data_out.T3_TB_Bending_1.scaled_data.*cos(TB_angles) + data_out.T3_TB_Bending_2.scaled_data.*sin(TB_angles);
        data_out.T3_TB_momentSS.scaled_data = data_out.T3_TB_Bending_1.scaled_data.*sin(TB_angles) - data_out.T3_TB_Bending_2.scaled_data.*cos(TB_angles);
        data_out.T3_TB_momentFA.avtype = 'scalar';
        data_out.T3_TB_momentSS.avtype = 'scalar';
        % tower top
        TT_angles = (yaw-headingTT).*pi/180;
        data_out.T3_TT_momentFA.scaled_data = data_out.T3_TT_Bending_1.scaled_data.*cos(TT_angles) - data_out.T3_TT_Bending_2.scaled_data.*sin(TT_angles);
        data_out.T3_TT_momentSS.scaled_data = data_out.T3_TT_Bending_1.scaled_data.*sin(TT_angles) + data_out.T3_TT_Bending_2.scaled_data.*cos(TT_angles);
        data_out.T3_TT_momentFA.avtype = 'scalar';
        data_out.T3_TT_momentSS.avtype = 'scalar';
        % mainshaft
        MS_angles = (data_out.T3_Rotor_Azimuth.scaled_data-180).*pi/180;
        data_out.T3_MS_yaw.scaled_data = data_out.T3_MS_Bending_1.scaled_data.*sin(MS_angles) - data_out.T3_MS_Bending_2.scaled_data.*cos(MS_angles);
        data_out.T3_MS_pitch.scaled_data = -data_out.T3_MS_Bending_1.scaled_data.*cos(MS_angles) - data_out.T3_MS_Bending_2.scaled_data.*sin(MS_angles);
        data_out.T3_MS_yaw.avtype = 'scalar';
        data_out.T3_MS_pitch.avtype = 'scalar';
        % blades
        fns = [fns;'T3_TB_momentFA';'T3_TB_momentSS';'T3_TT_momentFA';'T3_TT_momentSS';'T3_MS_yaw';'T3_MS_pitch';'T3_WindDir'];
        data_out.chanlist = fns;
        %% Calculate stats and store
        data_out = calcstats(data_out,'10min'); % calc 10min stats
        % check if every channel has a valid statistic (no NaNs)...i had to add this because i noticed that the raw stats 
        % produced by this script had variables that were not equal length at a certain location. This
        % ensures that all variables match together.
        test_nan = 0;
        for k = 1:length(fns)
            if isnan(data_out.(fns{k}).mean10min)
                test_nan = 1; 
            end
        end
        if test_nan==1
            clear test_nan
            fprintf('Skipped: %s \n',list_files(j).name)
            continue % skip iteration if statsfile has nan values
        end
        % create struct containing only stats data
        for k = 1:length(fns)
            for x = 1:length(statslist)
                try
                    t3_stats.(fns{k}).(statslist{x}) = [t3_stats.(fns{k}).(statslist{x}), data_out.(fns{k}).(statslist{x})];
                catch
                    t3_stats.(fns{k}).(statslist{x}) = data_out.(fns{k}).(statslist{x});
                end
            end
        end
        % add nans to blade data if not existent
        if ~isfield(data_out,'T3_B1_Flap')
            for x = 1:length(statslist)
                try
                    t3_stats.T3_B1_Flap.(statslist{x}) = [t3_stats.T3_B1_Flap.(statslist{x}), NaN];
                    t3_stats.T3_B1_Edge.(statslist{x}) = [t3_stats.T3_B1_Edge.(statslist{x}), NaN];
                    t3_stats.T3_B2_Flap.(statslist{x}) = [t3_stats.T3_B2_Flap.(statslist{x}), NaN];
                    t3_stats.T3_B2_Edge.(statslist{x}) = [t3_stats.T3_B2_Edge.(statslist{x}), NaN];
                catch
                    t3_stats.T3_B1_Flap.(statslist{x}) = NaN;
                    t3_stats.T3_B1_Edge.(statslist{x}) = NaN;
                    t3_stats.T3_B2_Flap.(statslist{x}) = NaN;
                    t3_stats.T3_B2_Edge.(statslist{x}) = NaN;
                end
            end
        end
        % save updated data_out struct
        save(strcat('..\FastDataMATLAB\T3_FT_formatted\',list_files(j).name),'data_out');
    end
    fprintf('Finished: %s \n',list_day{i})
end

% check if all variable stats are the same length
slength = length(t3_stats.(fns{1}).mean10min);
for k = 1:length(fns)
    for x = 1:length(statslist)
        if ~isequal(slength,length(t3_stats.(fns{k}).(statslist{x})))
            fprintf('WARNING: variable length differs: %s-%s \n',fns{k},statslist{x});
        end
    end
end

% save raw stats
save('results\T3_stats_raw.mat','t3_stats');
% save raw excel datetimes
%dlmwrite('T3_raw_stat_timestamps.csv',t3_stats.MS_Excel_Timestamp.mean10min,'precision',10);

diary off

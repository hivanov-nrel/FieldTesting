%% PLOTTING 

clear; dbstop if error; close all

%% LOAD DATA

% select database
withX = 'withAllLoads';

% load in data
load('results\T3_stats_raw.mat');
load(strcat('results\',withX,'\t3_stats_filt1.mat'));
load(strcat('results\',withX,'\t3_stats_filt2.mat'));
load(strcat('results\',withX,'\t3_stats_filt3.mat'));
load(strcat('results\',withX,'\t3_binned_raw.mat'));
load(strcat('results\',withX,'\t3_binned_filt1.mat'));
load(strcat('results\',withX,'\t3_binned_filt2.mat'));
load(strcat('results\',withX,'\t3_binned_filt3.mat'));
%load(strcat('results\',withX,'\t3_binned_filtnorm.mat'));
load(strcat('results\',withX,'\t3_binned_filtnormOn.mat'));
load(strcat('results\',withX,'\t3_binned_filtnormOff.mat'));
load('Mlife\DELs\DELs_filt2off.mat')
load('Mlife\DELs\DELs_filt2on.mat')
load('Mlife\DELs\binDELoff.mat')
load('Mlife\DELs\binDELon.mat')
load('Mlife\DELs\binDELnormOff.mat')
load('Mlife\DELs\binDELnormOn.mat')

% common variables
fns = fieldnames(stats_filt1);
statslist = ["mean10min","max10min","min10min"];
anglesF = fieldnames(stats_filt3.(fns{3}));
angles = 2.5:2.5:20;
windspeed = stats_binfilt3.T3_Active_Power.angle1.mean10min(2:16,1);

% eliminate zero bins for cleaner plots
for x = 3:length(fns)
    for y = 1:length(statslist)
        stats_binraw.(fns{x}).(statslist{y}) = stats_binraw.(fns{x}).(statslist{y})(find(stats_binraw.(fns{x}).(statslist{y})(:,2)),:);
        stats_binfilt1.(fns{x}).(statslist{y}) = stats_binfilt1.(fns{x}).(statslist{y})(find(stats_binfilt1.(fns{x}).(statslist{y})(:,2)),:);
        stats_binfilt2off.(fns{x}).(statslist{y}) = stats_binfilt2off.(fns{x}).(statslist{y})(find(stats_binfilt2off.(fns{x}).(statslist{y})(:,2)),:);
        stats_binfilt2on.(fns{x}).(statslist{y}) = stats_binfilt2on.(fns{x}).(statslist{y})(find(stats_binfilt2on.(fns{x}).(statslist{y})(:,2)),:);
    end
end

%% WIND ROSE
%[figure_handle,~,~,~,~] = WindRose(t3_stats.T3_WindDir.mean10min,t3_stats.T3_Wind_Speed.mean10min,'AngleNorth',0,'AngleEast',90);

%% CAPTURE MATRIX

[matx1,TI1,WS1] = gencaptmat(stats_filt1,'T3_Wind_Speed');
matrixplot(matx1,TI1,WS1)
[matx2,TI2,WS2] = gencaptmat(stats_filt2on,'T3_Wind_Speed');
matrixplot(matx2,TI2,WS2,'Steered')
[matx3,TI3,WS3] = gencaptmat(stats_filt2off,'T3_Wind_Speed');
matrixplot(matx3,TI3,WS3,'Baseline')
%savefig('plots\CAPTURE_MATRIX.fig')

%% CAPTURE MATRIX - 2D
figure;
hold on
bar(stats_binfilt2off.T3_Wind_Speed.mean10min(:,1),stats_binfilt2off.T3_Active_Power.mean10min(:,2))
bar(stats_binfilt2on.T3_Wind_Speed.mean10min(:,1),stats_binfilt2on.T3_Active_Power.mean10min(:,2),'FaceAlpha',0.5)
xlabel('Wind speed bin center [m/s]');ylabel('Num of 10-minute statistics');grid on; legend('Baseline','Steered')
xlim([0 25])
hold off

%% CAPTURE MATRIX - YAW ANGLE????????

figure;
clear offset_data
for y = 1:length(anglesF)
    offset_data(:,y) = stats_binfilt3.T3_Wind_Speed.(anglesF{y}).max10min(5:16,2);
end
h = bar3(offset_data,1);
shading interp; colorbar;
for i = 1:length(h)
    zdata = get(h(i),'Zdata');
    set(h(i),'Cdata',zdata)
    set(h,'EdgeColor','k')
end
for ii = 1:size(offset_data,1)
  for jj = 1:size(offset_data,2)
      if offset_data(ii,jj) == 0
          numdisplay = '';
      else
          numdisplay = num2str(offset_data(ii,jj));
      end
    text(jj,ii,offset_data(ii,jj),numdisplay,'FontSize',12,'HorizontalAlignment','center','Color',[0,0,0],'Position',[jj ii offset_data(ii,jj)+0.1 ])%,'Color',[0.8,0.8,0.8]
  end
end
view(-90,90)
yticklabels(windspeed(4:end)); xticklabels(angles);
ylabel('Wind Speed [m/s]'); xlabel('Yaw Offset');
ylabel(colorbar,'Number of 10-min averages')

%% TURBULENCE INTENSITY

figure
hold on
scatter(stats_filt2off.T3_Wind_Speed.mean10min,stats_filt2off.T3_TI.mean10min*100)
scatter(stats_filt2on.T3_Wind_Speed.mean10min,stats_filt2on.T3_TI.mean10min*100)
hold off
grid on; xlabel('Wind Speed [m/s]'); ylabel('TI [%]')
legend('Baseline','Steered')

figure;
set(gca,'DefaultLineLineWidth',1.2);
hold on
errorbar(stats_binfilt2off.T3_TI.mean10min(:,1),stats_binfilt2off.T3_TI.mean10min(:,3)*100,stats_binfilt2off.T3_TI.mean10min(:,4)*100,'MarkerSize',1)
errorbar(stats_binfilt2on.T3_TI.mean10min(:,1),stats_binfilt2on.T3_TI.mean10min(:,3)*100,stats_binfilt2on.T3_TI.mean10min(:,4)*100,'MarkerSize',1)
xlabel('Wind Speed [m/s]'); ylabel('TI [%]'); grid on; legend('Baseline','Steered');
hold off

%% BINNED STATS NORMALIZED PLOTS

variab = ["T3_B1_Flap","T3_B1_Edge","T3_MS_pitch","T3_MS_yaw","T3_MS_Torque","T3_TT_momentFA","T3_TT_momentSS","T3_TT_Torque","T3_TB_momentFA","T3_TB_momentSS"];
label = ["Normalized Blade Flap Moment [-]","Normalized Blade Edge Moment [-]","Normalized LSS M_y [-]","Normalized LSS M_z [-]","Normalized LSS Torque [-]","Normalized Tower Top F-A Moment [-]","Normalized Tower Top S-S Moment [-]","Normalized Tower Top Torque [-]","Normalized Tower Base F-A Moment [-]","Normalized Tower Base S-S Moment [-]"];
for x = 1:length(variab)
    figure;
    hold on
    errorbar(stats_binfiltnormOff.(variab{x}).mean10min(:,1),stats_binfiltnormOff.(variab{x}).mean10min(:,3),stats_binfiltnormOff.(variab{x}).mean10min(:,4),'Color',[0, 0.4470, 0.7410]);
    errorbar(stats_binfiltnormOn.(variab{x}).mean10min(:,1),stats_binfiltnormOn.(variab{x}).mean10min(:,3),stats_binfiltnormOn.(variab{x}).mean10min(:,4),'Color',[0.8500, 0.3250, 0.0980]);
    errorbar(stats_binfiltnormOff.(variab{x}).mean10min(:,1),stats_binfiltnormOff.(variab{x}).max10min(:,3),stats_binfiltnormOff.(variab{x}).max10min(:,4),'-^','Color',[0, 0.4470, 0.7410]);
    errorbar(stats_binfiltnormOn.(variab{x}).mean10min(:,1),stats_binfiltnormOn.(variab{x}).max10min(:,3),stats_binfiltnormOn.(variab{x}).max10min(:,4),'-^','Color',[0.8500, 0.3250, 0.0980]);
    errorbar(stats_binfiltnormOff.(variab{x}).mean10min(:,1),stats_binfiltnormOff.(variab{x}).min10min(:,3),stats_binfiltnormOff.(variab{x}).min10min(:,4),'-v','Color',[0, 0.4470, 0.7410]);
    errorbar(stats_binfiltnormOn.(variab{x}).mean10min(:,1),stats_binfiltnormOn.(variab{x}).min10min(:,3),stats_binfiltnormOn.(variab{x}).min10min(:,4),'-v','Color',[0.8500, 0.3250, 0.0980]);
    hold off
    xlim([4 15]); grid on; xlabel('Wind Speed [m/s]'); ylabel(label{x});
    legend('Baseline','Steered','Location','best')
end

%% DEL NORM COMPARE

variab = ["T3_B1_Flap","T3_B1_Edge","T3_MS_pitch","T3_MS_yaw","T3_MS_torque","T3_TT_momentFA","T3_TT_momentSS","T3_TT_Torque","T3_TB_momentFA","T3_TB_momentSS"];
label = ["Normalized Blade Flap Moment [-]","Normalized Blade Edge Moment [-]","Normalized LSS M_y [-]","Normalized LSS M_z [-]","Normalized LSS Torque [-]","Normalized Tower Top F-A Moment [-]","Normalized Tower Top S-S Moment [-]","Normalized Tower Top Torque [-]","Normalized Tower Base F-A Moment [-]","Normalized Tower Base S-S Moment [-]"];
for x = 1:length(variab)
    binplot_DELs(stats_DELnormOff,stats_DELnormOn,variab{x},'Wind Speed',label{x})
end

%% WAKE STEERING PLOTS
% figure;
% hold on
% bar(stats_binfilt2off.T3_Wind_Speed.mean10min(:,1),stats_binfilt2off.T3_Active_Power.mean10min(:,3))
% bar(stats_binfilt2on.T3_Wind_Speed.mean10min(:,1),stats_binfilt2on.T3_Active_Power.mean10min(:,3),'FaceAlpha',0.5)
% xlabel('Wind speed');ylabel('Power');grid on; legend('baseline','on')
% hold off
% 
% figure;
% hold on
% bar(stats_binfilt2off.T3_Wind_Speed.mean10min(:,1),stats_binfilt2off.T3_TB_momentFA.mean10min(:,3))
% bar(stats_binfilt2on.T3_Wind_Speed.mean10min(:,1),stats_binfilt2on.T3_TB_momentFA.mean10min(:,3),'FaceAlpha',0.5)
% xlabel('Wind speed');ylabel('TB_FA');grid on; legend('baseline','on')
% hold off

%% PLOT SCATTER

% scatplot(stats_filt1,'T3_Wind_Speed','T3_Active_Power','Wind Speed','Power')
% scatplot(stats_filt1,'T3_Wind_Speed','T3_TB_momentFA','Wind Speed','TB FAMoment')
% scatplot(stats_filt1,'T3_Wind_Speed','T3_TB_momentSS','Wind Speed','TB SSMoment')
%scatplot(stats_filt1,'T3_Wind_Speed','T3_TT_momentFA','Wind Speed','TT FAMoment')
%scatplot(stats_filt1,'T3_Wind_Speed','T3_TT_momentSS','Wind Speed','TT SSMoment')
% scatplot(stats_filt1,'T3_Wind_Speed','T3_TT_Torque','Wind Speed','TT torque')
% scatplot(stats_filt1,'T3_Wind_Speed','T3_MS_yaw','Wind Speed','MS yaw')
% scatplot(stats_filt1,'T3_Wind_Speed','T3_MS_pitch','Wind Speed','MS pitch')
% scatplot(stats_filt1,'T3_Wind_Speed','T3_MS_Torque','Wind Speed','MS torque')
% scatplot(stats_filt1,'T3_Wind_Speed','T3_B1_Flap','Wind Speed','B1 flap')
% scatplot(stats_filt1,'T3_Wind_Speed','T3_B1_Edge','Wind Speed','B1 edge')

% scatplot(stats_filt1,'T3_Wind_Speed','T3_MS_Bending_2','Wind Speed','MS2')
%scatplot(stats_filt1,'T3_Wind_Speed','T3_WindDir','Wind Speed','Dir')

% scatplot(stats_filt2on,'T3_Wind_Speed','T3_MS_yaw','Wind Speed','MSyaw')
% scatplot(stats_filt2off,'T3_Wind_Speed','T3_MS_yaw','Wind Speed','MSyaw')

%% PLOT BINNED DATA

% binned data
%binplot(stats_binfilt1,'T3_Active_Power','Wind Speed','Power')
%binplot(stats_binfilt1,'T3_MS_pitch','Wind Speed','MSpitch')
%binplot(stats_binfilt1,'T3_B1_Flap','Wind Speed','B1_flap')
%binplot(stats_binfilt2on,'T3_TT_Torque','Wind Speed','tt_torque')
%binplot(stats_binfilt2off,'T3_TT_Torque','Wind Speed','tt_torque')
% binplot(stats_binfilt2on,'T3_TB_momentFA','Wind Speed','tb foreaft')
% binplot(stats_binfilt2off,'T3_TB_momentFA','Wind Speed','tb foreaft')

%% BINNED COMPARISON PLOTS
% binplot_compare(stats_binfilt2off,stats_binfilt2on,'T3_MS_yaw','wind speed','ms yaw')
binplot_compare(stats_binfilt2off,stats_binfilt2on,'T3_MS_pitch','wind speed','ms pitch')
%binplot_compare(stats_binfilt2off,stats_binfilt2on,'T3_TB_momentSS','wind speed','tb ss')
%binplot_compare(stats_binfilt2off,stats_binfilt2on,'T3_TT_momentFA','wind speed','tt FA')
%binplot_compare(stats_binfilt2off,stats_binfilt2on,'T3_TT_Torque','wind speed','tt torque')
%binplot_compare(stats_binfilt2off,stats_binfilt2on,'T3_Active_Power','wind speed','power')


%% DELS BINNED COMPARE
variab = 'T3_MS_pitch';
binplot_DELs(binDELoff,binDELon,variab,'Wind Speed',variab)

%% DEL SCATTER COMPARE
variab = 'T3_B1_Edge';
figure;
hold on
scatter(DEL2off.windspeed,DEL2off.(variab),'^')
scatter(DEL2on.windspeed,DEL2on.(variab),'+')
xlabel('Wind Speed'); ylabel(variab);
grid on

%% YAW ANGLE OFFSETS 3D
variab = "T3_MS_yaw";
clear data
for y = 1:length(anglesF)
    data(:,y) = abs(stats_binfilt3.(variab).(anglesF{y}).max10min(2:16,3));
end
bar3(data,1)
%xlim([2 8])
xticklabels(angles);
yticklabels(windspeed); % ylim([0 21]); yticklabels(int2str(windspeed(2:2:21)));
ylabel('Wind Speed Bin Center [m/s]'); xlabel('Yaw Offset [deg]'); %zlabel(vars{v},'Interpreter','none')

%% YAW ANGLE OFFSETS SCATTER
% scatter plots
for x = 1:length(angles)
    scatplot_angle(stats_filt3,'T3_Wind_Speed','T3_MS_yaw',x,'Wind Speed','MS_yaw')
end

%% YAW ANGLE OFFSETS BIN PLOTTER
colors = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560];[0.4660 0.6740 0.1880];[0.3010 0.7450 0.9330];[0.6350 0.0780 0.1840]];
variab = "T3_MS_yaw";
figure
hold on
for x=1:length(angles)-1
    s1(x) = plot(stats_binfilt3.(variab).(anglesF{x}).max10min(2:16,1),stats_binfilt3.(variab).(anglesF{x}).max10min(2:16,3));
    s1(x).Color = colors(x,:);
    s1(x).LineWidth = 1;
    s2(x) = plot(stats_binfilt3.(variab).(anglesF{x}).min10min(2:16,1),stats_binfilt3.(variab).(anglesF{x}).min10min(2:16,3));
    s2(x).Color = colors(x,:);
    s2(x).LineWidth = 1;
end
hold off
xlabel('Wind Speed Bin Center [m/s]'); ylabel(variab); grid on;
xticks(stats_binfilt3.T3_Active_Power.angle1.mean10min(:,1))
legend(s1(:),string(angles(1:7)),'Location','best')






%% plotting functions
% scatter plots
function scatplot(stats,x,y,xlab,ylab)
figure;
hold on
sz = 25;
scatter(stats.(x).mean10min,stats.(y).max10min,sz,'^')
scatter(stats.(x).mean10min,stats.(y).mean10min,sz)
scatter(stats.(x).mean10min,stats.(y).min10min,sz,'v')
scatter(stats.(x).mean10min,stats.(y).stdev10min,sz,'+')
xlabel(xlab); ylabel(ylab); grid on; legend('max','mean','min','stdev','Location','best');
hold off
end

function scatplot_angle(stats,x,y,anglenum,xlab,ylab)
figure;
hold on
sz = 25;
angles = 2.5:2.5:20;
anglesF = fieldnames(stats.T3_Active_Power);
angle = anglesF{anglenum};
scatter(stats.(x).(angle).mean10min,stats.(y).(angle).max10min,sz,'^')
scatter(stats.(x).(angle).mean10min,stats.(y).(angle).mean10min,sz)
scatter(stats.(x).(angle).mean10min,stats.(y).(angle).min10min,sz,'v')
scatter(stats.(x).(angle).mean10min,stats.(y).(angle).stdev10min,sz,'+')
xlabel(xlab); ylabel(ylab); grid on; legend('max','mean','min','stdev','Location','best');
xlim([2 14])
title(strcat(num2str(angles(anglenum)),' deg offset'));
hold off
end

% binned plot data
function binplot(databin,y,xlab,ylab)
figure;
set(gca,'DefaultLineLineWidth',1.2);
hold on
errorbar(databin.(y).mean10min(:,1),databin.(y).max10min(:,3),databin.(y).max10min(:,4),'MarkerSize',1)
errorbar(databin.(y).mean10min(:,1),databin.(y).mean10min(:,3),databin.(y).mean10min(:,4),'MarkerSize',1)
errorbar(databin.(y).mean10min(:,1),databin.(y).min10min(:,3),databin.(y).min10min(:,4),'MarkerSize',1)
xlabel(xlab); ylabel(ylab); grid on; legend('max','mean','min','Location','best');
hold off
end

% bin plot compare
function binplot_compare(databin1,databin2,y,xlab,ylab)
figure;
set(gca,'DefaultLineLineWidth',1.2);
hold on
% data1
errorbar(databin1.(y).mean10min(:,1),databin1.(y).max10min(:,3),databin1.(y).max10min(:,4),'MarkerSize',1,'Color',[0 0.4470 0.7410])
errorbar(databin1.(y).mean10min(:,1),databin1.(y).mean10min(:,3),databin1.(y).mean10min(:,4),'MarkerSize',1,'Color',[0 0.4470 0.7410],'HandleVisibility','off')
errorbar(databin1.(y).mean10min(:,1),databin1.(y).min10min(:,3),databin1.(y).min10min(:,4),'MarkerSize',1,'Color',[0 0.4470 0.7410],'HandleVisibility','off')
% data2
errorbar(databin2.(y).mean10min(:,1),databin2.(y).max10min(:,3),databin2.(y).max10min(:,4),'MarkerSize',1,'Color',[0.8500 0.3250 0.0980])
errorbar(databin2.(y).mean10min(:,1),databin2.(y).mean10min(:,3),databin2.(y).mean10min(:,4),'MarkerSize',1,'Color',[0.8500 0.3250 0.0980],'HandleVisibility','off')
errorbar(databin2.(y).mean10min(:,1),databin2.(y).min10min(:,3),databin2.(y).min10min(:,4),'MarkerSize',1,'Color',[0.8500 0.3250 0.0980],'HandleVisibility','off')
xlabel(xlab); ylabel(ylab); grid on; legend('Baseline','Steered');
xlim([3 15])
hold off
end

function binplot_DELs(databin1,databin2,y,xlab,ylab)
figure;
set(gca,'DefaultLineLineWidth',1.2);
hold on
% data1
errorbar(databin1.(y)(:,1),databin1.(y)(:,2),databin1.(y)(:,3),'MarkerSize',1,'Color',[0 0.4470 0.7410])
% data2
errorbar(databin2.(y)(:,1),databin2.(y)(:,2),databin2.(y)(:,3),'MarkerSize',1,'Color',[0.8500 0.3250 0.0980])
xlabel(xlab); ylabel(ylab); grid on; legend('Baseline','Steered','Location','best');
xlim([4 15])
hold off
end


% capture matrixplots
function matrixplot(matx,TI,WS,tlabel)
matx = matx';
figure;
h = bar3(matx,1);
shading interp; colorbar;
for i = 1:length(h)
    zdata = get(h(i),'Zdata');
    set(h(i),'Cdata',zdata)
    set(h,'EdgeColor','k')
end
for ii = 1:size(matx,1)
  for jj = 1:size(matx,2)
      if matx(ii,jj) == 0
          numdisplay = '';
      else
          numdisplay = num2str(matx(ii,jj));
      end
    text(jj,ii,matx(ii,jj),numdisplay,'FontSize',12,'HorizontalAlignment','center','Color',[0,0,0],'Position',[jj ii matx(ii,jj)+0.1 ])%,'Color',[0.8,0.8,0.8]
  end
end
view(-90,90)
yticks(0.5:1:25.5); yticklabels(WS); xticks(1.5:1:15.5); xticklabels(TI);
ylabel('Wind Speed [m/s]'); xlabel('Turbulence Intensity [%]');
ylabel(colorbar,'Number of 10-min averages')
if nargin > 3
    title(tlabel)
end
end
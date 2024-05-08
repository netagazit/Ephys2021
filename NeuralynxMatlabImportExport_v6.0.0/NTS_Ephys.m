clear all; close all; clc
%% Define the variables
   Variables.FileName='2020-02-21 SUBLAT13-3 Tetrode 2 Unit 1_File_03';
for SelectVariables=1:1    %% Variables we choose
Variables.LickingBout=true;
Variables.WindowSizeSec=30; % window of peri-event plot 
Variables.MinimalboutIntervalSec=2;
Variables.TimeLimit=15; % time in minutes-limits the timeframe of the analysis from start to the specified time
Variables.Tagged=false;
Variables.PlotSpikeSorting=true;
Variables.MinTimestampsJelly=0; %limits for unit/ttl timestamps image
Variables.MaxTimestampsJelly=900; %limits for unit/ttl timestamps image
Variables.MinTimestampsChow=500; %limits for unit/ttl timestamps image
Variables.MaxTimestampsChow=900; %limits for unit/ttl timestamps image
Variables.CorrelationBins=30; % in seconds
Variables.ClusterCut=0.5;
Variables.TimeOfInterest=30; % for burst analysis time in seconds. 
Variables.PlotOffset=true;
Variables.TakeOutSingels=true;
Variables.LimitAxis=false;
end
for OtherVariables=1:1 %% Variables to extract from the file name
% add libreries to path
% for Neta's pc:
addpath(genpath('C:\Users\netas\Documents\MATLAB\GitHub\Obesity\NeuralynxMatlabImportExport_v6.0.0'));
addpath(genpath('C:\Users\netas\Documents\MATLAB\Han tetrode\in-vivo_ephys-master\in-vivo_ephys-master\Field Trip Neuralynx'));
% for Neta's new pc:
addpath(genpath('C:\Users\netas\Documents\MATLAB\GitHub\Backup_Works\NeuralynxMatlabImportExport_v6.0.0'));
addpath(genpath('C:\Users\netas\Documents\MATLAB\GitHub\Backup_Works\Field Trip Neuralynx'));
addpath 'C:\Users\netas\Documents\MATLAB\Obesity\InVivoEphys'
% for NLX PC:
addpath(genpath('C:\Users\admin\Documents\MATLAB\Events_spikes_v.2.1'));
addpath(genpath('C:\Users\admin\Documents\MATLAB\Obesity\November 2021\Field Trip Neuralynx'));

% find the folder location
if exist('D:\CheetahData\NG\Data\Temp Storage')
Variables.ComputerDir='D:\CheetahData\NG\Data\Temp Storage';
elseif exist('D:\Single unit data and analysis\SUBLAT')
Variables.ComputerDir='D:\Single unit data and analysis\SUBLAT';
elseif exist('D:')
Variables.ComputerDir='D:';
elseif exist('E:')
Variables.ComputerDir='E:';
else % if you cant find any data drives, use the example recording on this pc
Variables.ComputerDir='C:\Users\netas\Documents\Lammel lab\Obesity paper\EXPERIMENTS\Single units\Data for testing';
end
%% Details
ExtractName=char(Variables.FileName);
Variables.TetrodeNumber=str2num(ExtractName(31));
Variables.UnitNumber=str2num(ExtractName(38));
Variables.FileNumber=str2num(ExtractName(end));
Variables.MouseName= ExtractName(12:21);
Variables.Date= ExtractName(1:10);
%% diet type
REGlist={'SUBLAT13-3', 'SUBLAT13-9','SUBLAT11-7','SUBLAT18-9','SUBLAT19-9','SUBLAT18-5','SUBLAT21-5','SUBLAT20-7'};
for u=1:length(REGlist)
    if Variables.MouseName== REGlist{u}
        Variables.DietType='REG';
        break
    else
       Variables.DietType='HFD'; 
    end
end
%% folder actual name
FileDir=dir([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName]);FileDir=FileDir(3:end);
try
FileDir=extractfield(FileDir,'name')';
catch
    disp(['requested folder not found']);
   
end
if isempty(FileDir)
    disp(['requested folder not found']);
end
Variables.FolderName=char(FileDir(Variables.FileNumber));
Variables.FileDir=FileDir; clear FileDir;
%% determine the condition
ConditionList={'Jelly 1', 'Jelly 2', 'Empty','Chow', '5ms_2hz', 'Laser_Burst-5ms_2hz','Laser_Single-5ms_2hz','Laser_Single-1ms_2hz','Laser_Single-5ms_1hz','Jelly - Exposure','Jelly - OFF','Jelly','jelly','chow','Empty','Laser'};
for u=1:length(ConditionList)
    if contains(Variables.FolderName,ConditionList{u})
  Variables.Condition=ConditionList{u};
  break
    end
end; clear u
%% determine TTL port to analyse
% determine TTL type to analyse (Lickometer / Laser)
if ~isempty(strfind(Variables.FolderName,'Laser'))
Variables.TTLType='Laser'; % Lickometer / Laser
else
Variables.TTLType='Lickometer'; % Lickometer / Laser
end
if contains(Variables.TTLType,'Lickometer')
Variables.TTLPort=11; % TTL signal for lickometer enters port 1
elseif contains(Variables.TTLType,'Laser')
Variables.TTLPort=10; % TTL signal for lickometer enters port 0
end 
end %% other variables:
%%
%% Run InVivoEphysShort function to get basic data from folder
FileFolder=[Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\',Variables.FolderName,'\TT',num2str(Variables.TetrodeNumber),'_s.ntt'];
EphysObj = InVivoEphysShort([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\',Variables.FolderName,'\TT',num2str(Variables.TetrodeNumber),'_s.ntt'],Variables);
EphysObj.Variables=Variables;
for ExtractVariablesFromData=1:1
%% figure out which TTL to take, if we find laser pulses we take laser TTL, otherwise its Lickometer
% check if selected TTL type exists
for TTLs=1:1
ExtractTTL={EphysObj.TTL.type};
for w=1:length(ExtractTTL)
    if ExtractTTL{w}==10
        disp('Laser TTLs Analysed')
        break
     elseif ExtractTTL{w}==11
    disp('Licking TTLs Analysed')
    break
    else
 disp('relevant TTLs not found')
    end
end
for i=1:length(EphysObj.TTL)
Test=EphysObj.TTL(i).type; 
if contains(Variables.TTLType,'Laser')
if Test==10
TTL_Line=i;
EphysObj.Variables.TTL_Line=i;
end
elseif contains(Variables.TTLType,'Lickometer')
if Test==11
TTL_Line=i;
EphysObj.Variables.TTL_Line=i;

end
end
end
end; clear TTLs
if exist('TTL_Line')~=0
end
end
%% get the cell activity 
% get cell timestamps relatively to first recorded timestamp
CellTimeStamp_us=((EphysObj.raw_data.TimeStamp(EphysObj.raw_data.CellNumber==EphysObj.Variables.UnitNumber))-double(EphysObj.raw_data.hdr.FirstTimeStamp));
% Calculate the time interval between consecutive timestamps
Cell_time_interval_us = diff(CellTimeStamp_us);
% Convert time interval to seconds
Cell_time_interval_s = double(Cell_time_interval_us) / 1000000; % 1 microsecond = 1e-6 seconds
%Calculate frequency in Hz
Cell_frequency_hz = 1 ./ Cell_time_interval_s;
%% get the Piezo TTL activity
% get TTL timestamps relatively to first recorded timestamp
TTLTimeStamp_us=(EphysObj.TTL(w).on-double(EphysObj.raw_data.hdr.FirstTimeStamp)); %timestamp in microseconds
% Calculate the time interval between consecutive timestamps
TTL_time_interval_us = diff(TTLTimeStamp_us);
% Convert time interval to seconds
TTL_time_interval_s = double(TTL_time_interval_us) / 1000000; % 1 microsecond = 1e-6 seconds
%Calculate frequency in Hz
TTL_frequency_hz = 1 ./ TTL_time_interval_s;TTL_frequency_hz=[0 TTL_frequency_hz];






%% Plot TTLs and cell activity on a timeline
Xaxis=[0:(EphysObj.raw_data.hdr.LastTimeStamp-EphysObj.raw_data.hdr.FirstTimeStamp)];
scatter(CellTimeStamp_us,[ones(1,length(CellTimeStamp_us))],'|');
hold on
scatter(TTLTimeStamp_us,[2*ones(1,length(TTLTimeStamp_us))],'|');
ylim([0 10])




%% create variable to save
UnitData.FilePath=[EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\',EphysObj.Variables.FolderName];
UnitData.MouseName=EphysObj.Variables.MouseName;
UnitData.Condition=EphysObj.Variables.Condition;
UnitData.UnitTimestampsZeroed=CellTimeStamp_us;
UnitData.TTLTimestampsZeroed=TTLTimeStamp_us;
UnitData.FolderName=EphysObj.Variables.FolderName;  

%% this code uses the following function: InVivoEphys.m , Field Trip Neuralynx, NeuralynxMatlabImportExport_v6.0.0, burst.m. 
 clearvars -Except Variables BadCount FileStructure ListOfFiles Count  %% Count=1 % ListOfFiles={} % clear all
% clear all
a='2021-10-29_SUBLAT18-5_Tetrode_1_Unit_1_File_03'; % 2, 3, 4, 5, 2020-03-09_SUBLAT11-7_Tetrode_3_Unit_1_File_3
Count=1;BadCount=1;if ~exist('ListOfFiles'); ListOfFiles={a};end
for AddPathAndDir=1:1%% add to pathc
% add libraries of New PC
addpath(genpath('C:\Users\netas\OneDrive\Documents\GitHub\Ephys2021\DeepLabCut'));
addpath(genpath('C:\Users\netas\OneDrive\Documents\GitHub\Ephys2021\Field Trip Neuralynx'));
addpath(genpath('C:\Users\netas\OneDrive\Documents\GitHub\Ephys2021\NeuralynxMatlabImportExport_v6.0.0'));
% add libraries of old PC
addpath(genpath('C:\Users\netas\Documents\GitHub\InVivo_Ephys\NeuralynxMatlabImportExport_v6.0.0'));
addpath(genpath('C:\Users\netas\Documents\GitHub\InVivo_Ephys\Field Trip Neuralynx'));
addpath(genpath('C:\Users\netas\Documents\GitHub\InVivo_Ephys\DeepLabCut'));

%% find files location
if exist('D:')
Variables.ComputerDir='D:';
else % if you cant find any data drives, use the example recording on this pc
Variables.ComputerDir='C:\Users\netas\Documents\Lammel lab\Obesity paper\EXPERIMENTS\Single units\Data for testing';
end
end; clear AddPathAndDir
CountGood=0;
CountBad=0;
for i=1:length(ListOfFiles)
    try
clear EphysObj
for GetFileName=1:1
% the following code is for running a lot of files one by one
RunBatch=true;
if RunBatch
if ~exist('Count')
Count=1;NameList={};CollectAllData=nan(100,100);AllEphysObjs={};CollectAllRsquared=[];CollectAllRsquaredNoZeros=[];
DataBase=[];LaserCollect=[];LickometerCollect=[];
end   
% namelist=NameList(:,2);
ExtractName=char(ListOfFiles(i));
Variables.TetrodeNumber=str2num(ExtractName(31));
Variables.UnitNumber=str2num(ExtractName(38));
Variables.FileNumber=str2num(ExtractName(end));
Variables.MouseName= ExtractName(12:21);
Variables.Date= ExtractName(1:10);
Variables.UnitName=char(ListOfFiles(i)); %old Jelly 3, new Jelly 4
Variables.ComputerDir='D:';
Variables.Tagged=true;
% end
end
end ; clear GetFileName
%% if files are above 9 insert manually
for FileNumber=str2double([ExtractName(45) ExtractName(46)]) 
for Vars=1:1
if ~exist('Variables')
Variables.Date=a(1:10);               %'2021-08-29'; 
Variables.MouseName= a(12:21)           %'SUBLAT19-9'; 	
Variables.TetrodeNumber=str2double(a(31))               %[4];
Variables.UnitNumber=str2double(a(38))%[1];
Variables.FileNumber=FileNumber; %old Jelly 3, new Jelly 4
Variables.UnitName=a; %old Jelly 3, new Jelly 4
end
end
for SelectVars=1:1
%% other variables:
% for scatter plots
Variables.PlotScatter=true;
Variables.PlotOnset=true;
Variables.WindowSizeSec=25; % window of 10 sec peri-event plot (0.1 for laser, 20 for food)
Variables.RasterAll =false;%[1:8 10]
Variables.PlotZscore=false;
Variables.SortByMotifLength=false;
Variables.DisplayPlot=true;
Variables.UseAnova=true;


%% VARIABLES RELATED TO ANALYSIS
Variables.MinimalboutIntervalSec=4; % interval between bouts
Variables.BoutLengthSecondsLimit=0.5; % minimal bout length
Variables.TimeLimit=15; % time in minutes-limits the timeframe of the analysis from start to the specified time
%% Plot related variables 
%% VARIABLES THAT WE CAN EXTRACT FROM FILE NAME
REGlist={'SUBLAT13-3', 'SUBLAT13-9','SUBLAT11-7','SUBLAT18-9','SUBLAT19-9','SUBLAT18-5','SUBLAT21-5','SUBLAT20-7'};
for u=1:length(REGlist)
    if Variables.MouseName== REGlist{u}
        Variables.DietType='REG';
        break
    else
       Variables.DietType='HFD'; 
    end
end
%% Extract files from mouse name folder
FileDir=dir([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName]);FileDir=FileDir(3:end);
try
FileDir=extractfield(FileDir,'name')';
catch
    disp(['requested folder not found']); 
end
if isempty(FileDir)
    disp(['requested folder not found']);
end
Variables.FolderName=char(FileDir(Variables.FileNumber))
Variables.FileDir=FileDir; clear FileDir;
%% determine TTL type to analyse (Lickometer / Laser)
if ~isempty(strfind(Variables.FolderName,'Laser'))
Variables.TTLType='Laser'; % Lickometer / Laser
else
Variables.TTLType='Lickometer'; % Lickometer / Laser
end
%% determine the condition
ConditionList={'Jelly1','Jelly2','Jelly 1', 'Jelly 2', 'Empty','Chow', '5ms_2hz', 'Laser_Burst-5ms_2hz','Laser_Single-5ms_2hz','Laser_Single-1ms_2hz','Laser_Single-10ms_2hz','Laser_Single-5ms_1hz','Laser_Single-5ms_5hz','Laser_Single-5ms_10hz','Laser_Single-5ms_20hz','Jelly - Exposure','Jelly - OFF','Jelly','jelly','chow','Empty','Laser'};
for u=1:length(ConditionList)
    if contains(Variables.FolderName,ConditionList{u})
  Variables.Condition=ConditionList{u};
  break
    end
end
%% determine TTL port to analyse
for TTLs=1:1 
if contains(Variables.TTLType,'Lickometer')
Variables.TTLPort=11; % TTL signal for lickometer enters port 1
elseif contains(Variables.TTLType,'Laser')
Variables.TTLPort=10; % TTL signal for lickometer enters port 0
end
end; clear TTLs
Variables.FoodBR=datetime(Variables.Date)<datetime('2021-01-01');

end; clear SelectVars% for
%% Run EphysAnalysis function to get basic data from folder
EphysObj = EphysAnalysis([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\',Variables.FolderName,'\TT',num2str(Variables.TetrodeNumber),'_s.ntt'],Variables);
EphysObj.Variables=Variables;
%% get the cell activity
CellTimeStamp=(EphysObj.timestamps(EphysObj.cells==Variables.UnitNumber)-double(EphysObj.raw_data.TimeStamp(1))); %timestamp in microseconds
%% figure out which TTL to take, if we find laser pulses we take laser TTL, otherwise its Lickometer
% check if selected TTL type exists
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
%% get the TTL activity
TTLTimeStamp=(EphysObj.TTL(w).on-double(EphysObj.raw_data.TimeStamp(1))); %timestamp in microseconds
%% create variable to save
UnitData.FilePath=[EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\',EphysObj.Variables.FolderName];
UnitData.MouseName=EphysObj.Variables.MouseName;
UnitData.Condition=EphysObj.Variables.Condition;
UnitData.UnitTimestampsZeroed=CellTimeStamp; 
UnitData.TTLTimestampsZeroed=TTLTimeStamp;
UnitData.FolderName=EphysObj.Variables.FolderName;  
%% Analyse video based on DLC
[EphysObj] = VideoAnalyse(EphysObj,true,20,4); %(EphysObj,Plot(true/false),Threshold,MinimalBoutTime(seconds))
UnitData.Motifs=EphysObj.Video.Motifs;
UnitData.Variables=EphysObj.Variables;
UnitData.RawData=EphysObj.raw_data;
%% Make scatter plots for each motif
if Variables.PlotScatter
[ArrayForScatter,CollectFrequency]=BoutAnalysis(EphysObj);
end
%% Make scatter plots for each motif
[Obj2Save]=MultipleCompHeatmap(EphysObj);

%% save the file
if EphysObj.Variables.Tagged; IsTagged='Tagged';else; IsTagged='NotTagged';end
if EphysObj.Variables.UseAnova; IsTest='ANOVA';else; IsTest='kruskalwallis';end
if length(EphysObj.Variables.Condition)>4 ; IsCondition='JELLY' ; else IsCondition='CHOW'; end
if EphysObj.Variables.UseAnova; IsTest='ANOVA';else; IsTest='kruskalwallis';end


save(['D:\SummaryMay2024\',EphysObj.Variables.DietType,' ',IsCondition,'\',IsTagged,'\',IsTest,'\',EphysObj.Variables.UnitName,'_Obj.mat'],'Obj2Save');
else
save(['D:\SummaryMay2024\Tagged\kruskalwallis\',EphysObj.Variables.UnitName,'_Obj.mat'],'Obj2Save');
end
else
    if EphysObj.Variables.UseAnova
save(['D:\SummaryMay2024\NotTagged\ANOVA\',EphysObj.Variables.UnitName,'_Obj.mat'],'Obj2Save');
else
save(['D:\SummaryMay2024\NotTagged\kruskalwallis\',EphysObj.Variables.UnitName,'_Obj.mat'],'Obj2Save');
end
end
clearvars -EXCEPT CountGood CountBad Obj2Save BadCount CollectBadFiles EphysObj FileStructure ListOfFiles CollectAllData NameList Count AllEphysObjs EphysObj CollectAllRsquared CollectAllRsquaredNoZeros  DataBase  LaserCollect LickometerCollect




%%    
    end

%     for FindIssues=1:1
%     CollectBadFiles(BadCount).ListName=ListOfFiles(i);
%         CollectBadFiles(BadCount).EphysObjName=EphysObj.Variables.FolderName;
% 
%      CollectBadFiles(BadCount).Diet=EphysObj.Variables.DietType;
%       CollectBadFiles(BadCount).Food=EphysObj.Variables.Condition;
% %check if ttl exists 
% ExtractTTL={EphysObj.TTL.type};
% for w=1:length(ExtractTTL)
%    if ExtractTTL{w}==11
% CollectBadFiles(BadCount).TTLexists=true;
%     break    
%     else
%     CollectBadFiles(BadCount).TTLexists=false;
%     end
% end
% % check if DLC file really exists 
% TestDir=dir([Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\Video']);
% TestDir=extractfield(TestDir,'name');TestDir=TestDir(3:end)';
% DLCcount=1;clear DLCFiles
% for TestFile=1:length(TestDir)
%     if contains(lower(TestDir{TestFile, 1}),'dlc')
%     CollectBadFiles(BadCount).AnyDLC=true;    
%     DLCFiles(DLCcount)={TestDir{TestFile, 1}};
%     DLCcount=DLCcount+1;
%     end
% end
%      if ~exist('DLCFiles')
%   CollectBadFiles(BadCount).AnyDLC=false;
%   CollectBadFiles(BadCount).SpecificDLC=false; 
%      else
%   DLCFiles=DLCFiles';
%      end
% if CollectBadFiles(BadCount).AnyDLC
% for TestFile=1:length(DLCFiles)
% if contains(lower(DLCFiles{TestFile, 1}),lower(EphysObj.Variables.Condition(1:3)))
% CollectBadFiles(BadCount).SpecificDLC=true;
% else
% CollectBadFiles(BadCount).SpecificDLC=false;    
% end    
% end
% end
% BadCount=BadCount+1;
% clearvars -EXCEPT BadCount CollectBadFiles  ListOfFiles Count i  
%     end
continue
    
    catch
    CountBad=CountBad+1;
    continue
    end
end
% clear all

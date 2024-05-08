close all;clear all;clc
%% Run a single file
% FileFolder='C:\Users\netas\Documents\Lammel lab\Obesity paper\EXPERIMENTS\OF\NTS FLOX HFD\Cropped';
FileFolder='C:\Users\netas\Documents\Lammel lab\Obesity paper\EXPERIMENTS\OF\TO analyse';
FileName='BKNTSOE6_1_croppedDLC_resnet50_OpenFieldJun29shuffle1_100000.csv' ; 
% FileName='BKNTSOE10_1_croppedDLC_resnet50_OpenFieldJun29shuffle1_100000.csv' ; 
DisplayPlot=true;
Threshold=20;
[UnitData] = OFVideoAnalyse(FileFolder,FileName,DisplayPlot,Threshold) ;

%% Run batch files
FileFolder='C:\Users\netas\Documents\Lammel lab\Obesity paper\EXPERIMENTS\OF\TO analyse';
% FileFolder='C:\Users\netas\Documents\Lammel lab\Obesity paper\EXPERIMENTS\OF\NTS FLOX HFD\Cropped';
DLCFileDir=dir(FileFolder);
DLCFileDir=extractfield(DLCFileDir,'name');DLCFileDir=DLCFileDir(3:end)';
CollectCSV=[];
for FileNumber=1:length(DLCFileDir)
    if ~isempty(strfind(DLCFileDir{FileNumber, 1},'csv'))
    CollectCSV=[CollectCSV,FileNumber];
    end
end
DLCFiles=DLCFileDir(CollectCSV);
for FileNumber=1:length(DLCFiles)
    try
DisplayPlot=false;
Threshold=20;
FileName=DLCFiles{FileNumber, 1} ; 
[UnitData] = OFVideoAnalyse(FileFolder,FileName,DisplayPlot,Threshold) ;
    catch
    disp([FileName,' not analysed'])
        continue
    end
end

%% collect the data from all files
FileDir=dir('C:\Users\netas\Documents\Lammel lab\Obesity paper\EXPERIMENTS\OF\Matlab files');
FileDir=extractfield(FileDir,'name');FileDir=FileDir(3:end)';
for FileNumber=1:length(FileDir)
DataBase(FileNumber).MatFileName=FileDir{FileNumber, 1};  
load(['C:\Users\netas\Documents\Lammel lab\Obesity paper\EXPERIMENTS\OF\Matlab files\',FileDir{FileNumber, 1}])
DataBase(FileNumber).RawData=UnitData.Motifes;
DataBase(FileNumber).NumberOfBouts=extractfield(UnitData.Motifes,'NumberOfBouts');
DataBase(FileNumber).LengthSeconds=extractfield(UnitData.Motifes,'SumOfLengthSeconds');
DataBase(FileNumber).TotalFrameTime=extractfield(UnitData.Motifes,'TotalFrameTime');
DataBase(FileNumber).TotalFrameTimePrecent=100*extractfield(UnitData.Motifes,'TotalFrameTime')/sum(extractfield(UnitData.Motifes,'TotalFrameTime'));
end
SummaryData.TotalFrameTime=nan(length(DataBase),length(UnitData.Motifes)-1);
SummaryData.NumberOfBouts=nan(length(DataBase),length(UnitData.Motifes)-1);
SummaryData.LengthSeconds=nan(length(DataBase),length(UnitData.Motifes)-1);
for UnitNumber=1:length(DataBase)
SummaryData.TotalFrameTime(UnitNumber,1:length(DataBase(UnitNumber).TotalFrameTime))=DataBase(UnitNumber).TotalFrameTime;
SummaryData.NumberOfBouts(UnitNumber,:)=DataBase(UnitNumber).NumberOfBouts;
SummaryData.LengthSeconds(UnitNumber,:)=DataBase(UnitNumber).LengthSeconds;
SummaryData.TotalFrameTime(UnitNumber,:)=DataBase(UnitNumber).TotalFrameTime;
SummaryData.TotalFrameTimePrecent(UnitNumber,:)=DataBase(UnitNumber).TotalFrameTimePrecent;
end
SummaryData.ListOfConditions={'Rearing';'Stopping';'Walking';'Speeding';'Running';'RightTurn';'LeftTurn'}';








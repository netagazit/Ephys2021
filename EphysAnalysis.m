classdef EphysAnalysis<handle
    properties
        % Data properties of EphysObj:
        Variables       % Stores the variables
        raw_data        % Stores imported InVivoEphys (struct)
        Events          % Stores event data 
        ContData        % Stores continuous data in the 4 channels assosiated with the tetrode file
        data            % Data after manipulations (c,d,n)
        attributes      % Attributes of the sweeps in data (e.g. 'height')
        timestamps      % Timestamps of the waveforms in data (n)
        cells           % Clustered units(n);
        nr_cells        % The total number of clustered units
        TTL             % Imported TTL arrays (struct)
        Video           % Imports Video data
        settings        % Struct with settings.

    end
    
%% Standard methods
    methods
        function EphysObj = EphysAnalysis(filename, Variables)
            %InVivoEphys Construct an instance of this class
                                         
            % Uses Fieldtrip to import a NTT file
            EphysObj.raw_data = read_neuralynx_ntt(filename, 1, inf);
            % Collect the raw
            EphysObj.data = EphysObj.raw_data.dat;
            % Collect clustered cells
            EphysObj.cells = EphysObj.raw_data.CellNumber;
            % Find the total number of clusters
            EphysObj.nr_cells = max(EphysObj.cells);
            % Fill out the timestamps
            EphysObj.timestamps = EphysObj.raw_data.TimeStamp;
            % Find the recording interval
            EphysObj.settings.recording_interval =...
                [EphysObj.raw_data.hdr.FirstTimeStamp, EphysObj.raw_data.hdr.LastTimeStamp];
            % Fill out the settings
            EphysObj.settings.interval = [1, 32];
            EphysObj.settings.preferred_electrodes = [1 2 3]; % three dimensions
            EphysObj.settings.colorlist = ['krmygcgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbmgrbmycgrbm'];
            EphysObj.settings.remove_check = true; % Prevents removal of sweeps without user confirmation
try            
            % Find the DBitVolts (from the header char...)
            header = EphysObj.raw_data.hdr.Header;
            for i = 1:5000
                if strcmp('ADBitVolts', header(i:i+9))
                    break;
                end
            end

            EphysObj.settings.ADBitVolts(1) = str2num(header(i+11:i+36));
            EphysObj.settings.ADBitVolts(2) = str2num(header(i+38:i+63));
            EphysObj.settings.ADBitVolts(3) = str2num(header(i+65:i+90));
            EphysObj.settings.ADBitVolts(4) = str2num(header(i+92:i+117));
catch
[ntt] = read_neuralynx_ntt([filename(1:end-6),'.ntt']);
EphysObj.raw_data.hdr.Heade=ntt.hdr.Header;
header=EphysObj.raw_data.hdr.Heade;  
 % Find the DBitVolts (from the header char...)
            for i = 1:5000
                if strcmp('ADBitVolts', header(i:i+9))
                    break;
                end
            end

            EphysObj.settings.ADBitVolts(1) = str2num(header(i+11:i+36));
            EphysObj.settings.ADBitVolts(2) = str2num(header(i+38:i+63));
            EphysObj.settings.ADBitVolts(3) = str2num(header(i+65:i+90));
            EphysObj.settings.ADBitVolts(4) = str2num(header(i+92:i+117));
 end
            % Collect the recording date
            for i = 1:5000
                if strcmp('(m/d/y)', header(i:i+6))
                    break;
                end
            end
            EphysObj.raw_data.date = header(i+9:i+18);
            
            % Collect the recording time
            for i = 1:5000
                if strcmp('(h:m:s.ms)', header(i:i+9))
                    break;
                end
            end
            EphysObj.raw_data.time = header(i+11:i+22);
            
            % Figure out which electrodes are not recording anything at all
            % and label them 'false' in the working electrodes property
            for i=1:4
                if sum(sum(EphysObj.data(i,:,:)))==0
                    EphysObj.settings.working_electrodes(i) = false;
                else
                    EphysObj.settings.working_electrodes(i) = true;
                end
            end
            
            % Find the Input range(from the header char...)
            for i = 1:5000
                if strcmp('InputRange', header(i:i+9)); break; end
            end
            new_char = header(i+11:i+30);
            for i=1:length(header)
                if strcmp(new_char(i),'-'); break; end
            end
            EphysObj.settings.InputRange = str2num(new_char(1:i-1));
            
            % Find the Threshold Value(from the header char...)
            for i = 1:5000
                if strcmp('ThreshVal', header(i:i+8)); break; end
            end
            new_char = header(i+10:i+30);
            for i=1:length(header)
                if strcmp(new_char(i),'-'); break; end
            end
            EphysObj.settings.ThresVal = str2num(new_char(1:i-1));
           
            % Fill out the attributes
            EphysObj.find_attributes;
            
            % Set the cutoff for display performance issues
            EphysObj.settings.disp_cutoff = 40000;
            
            % Figure out the prefered electrodes (which the most variation in height)
            for i=1:4
                devs(i) = std(EphysObj.attributes(i).height);
            end
            [~, show_electrodes] = sort(devs, 'descend');
            show_electrodes = sort(show_electrodes(1:3));
            EphysObj.settings.preferred_electrodes = show_electrodes;
 %% Get Event (TTLs) file 
            % Load any associated event file
            path = fileparts(filename);
            if ~isempty(path)
                event_file = [path '\Events.nev'];
            else
                event_file = 'Events.nev';
            end
            if isfile(event_file)
Events = read_neuralynx_nev(event_file);
%% find event port for each TTL
for i=1:length(Events)
% go through the "ON" Events and figure out the port of each TTL
% event port 10 is port 0 11 is port 1 atc.
if contains(Events(i).EventString,' value (0x0000)')
Events(i).EventPort=nan;
elseif contains(Events(i).EventString,'port 0')
Events(i).EventPort=10;
elseif contains(Events(i).EventString,'port 1')
Events(i).EventPort=11;
elseif contains(Events(i).EventString,'port 2')
Events(i).EventPort=12;
elseif contains(Events(i).EventString,'port 3')
Events(i).EventPort=13;
else
Events(i).EventPort=nan;
end %if
end; clear i % end for
%% mark the TTLs not from the selected port as false
for i=flip(1:length(Events))
if Events(i).EventPort==Variables.TTLPort
Events(i).SelectedPort=true;
else
Events(i).SelectedPort=false;
end
end; clear i
EphysObj.Events=Events;
          % Load the file
                disp ('Event file found... importing TTL.')
                events =Events; %
                %%
                for i=1:length(events)
                    TimestampsEvents(i) = events(i).TimeStamp;
                     ev_indexer(i) = events(i).EventPort; %  Endexer(i) = events(i).TTLValue;
                end
                % Make array with all possible TTL types
                TTL_types = unique(ev_indexer);
                TTL_types(isnan(TTL_types))=[];    
                TTL_types = TTL_types(TTL_types~=0);
 if length(TTL_types)==1 
                for i=1:length(TTL_types)
                    EphysObj.TTL(i).type = TTL_types(i);
                    indexer = (ev_indexer == TTL_types(i));
                    EphysObj.TTL(i).on = TimestampsEvents(indexer);
                    EphysObj.TTL(i).off = TimestampsEvents([false, indexer(1:end-1)]);
                    EphysObj.TTL(i).ave_duration = mean(EphysObj.TTL(i).off - EphysObj.TTL(i).on);
                end
 else
                for i=1:length(TTL_types) %length(TTL_types)-1
                    EphysObj.TTL(i).type = TTL_types(i);
                    indexer = (ev_indexer == TTL_types(i));
                    EphysObj.TTL(i).on = TimestampsEvents(indexer);
                    EphysObj.TTL(i).off = TimestampsEvents([false, indexer(1:end-1)]);
                    EphysObj.TTL(i).ave_duration = mean(EphysObj.TTL(i).off - EphysObj.TTL(i).on);
                end
end 
            end
  %% now get continuous data
try
EphysObj.ContData.Ch1=read_neuralynx_ncs([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\',Variables.FolderName,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+1),'.ncs']) ;EphysObj.ContData.Ch1.ChannelNumber=4*(Variables.TetrodeNumber-1)+1;
EphysObj.ContData.Ch2=read_neuralynx_ncs([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\',Variables.FolderName,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+2),'.ncs']) ;EphysObj.ContData.Ch2.ChannelNumber=4*(Variables.TetrodeNumber-1)+2;
EphysObj.ContData.Ch3=read_neuralynx_ncs([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\',Variables.FolderName,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+3),'.ncs']) ;EphysObj.ContData.Ch3.ChannelNumber=4*(Variables.TetrodeNumber-1)+3;
EphysObj.ContData.Ch4=read_neuralynx_ncs([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\',Variables.FolderName,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+4),'.ncs']) ;EphysObj.ContData.Ch4.ChannelNumber=4*(Variables.TetrodeNumber-1)+4;
catch
[EphysObj.ContData.Ch1.TimeStamp, ~, ~, ~, ~, ~] = Nlx2MatCSC([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\',Variables.FolderName,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+1),'.ncs'], [1 1 1 1 1], 1, 1, []); EphysObj.ContData.Ch1.ChannelNumber=4*(Variables.TetrodeNumber-1)+1;
[EphysObj.ContData.Ch2.TimeStamp, ~, ~, ~, ~, ~] = Nlx2MatCSC([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\',Variables.FolderName,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+2),'.ncs'], [1 1 1 1 1], 1, 1, []); EphysObj.ContData.Ch2.ChannelNumber=4*(Variables.TetrodeNumber-1)+2;
[EphysObj.ContData.Ch3.TimeStamp, ~, ~, ~, ~, ~] = Nlx2MatCSC([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\',Variables.FolderName,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+3),'.ncs'], [1 1 1 1 1], 1, 1, []); EphysObj.ContData.Ch3.ChannelNumber=4*(Variables.TetrodeNumber-1)+3;
[EphysObj.ContData.Ch4.TimeStamp, ~, ~, ~, ~, ~] = Nlx2MatCSC([Variables.ComputerDir,'\',Variables.Date,'\',Variables.MouseName,'\',Variables.FolderName,'\CSC',num2str(4*(Variables.TetrodeNumber-1)+4),'.ncs'], [1 1 1 1 1], 1, 1, []); EphysObj.ContData.Ch4.ChannelNumber=4*(Variables.TetrodeNumber-1)+4;
end
EphysObj.settings.FirstTimestamp=EphysObj.ContData.Ch3.TimeStamp(1);
EphysObj.Variables.TimeLimitStamp=EphysObj.settings.FirstTimestamp+Variables.TimeLimit*60*1000000;
for i=1:length(TTL_types)
EphysObj.TTL(i).on=EphysObj.TTL(i).on(EphysObj.TTL(i).on<EphysObj.Variables.TimeLimitStamp);
EphysObj.TTL(i).off=EphysObj.TTL(i).off(EphysObj.TTL(i).off<EphysObj.Variables.TimeLimitStamp);
end

        end
 %%      
        function find_attributes(EphysObj)
            % Find attributes updates all the sweep attributes
            
            % Analyse the following interval
            start = EphysObj.settings.interval(1);
            stop = EphysObj.settings.interval(2);
            
            % get the data from the 4 electrodes
            for i=1:4
                temp_data = squeeze(EphysObj.data(i,start:stop,:)).*EphysObj.settings.ADBitVolts(i)*10^6;
                
                % get the peak
                EphysObj.attributes(i).peak = max(temp_data);
                
                % get the valley
                EphysObj.attributes(i).valley = min(temp_data);
                
                % get the height
                EphysObj.attributes(i).height = EphysObj.attributes(i).peak - EphysObj.attributes(i).valley;
              
                  end
            
        end
           
 
        function [EphysObj] = VideoAnalyse(EphysObj,DisplayPlot,Threshold,MinimalBoutTime) 
Events = EphysObj.Events; 
% EphysObj.Variables.FoodBR=datetime(EphysObj.Variables.Date)<datetime('2021-01-01');
%% find video file
VideoDirectory=[EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\Video'];
VideoDirectory=dir(VideoDirectory);
RelevantVideoFiles=[];
Condition=EphysObj.Variables.Condition(find(~isspace(EphysObj.Variables.Condition)));
  Condition=strrep(Condition, '-', ' ');
  Condition=strrep(Condition, '_', ' ');
  Condition=lower(Condition(find(~isspace(Condition))));
for FindRelevantFiles=3:length(VideoDirectory)
   VideoDirectoryTemp=VideoDirectory(FindRelevantFiles,1).name;
   % cange to lowercase and take out '-', '_' and spaces
   VideoDirectoryTemp=strrep(VideoDirectoryTemp, '-', ' ');
   VideoDirectoryTemp=strrep(VideoDirectoryTemp, '_', ' ');
   VideoDirectoryTemp=lower(VideoDirectoryTemp(find(~isspace(VideoDirectoryTemp))));
     if ~isempty(strfind(VideoDirectoryTemp,Condition))
    if ~isempty(strfind(VideoDirectoryTemp,'avi'))
        if ~isempty(strfind(VideoDirectoryTemp,'cropped'))
        EphysObj.Video.VideoAVI=VideoDirectory(FindRelevantFiles,1).name;
        else
        EphysObj.Video.VideoAVIOriginal=VideoDirectory(FindRelevantFiles,1).name;
        end
    elseif ~isempty(strfind(VideoDirectoryTemp,'csv'))
        if ~isempty(strfind(VideoDirectoryTemp,'dlc'))
           EphysObj.Video.VideoDLC=VideoDirectory(FindRelevantFiles,1).name;         
        else
        EphysObj.Video.VideoCsv=VideoDirectory(FindRelevantFiles,1).name;
        end
    elseif ~isempty(strfind(VideoDirectoryTemp,'mp4'))
        EphysObj.Video.VideoMp4=VideoDirectory(FindRelevantFiles,1).name;
    end
else
    continue
end

end
% look at the track from bioobserve
try
EphysObj.Video.Track = biobserve_import([EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\Video\',EphysObj.Video.VideoCsv]);
catch
end
% other variables
timestamp_interval = 5; %Sec - this is the interval of the TTL signal from the bioobserve to NLX
try
Info=readtable([EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\Video\',EphysObj.Video.VideoDLC]);
catch
Info=readtable([EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\Video\',EphysObj.Video.VideoCsv]);
end 
InfoArray=table2array(Info);
% read the video file
VideoInfo=VideoReader([EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\Video\',EphysObj.Video.VideoAVI]);EphysObj.Video.VideoInfo=VideoInfo;
FrameTimes=linspace(0,EphysObj.Video.VideoInfo.Duration,EphysObj.Video.VideoInfo.NumFrames); % X values for frame basd analysis
EphysObj.Video.FrameEdgesSec=linspace(0,FrameTimes(end)*2-mean([FrameTimes(end),FrameTimes(end-1)]),EphysObj.Video.VideoInfo.NumFrames+1);
EphysObj.Video.FrameEdgesMin=EphysObj.Video.FrameEdgesSec/60;
%% work on the DeepLabCut output
% Take out data that has low accuracy
Temp=InfoArray(:,2:3);Temp(InfoArray(:,4)<0.85,1:2)=nan;Struct.Nose=Temp;clear Temp 
Temp=InfoArray(:,5:6);Temp(InfoArray(:,7)<0.85,1:2)=nan;Struct.Implant=Temp;clear Temp
Temp=InfoArray(:,8:9);Temp(InfoArray(:,10)<0.85,1:2)=nan;Struct.Body=Temp;clear Temp
Temp=InfoArray(:,11:12);Temp(InfoArray(:,13)<0.85,1:2)=nan;Struct.Tail=Temp;clear Temp
Temp=InfoArray(:,14:15);Temp(InfoArray(:,16)<0.85,1:2)=nan;Struct.TopLeft=Temp;clear Temp
Temp=InfoArray(:,17:18);Temp(InfoArray(:,19)<0.85,1:2)=nan;Struct.BottomRight=Temp;clear Temp
EphysObj.Video.DLC_Struct=Struct; % save to EphysObj
% collect data by frames
NewArray=nan(length(InfoArray),12);
NewArray(:,1:2)=Struct.Nose;
NewArray(:,3:4)=Struct.Implant;
NewArray(:,5:6)=Struct.Body;
NewArray(:,7:8)=Struct.Tail;
NewArray(:,9:10)=Struct.TopLeft;
NewArray(:,11:12)=Struct.BottomRight;
% define the arena
MaxX=max([NewArray(:,[1 3 5 7 9 11])]);MaxX=max(MaxX);
MinX=min([NewArray(:,[1 3 5 7 9 11])]);MinX=min(MinX);
MaxY=max(NewArray(:,[2 4 6 8 10 12]));MaxY=max(MaxY);
MinY=min(NewArray(:,[2 4 6 8 10 12]));MinY=min(MinY);
ArenaLength=MaxX-MinX;
ArenaLengthFactor=33/ArenaLength;
ArenaWidth=MaxY-MinY;
ArenaWidthFactor=20/ArenaWidth;
ArenaFactor=mean([ArenaLengthFactor,ArenaWidthFactor]);
ArenaAll(1,1:2)=[MinX MaxY];Arena.BottomLeft=[MinX MaxY];
ArenaAll(3,1:2)=[MinX MinY];Arena.TopLeft=[MinX MinY];
ArenaAll(2,1:2)=[MaxX MaxY];Arena.BottomRight=[MaxX MaxY];
ArenaAll(4,1:2)=[MaxX MinY];Arena.TopRight=[MaxX MinY];
EphysObj.Video.ArenaAll=ArenaAll;% save to EphysObj
% define BottomRight location
CenterBottomRightX=median(Struct.BottomRight(:,1),'omitnan');
CenterBottomRightY=median(Struct.BottomRight(:,2),'omitnan');
EphysObj.Video.PlateBRRelative=[CenterBottomRightX/MaxX,CenterBottomRightY/MaxY];
CenterTopLeftX=median(Struct.TopLeft(:,1),'omitnan'); %%CenterTopLeftX=MaxX-nanmedian(Struct.BottomRight(:,1));
CenterTopLeftY=median(Struct.TopLeft(:,2),'omitnan'); %%CenterTopLeftY=MaxY-nanmedian(Struct.BottomRight(:,2));
EphysObj.Video.PlateTLRelative=[CenterTopLeftX/MaxX,CenterTopLeftY/MaxY];
EphysObj.Video.DistanceBetweenPlates=sqrt((CenterTopLeftX-CenterBottomRightX).^2 + (CenterTopLeftY-CenterBottomRightY).^2);
EphysObj.Video.DistanceBetweenPlatesXY=[CenterTopLeftX-CenterBottomRightX,CenterTopLeftY-CenterBottomRightY];
% correct to the case one of the plates is detected in the wrong location
for Plot=1:1
    if DisplayPlot
      figure;scatter(Arena.TopLeft(1),-1*Arena.TopLeft(2)); hold on
      scatter(Arena.TopRight(1),-1*Arena.TopRight(2)); hold on
       scatter(Arena.BottomLeft(1),-1*Arena.BottomLeft(2)); hold on
      scatter(Arena.BottomRight(1),-1*Arena.BottomRight(2)); hold on
      scatter(CenterBottomRightX,-1*CenterBottomRightY); hold on
      scatter(CenterTopLeftX,-1*CenterTopLeftY); hold on
legend('ArenaTL','ArenaTR','ArenaBL','ArenaBR','PlateBR', 'PlateTL')
      
    end
end
for FindWrongPlate=1:1
if EphysObj.Video.DistanceBetweenPlates<100
% check which plate is in the right place
% since the plates are too close, we need to find which one is correct. the one closer to the relevant corner
DistanceTLplate2TLCorner= sqrt((Arena.TopLeft(1)-CenterTopLeftX).^2 + (Arena.TopLeft(2)-CenterTopLeftY).^2);
DistanceBRplate2TLCorner= sqrt((Arena.TopLeft(1)-CenterBottomRightX).^2 + (Arena.TopLeft(2)-CenterBottomRightY).^2);
DistanceTLplate2BRCorner= sqrt((Arena.BottomRight(1)-CenterTopLeftX).^2 + (Arena.BottomRight(2)-CenterTopLeftY).^2);
DistanceBRplate2BRCorner= sqrt((Arena.BottomRight(1)-CenterBottomRightX).^2 + (Arena.BottomRight(2)-CenterBottomRightY).^2);
% deal with a wrong plate assingment
if mean([DistanceBRplate2BRCorner,DistanceTLplate2BRCorner])-mean([DistanceTLplate2TLCorner,DistanceBRplate2TLCorner])>0
ReferanceCorner=Arena.BottomRight; %disp('closer2TL corner')
% what plate is closer
if DistanceTLplate2TLCorner-DistanceBRplate2TLCorner>0
    Plate2Modify=[CenterTopLeftX,CenterTopLeftY]; %    disp('BR plate closer to TL corner')
           BottomRightPlate=false;
else
       Plate2Modify=[CenterBottomRightX,CenterBottomRightY];%   disp('TL plate closer to TL corner')
       BottomRightPlate=true;
end
else
ReferanceCorner=Arena.TopLeft;%disp('closer2BR corner')
if DistanceTLplate2BRCorner-DistanceBRplate2BRCorner>0
       Plate2Modify=[CenterTopLeftX,CenterTopLeftY];%    disp('BR plate closer to BR corner')
       BottomRightPlate=false;
else
       Plate2Modify=[CenterBottomRightX,CenterBottomRightY];%    disp('TL plate closer to BR corner')
       BottomRightPlate=true;
end
end % if
if BottomRightPlate % need to modify coordinates of bottom right
CenterBottomRightX=0.9*ReferanceCorner(1);
CenterBottomRightY=0.8*ReferanceCorner(2);
if DisplayPlot
    figure
scatter(CenterBottomRightX,-1*CenterBottomRightY); end
legend('ArenaTL','ArenaTR','ArenaBL','ArenaBR','PlateBR', 'PlateTL','ModifiedBRPlateLocation')
else % need to modify coordinates of top left
CenterTopLeftX=0.9*ReferanceCorner(1);
CenterTopLeftY=0.8*ReferanceCorner(2);
if DisplayPlot
scatter(CenterTopLeftX,-1*CenterTopLeftY); end
end
end % if
end % for

%% calculate all the distances between the markers. use the CenterBottomRight/TopLeft
NoseArray=nan(size(NewArray,1),size(NewArray,2)/2);
NewArrayX=NewArray(:,[1:2:11]); % get x values of array
NewArrayY=NewArray(:,[2:2:12]); % get Y values of array
EphysObj.Video.DLC_Raw=NewArray; %save to EphysObj
%% calculate distances between points
fields = fieldnames(Struct);
for FieldNumber=1:length(fields)
DistancesArray(FieldNumber).FieldName=fields(FieldNumber);
DistancesArray(FieldNumber).Distances=sqrt((NewArrayX-NewArrayX(:,FieldNumber)).^2 + (NewArrayY-NewArrayY(:,FieldNumber)).^2);
DistancesArray(FieldNumber).Distances=DistancesArray(FieldNumber).Distances(:,1:4);
end
DistancesArray(7).FieldName='CenterBottomRight';
DistancesArray(7).Distances=sqrt((NewArrayX-CenterBottomRightX).^2 + (NewArrayY-CenterBottomRightY).^2);
DistancesArray(8).Distances=sqrt((NewArrayX-CenterTopLeftX).^2 + (NewArrayY-CenterTopLeftY).^2);
DistancesArray(8).FieldName='CenterTopLeft';
EphysObj.Video.DistanceArray=DistancesArray;
%% compare FR and video
% Align the times
VideoTTL=[EphysObj.TTL.type];LocationVideoTTL=find(VideoTTL==13); %find the releant TTL row
EphysObj.Video.VideoTTLTimestamps=EphysObj.TTL(LocationVideoTTL).on; %save to EphysObj
EphysObj.Video.VideoStartTimestamp=EphysObj.Video.VideoTTLTimestamps(1)-5*1000000;% get the video start timestamp (5 sec before the first video ttl)
% find the NLX first and last timestamps, save to EphysObj
TTLString={EphysObj.Events.EventString}; % collect events strings
RecordingStart=find(ismember(TTLString,'Starting Recording'));EphysObj.Video.RecordingFirstTimeStamp=EphysObj.Events(RecordingStart).TimeStamp;
RecordingEnd=find(ismember(TTLString,'Stopping Recording'));EphysObj.Video.RecordingLastTimeStamp=EphysObj.Events(RecordingEnd).TimeStamp;
% determine what started first, the video or the recording
if EphysObj.Video.RecordingFirstTimeStamp>EphysObj.Video.VideoTTLTimestamps;RecordedFirst='Video'; else;RecordedFirst='NLX';end
UnitData=EphysObj.timestamps(EphysObj.cells==EphysObj.Variables.UnitNumber); % collect unit timestamps
% get easy to understand information about the experiment:
for SimpleNumbers=1:1
EphysObj.Video.SimpleNumbers.RecordingStart=(EphysObj.Video.RecordingFirstTimeStamp-EphysObj.Video.RecordingFirstTimeStamp)/1000000;
EphysObj.Video.SimpleNumbers.RecordingStopped=(EphysObj.Video.RecordingLastTimeStamp-EphysObj.Video.RecordingFirstTimeStamp)/1000000;
EphysObj.Video.SimpleNumbers.VideoStartTimestamp=(EphysObj.Video.VideoStartTimestamp-EphysObj.Video.RecordingFirstTimeStamp)/1000000;
EphysObj.Video.SimpleNumbers.VideoLastTTLTimestamp=(EphysObj.Video.VideoTTLTimestamps(end)-EphysObj.Video.RecordingFirstTimeStamp)/1000000;
EphysObj.Video.SimpleNumbers.VideoExpectedLastTimestamp=(EphysObj.Video.VideoStartTimestamp+EphysObj.Video.VideoInfo.Duration*1000000-EphysObj.Video.RecordingFirstTimeStamp)/1000000;
EphysObj.Video.SimpleNumbers.FirstUnitTimestamp=(UnitData(1)-EphysObj.Video.RecordingFirstTimeStamp)/1000000;
EphysObj.Video.SimpleNumbers.LastUnitTimestamp=(UnitData(end)-EphysObj.Video.RecordingFirstTimeStamp)/1000000;
end
% take out the parts of the NLX recording that didn't have video recording 
VideoOverlapData=UnitData(UnitData>EphysObj.Video.VideoStartTimestamp);
EphysObj.Video.VideoOverlapData=VideoOverlapData(VideoOverlapData<(EphysObj.Video.VideoStartTimestamp+EphysObj.Video.VideoInfo.Duration*1000000));
EphysObj.Video.NormData=EphysObj.Video.VideoOverlapData-EphysObj.Video.VideoStartTimestamp;% normalize unit data to video timeframe
% get the sensor data to plot
BinsForTimeHistogram=1000000*[0:1:(EphysObj.Video.VideoInfo.Duration)] ; % bin the data from 0 to the recording length in seconds
BinsForTimeHistogramFrame=1000000*[0:1/EphysObj.Video.VideoInfo.FrameRate:(EphysObj.Video.VideoInfo.Duration)] ; % bin the data from 0 to the recording length in seconds

if ~isempty(find([EphysObj.TTL.type]==11))
EphysObj.Video.NormSensorEvents=double(EphysObj.TTL(find([EphysObj.TTL.type]==11)).on)-double(EphysObj.Video.VideoStartTimestamp); % normalize sensor data to video timeframe
[SensorN,~,~] = histcounts(EphysObj.Video.NormSensorEvents,BinsForTimeHistogram);EphysObj.Video.SensorN=SensorN;% get frequency distrebution for sensor
[SensorNPerFrame,~,~] = histcounts(EphysObj.Video.NormSensorEvents,BinsForTimeHistogramFrame);EphysObj.Video.SensorN=SensorN;% get frequency distrebution for sensor
EphysObj.Video.SensorNPerFrame=SensorNPerFrame;
else
    disp('Requested TTL type not found, no sensor hits')
end
[DataN,EDGES,~] = histcounts(EphysObj.Video.NormData,BinsForTimeHistogram); EphysObj.Video.DataN=DataN;EphysObj.Video.EDGES=EDGES;% get frequency distrebution for data
[DataNPerFrame,EDGESPerFrame,~] = histcounts(EphysObj.Video.NormData,BinsForTimeHistogramFrame); EphysObj.Video.DataN=DataN;EphysObj.Video.EDGES=EDGES;% get frequency distrebution for data
EphysObj.Video.DataNPerFrame=DataNPerFrame;
Xvalues=linspace(0,EDGES(end),length(EDGES)-1);% get x centers to plot the FR and Sensor rate
XvaluesPerFrame=linspace(0,EDGESPerFrame(end),length(EDGESPerFrame)-1);% get x centers to plot the FR and Sensor rate
EphysObj.Video.XvaluesPerFrame=XvaluesPerFrame;
for Plot=1:1
    if DisplayPlot
%% Plot the general track
figure 
plot(Struct.Nose(:,1)',-1*Struct.Nose(:,2)');hold on
plot(Struct.Implant(:,1),-1*Struct.Implant(:,2));hold on
plot(Struct.Body(:,1),-1*Struct.Body(:,2));hold on
plot(Struct.Tail(:,1),-1*Struct.Tail(:,2));hold on
plot(Struct.TopLeft(:,1),-1*Struct.TopLeft(:,2));hold on
plot(Struct.BottomRight(:,1),-1*Struct.BottomRight(:,2));hold on
% xlim([0 MaxX]); ylim([0 -MaxY]);
legend('Nose','Implant','Body','Tail','TopLeft','BottomRight')
    end % if
end %for plot
figure
%% figure out the implant position relative to the plates
for i=[1:4]
EphysObj.Video.DistanceArray(i).AverageTL=mean(EphysObj.Video.DistanceArray(8).Distances(:,i),'omitnan');
EphysObj.Video.DistanceArray(i).AverageBR=mean(EphysObj.Video.DistanceArray(7).Distances(:,i),'omitnan');
end% test what plate did the mouse spent more time next to 
% collect the spikes for diffrent behaviors
ImplantDistances.names={'Body','Tail','BR','TL'};
ImplantDistances.KeyPoints(:,1)=(EphysObj.Video.DistanceArray(3).Distances(:,2)); % distance of implant from body
ImplantDistances.KeyPoints(:,2)=(EphysObj.Video.DistanceArray(4).Distances(:,2)); % distance of implant from tail
ImplantDistances.KeyPoints(:,3)=(EphysObj.Video.DistanceArray(7).Distances(:,2)); % distance of implant from BR
ImplantDistances.KeyPoints(:,4)=(EphysObj.Video.DistanceArray(8).Distances(:,2)); % distance of implant from TL
Times=FrameTimes/60;
EphysObj.Video.Times=Times; % in minutes
% collect the times where there is threshold crossing
RearingFrames=find([EphysObj.Video.DistanceArray(3).Distances(:,2)]<(Threshold/2)); EphysObj.Video.RearingFrames=RearingFrames;% small distance of body and head indicates rearing
Close2TLFrames=find([EphysObj.Video.DistanceArray(8).Distances(:,2)]<Threshold); % small distance indicates proximity to plate
Close2TLFrames=setdiff(Close2TLFrames,RearingFrames); EphysObj.Video.Close2TLFrames=Close2TLFrames;% if the mouse is rearing it is not eating
Close2BRFrames=find([EphysObj.Video.DistanceArray(7).Distances(:,2)]<Threshold); % small distance indicates proximity to plate
Close2BRFrames=setdiff(Close2BRFrames,RearingFrames);EphysObj.Video.Close2BRFrames=Close2BRFrames;% if the mouse is rearing it is not eating
ExcludeFrames=[RearingFrames;Close2TLFrames;Close2BRFrames];
%% calculate angele between points
for GetAngles=1:1
XYBody = [NewArray(:,5) NewArray(:,6)]; % body
XYHead = [NewArray(:,3) NewArray(:,4)]; % implant
XYTail = [NewArray(:,7) NewArray(:,8)]; % tail origin
VectorHB=(XYHead-XYBody); VectorBT=(XYBody-XYTail);
for Temp=1:length(VectorHB)
BodyHeadAngle(Temp)= (180/pi)*atan2(VectorHB(Temp,1)*VectorBT(Temp,2)-VectorBT(Temp,1)*VectorHB(Temp,2),dot(VectorHB(Temp,:),VectorBT(Temp,:)));
end
BodyHeadAngle=BodyHeadAngle';EphysObj.Video.Angle=BodyHeadAngle';

for FrameNumber=1:length(NewArray)-1
H0=XYHead(FrameNumber,1:2);
H1=XYHead(FrameNumber+1,1:2);
B0=XYBody(FrameNumber,1:2);  
B1=XYBody(FrameNumber+1,1:2); 
FrameAngle(FrameNumber) = (180/pi)*(atan2(abs(det([H1-B1;H0-B0])),dot(H1-B1,H0-B0))); % correct
V0=[XYHead(FrameNumber,1)-XYBody(FrameNumber,1),XYHead(FrameNumber,2)-XYBody(FrameNumber,2)];
V1=[XYHead(FrameNumber+1,1)-XYBody(FrameNumber+1,1),XYHead(FrameNumber+1,2)-XYBody(FrameNumber+1,2)];
SignedFrameAngle(FrameNumber)=atan2d(V0(1)*V1(2)-V0(2)*V1(1),V0(1)*V1(1)+V0(2)*V0(2));
FrameAngle(FrameNumber) = (180/pi)*(atan2(abs(det([H1-B1;H0-B0])),dot(H1-B1,H0-B0))); % correct
end; FrameAngle=FrameAngle';SignedFrameAngle=SignedFrameAngle';
EphysObj.Video.FrameAngle=FrameAngle;EphysObj.Video.SignedFrameAngle=SignedFrameAngle';
% get a sample of optional angles
Angles2Test=linspace(-180,180,13);
BodyHeadAngleNOExcluded=BodyHeadAngle;
BodyHeadAngleNOExcluded(ExcludeFrames)=nan;
for AngleRange=1:length(Angles2Test)-1
GroupAngle(AngleRange).Frame=(find(Angles2Test(AngleRange)<BodyHeadAngleNOExcluded&BodyHeadAngleNOExcluded<Angles2Test(AngleRange+1)));
GroupAngle(AngleRange).Values=(BodyHeadAngleNOExcluded(GroupAngle(AngleRange).Frame));
end
EphysObj.Video.GroupAngle=GroupAngle;
LeftTurnFrames=GroupAngle(9).Frame;%90
RightTurnFrames=GroupAngle(4).Frame;%-90
end %GetAngles
EphysObj.Video.RightTurnFrames=RightTurnFrames;
EphysObj.Video.LeftTurnFrames=LeftTurnFrames;
% work on the rest of the frames to get velocity
AllOfFrames=[1:1:VideoInfo.NumFrames]';
ExcludeFrames=[RearingFrames;Close2TLFrames;Close2BRFrames;RightTurnFrames;LeftTurnFrames];
%%
for GetVelocity=1:1
% calculate the velocity in the frames to decide if the mouse is walking or standing
BodyLocationBeforeFrame=EphysObj.Video.DLC_Struct.Body;
BodyLocationAfterFrame=[EphysObj.Video.DLC_Struct.Body(2:end,:); EphysObj.Video.DLC_Struct.Body(end,:)];
% calculate the distacne of body position between frames
BodyDistanceBetweenFrames=sqrt((BodyLocationAfterFrame(:,1)-BodyLocationBeforeFrame(:,1)).^2 + (BodyLocationAfterFrame(:,2)-BodyLocationBeforeFrame(:,2)).^2);% distance in pixel
BodyDistanceBetweenFramesNorm=BodyDistanceBetweenFrames*ArenaFactor*15; % cm/sec %15 is the sampling rate
%%
BodyDistanceBetweenFramesExcluded=BodyDistanceBetweenFramesNorm;
BodyDistanceBetweenFramesExcluded(ExcludeFrames')=nan;
StoppingFrames=find(BodyDistanceBetweenFramesExcluded<1); % threshold for low velocity
WalkingFrames=find(BodyDistanceBetweenFramesExcluded>=1&BodyDistanceBetweenFramesExcluded<13);
RunningFrames=find(BodyDistanceBetweenFramesNorm>=13&BodyDistanceBetweenFramesExcluded<70);% threshold for high velocity
RestOfFrames=setdiff(AllOfFrames,[RightTurnFrames;LeftTurnFrames;StoppingFrames;RunningFrames;RearingFrames;Close2TLFrames;Close2BRFrames;WalkingFrames]);
end
EphysObj.Video.RunningFrames=RunningFrames;
EphysObj.Video.StoppingFrames=StoppingFrames;
EphysObj.Video.WalkingFrames=WalkingFrames;
EphysObj.Video.RestOfFrames=RestOfFrames;
%% save all frames to one variable
AllFrames(1)={Close2TLFrames};
AllFrames(2)={Close2BRFrames};
AllFrames(3)={RearingFrames};
AllFrames(4)={WalkingFrames};
AllFrames(5)={RunningFrames};
AllFrames(6)={StoppingFrames};
AllFrames(7)={RightTurnFrames};
AllFrames(8)={LeftTurnFrames};
AllFrames(9)={RestOfFrames};
AllFrames=AllFrames';
%% get the relevant times from the frames
Close2TLTimes=Times(Close2TLFrames);
Close2BRTimes=Times(Close2BRFrames); 
RearingTimes=Times(RearingFrames); 
WalkingTimes=Times(WalkingFrames); 
RunningTimes=Times(RunningFrames);
StoppingTimes=Times(StoppingFrames);
RightTurnTimes=Times(RightTurnFrames);
LeftTurnTimes=Times(LeftTurnFrames);
RestOfFramesTimes=Times(RestOfFrames);
EphysObj.Video.AllFrames=AllFrames;
% for Plot=1:1
%     if DisplayPlot
%         % plot velocity
%         figure;
% yyaxis left
% plot(EphysObj.Video.Times,(EphysObj.Video.BodyDistanceBetweenFramesNorm/60));hold on
% ylim([0 5])
% yyaxis right
% plot(linspace(0,EphysObj.Video.EDGES(end),length(EphysObj.Video.EDGES)-1)/60000000,(EphysObj.Video.DataN),'-k');hold on;  % plot data in black
% ylim([0 60])
% legend('Velocity','Unit Activity')
% title ('Velocity vs unit activity')
xlim([0 15])
% plot the activity on one timeline
BehaviorFigure=figure;
% if BR is food plate
if EphysObj.Variables.FoodBR
scatter(Close2BRTimes,zeros(1,length(Close2BRTimes)),'|','MarkerEdgeColor',[0.91, 0.153, 0.557]);hold on % pink
scatter(Close2TLTimes,zeros(1,length(Close2TLTimes)),'|','MarkerEdgeColor',[0.91 0.827 0.153]);hold on
else
scatter(Close2TLTimes,zeros(1,length(Close2TLTimes)),'|','MarkerEdgeColor',[0.91, 0.153, 0.557]);hold on
scatter(Close2BRTimes,zeros(1,length(Close2BRTimes)),'|','MarkerEdgeColor',[0.91 0.827 0.153]);hold on % pink
end
scatter(RearingTimes,zeros(1,length(RearingTimes)),'|', 'MarkerEdgeColor',[0.447 0.569 0.475]);hold on
scatter(WalkingTimes,zeros(1,length(WalkingTimes)),'|','MarkerEdgeColor',[0.98 0.592 0]);hold on
scatter(RunningTimes,zeros(1,length(RunningTimes)),'|','MarkerEdgeColor',[0.369 0.408 0.8]);hold on
scatter(StoppingTimes,zeros(1,length(StoppingTimes)),'|','MarkerEdgeColor',[0.231 0.741 0.871]);hold on
scatter(RightTurnTimes,zeros(1,length(RightTurnTimes)),'|','MarkerEdgeColor',[0.094 0.871 0.224]);hold on
scatter(LeftTurnTimes,zeros(1,length(LeftTurnTimes)),'|','MarkerEdgeColor',[0.859 0.631 0.173]);hold on
% scatter(RestOfFramesTimes,ones(1,length(RestOfFramesTimes)),'|','k');hold on
if EphysObj.Variables.PlotZscore
try plot(Xvalues/60000000,zscore(DataN,1,'all'),'-k');hold on; catch;end % plot data in black
try plot(Xvalues/60000000,(-1*zscore(SensorN,1,'all')),'-r');hold on; catch;end % plot sensor in red

else
try plot(Xvalues/60000000,(DataN),'-k');hold on; catch;end % plot data in black
try plot(Xvalues/60000000,(-1*SensorN),'-r');hold on; catch;end % plot sensor in red
end


legend('Food Plate','Empty Plate','Rearing','Walking','Running','Stopping','RightTurn','LeftTurn','Sensor','Unit')
xlabel('Time (minutes)'); ylabel('Sensor frequency (hz)');
% ylim([-60 60])
xlim([0 15])
% %% Plot a selected frame with markers
% count=1; figure
% for MotifNumber=1:length(AllFrames)-1
%     TempFrames=AllFrames{MotifNumber, 1} ; 
%     SelectedFrames(MotifNumber,1:5)=[randsample(TempFrames',5)];
%  for Frame2Plot=SelectedFrames(MotifNumber,1:5) 
%      subplot(length(AllFrames)-1,5,count)
% scatter(Struct.Nose(Frame2Plot,1),-1*Struct.Nose(Frame2Plot,2));hold on
% scatter(Struct.Implant(Frame2Plot,1),-1*Struct.Implant(Frame2Plot,2));hold on
% scatter(Struct.Body(Frame2Plot,1),-1*Struct.Body(Frame2Plot,2));hold on
% scatter(Struct.Tail(Frame2Plot,1),-1*Struct.Tail(Frame2Plot,2));hold on
% scatter(CenterBottomRightX,-1*CenterBottomRightY);hold on
% scatter(CenterTopLeftX,-1*CenterTopLeftY);hold on
% xlim([0 MaxX]); ylim([-1*MaxY 0]);
% title(num2str(Times(Frame2Plot)))%subtitle(num2str(FrameTimes(Frame)))
% count=count+1;
%  end
%  clear TempFrames
% end
%  legend('Nose','Implant','Body','Tail','CenterBottomRight','CenterTopLeft')
%  clear Frame2Plot
% %% show the corresponding frames
% VideoInfo=VideoReader([EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\Video\',EphysObj.Video.VideoAVI]);EphysObj.Video.VideoInfo=VideoInfo;
%   try
% VideoInfo=VideoReader([EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\Video\',EphysObj.Video.VideoAVI]);EphysObj.Video.VideoInfo=VideoInfo;
% figure; count=1;
% for MotifNumber=1:length(AllFrames)-1
%     TempFrames=AllFrames{MotifNumber, 1} ; 
% for Frame2Plot=SelectedFrames(MotifNumber,1:5)
%     subplot(length(AllFrames)-1,5,count)
%     imshow(read(VideoInfo,Frame2Plot))
% title([num2str(EphysObj.Video.Times(Frame2Plot))])
% count=count+1;
% end
% clear TempFrames
% end
%   catch
%       VideoInfo=VideoReader([EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\Video\',EphysObj.Video.VideoAVI]);EphysObj.Video.VideoInfo=VideoInfo;
% figure; count=1;
% for MotifNumber=1:length(AllFrames)-1
%     TempFrames=AllFrames{MotifNumber, 1} ; 
% for Frame2Plot=SelectedFrames(MotifNumber,1:5) 
%     subplot(length(AllFrames)-1,5,count)
%     imshow(read(VideoInfo,Frame2Plot))
% title([num2str(EphysObj.Video.Times(Frame2Plot))])
% count=count+1;
% end
% clear TempFrames
% end
%   end
% clear Frame2Plot
% end %if
% end % for plot
%%
for GetFR=1:1
%% collect spikes
SpikeData=EphysObj.Video.NormData;SpikeData=double(SpikeData)/60000000; %data in minutes
%ListOfMotifs={Close2TLFrames,Close2BRFrames,RearingFrames,WalkingFrames,RunningFrames,StoppingFrames,RightTurnFrames,LeftTurnFrames,RestOfFrames};
try SensorTimes=double((EphysObj.TTL(1).on-double(EphysObj.raw_data.TimeStamp(1))))/60000000;catch; end %timestamps in seconds
if EphysObj.Variables.FoodBR
MotifsNames={'Empty','Food','Rearing','Walking','Running','Stopping','RightTurn','LeftTurn','other','Sensor','NoSensor'};
else
MotifsNames={'Food','Empty','Rearing','Walking','Running','Stopping','RightTurn','LeftTurn','other','Sensor','NoSensor'};
end
% create a table to store the FR in each event
EachEventFR=nan(100,11);
% get the bouts for each motif
for SelectedMotif=1:length(EphysObj.Video.AllFrames)+1
    try
    if SelectedMotif<length(EphysObj.Video.AllFrames)+1
MotifTimes=EphysObj.Video.AllFrames{SelectedMotif,1 };
MotifTimes=EphysObj.Video.Times(MotifTimes);
Motifs(SelectedMotif).Name=char(MotifsNames{1,SelectedMotif});
    elseif SelectedMotif==length(EphysObj.Video.AllFrames)+1 % add the sensor data
        try
MotifTimes=double(SensorTimes);
Motifs(SelectedMotif).Name=char(MotifsNames{1,SelectedMotif});
        catch
        end
    end
%% calculate time diffrance between adjecent frames
TimeIntervals=([MotifTimes(2:end),MotifTimes(end)]-MotifTimes);% calculate the intervals between bout start-start
BoatsStartTimesFrame=[1 find(TimeIntervals>MinimalBoutTime/60)+1];% take bouts that have a greater interval thatn 
BoatsStartTimes=MotifTimes(BoatsStartTimesFrame);% get bout start time
BoatsEndTimesFrame=[find(TimeIntervals>MinimalBoutTime/60) length(TimeIntervals)];
BoatsEndTimes=MotifTimes(BoatsEndTimesFrame);% get bout end time
if length(BoatsStartTimesFrame)>1
for i=1:length(BoatsStartTimesFrame)
Bout(i).StartTime=BoatsStartTimes(1,i);% time in minutes
Bout(i).EndTime=BoatsEndTimes(1,i);% time in minutes
Bout(i).BoutLength=Bout(i).EndTime-Bout(i).StartTime;% time in minutes
Bout(i).BoutLengthSeconds=60*double(Bout(i).BoutLength);% time in seconds
% take the spikes in the range of each boat 
SpikeData=EphysObj.Video.NormData';SpikeData=double(SpikeData)/60000000;
Bout(i).SpikesInRange=SpikeData(find(SpikeData>BoatsStartTimes(1,i) & SpikeData<BoatsEndTimes(1,i)));
Bout(i).LocationOfSpikesInRange=find(SpikeData>BoatsStartTimes(1,i) & SpikeData<BoatsEndTimes(1,i));
% calculated frequency in Hz by deviding the total boat length in seconds
Bout(i).FRinBoat=length(Bout(i).SpikesInRange)/Bout(i).BoutLengthSeconds;
% if Bout(i).BoutLengthSeconds<EphysObj.Variables.BoutLengthSecondsLimit; Bout(i).FRinBoat=nan; end % remove the data of too short bouts (less than 1 sec)
end
elseif length(BoatsStartTimesFrame)==1
for i=1:1 % if the first bout is happening when recording starts
Bout(i).StartTime=BoatsStartTimes(1,i);% time in seconds
Bout(i).EndTime=BoatsEndTimes(1,i);% time in seconds
Bout(i).BoutLength=Bout(i).EndTime-Bout(i).StartTime;% time in minutes
Bout(i).BoutLengthSeconds=60*Bout(i).BoutLength;% time in seconds
% take the spikes in the range of each boat 
SpikeData=EphysObj.Video.NormData';SpikeData=double(SpikeData)/60000000;
Bout(i).SpikesInRange=SpikeData(find(SpikeData>BoatsStartTimes(1,i) & SpikeData<BoatsEndTimes(1,i)));
Bout(i).LocationOfSpikesInRange=find(SpikeData>BoatsStartTimes(1,i) & SpikeData<BoatsEndTimes(1,i));
% calculated frequency in Hz by deviding the total boat length in seconds
Bout(i).FRinBoat=length(Bout(i).SpikesInRange)/Bout(i).BoutLengthSeconds;
if Bout(i).BoutLengthSeconds<EphysObj.Variables.BoutLengthSecondsLimit; Bout(i)=[];end
end
else
Bout=nan; SpikeData=nan;  
end

% Make the data into a table
EachEventFR(1:length(Bout),SelectedMotif)=extractfield(Bout,'FRinBoat')';
% continue with the code
Motifs(SelectedMotif).Bout=Bout;
Motifs(SelectedMotif).MeanFR4EachBout=mean([Motifs(SelectedMotif).Bout.FRinBoat],'omitnan');
Motifs(SelectedMotif).SEM = std([Motifs(SelectedMotif).Bout.FRinBoat],'omitnan')/sqrt(length([Motifs(SelectedMotif).Bout.FRinBoat]));  
Motifs(SelectedMotif).NumberOfBouts = length(Motifs(SelectedMotif).Bout);
% collect all the spikes in the range and bout length for selected motife
TotalLength=extractfield(Bout,'BoutLength');% in minutes
Motifs(SelectedMotif).SumOfLengthSeconds=60*sum(TotalLength); %in seconds
Motifs(SelectedMotif).Frames=EphysObj.Video.AllFrames{SelectedMotif,1 }; %in seconds
Motifs(SelectedMotif).TotalFrameTime=0.0667*length(Motifs(SelectedMotif).Frames); %in seconds
Motifs(SelectedMotif).TotalSpikes=extractfield(Bout,'SpikesInRange');% in minutes
Motifs(SelectedMotif).TotalSpikesNumber=length(Motifs(SelectedMotif).TotalSpikes);
Motifs(SelectedMotif).FR=Motifs(SelectedMotif).TotalSpikesNumber/Motifs(SelectedMotif).SumOfLengthSeconds;
clearvars -EXCEPT DisplayPlot ListOfFiles EachEventFR Threshold EphysObj Motifs SelectedMotif SpikeData ListOfMotifs MotifsNames NoSensorTimes SensorTimes MinimalBoutTime
catch
continue
end
end % for selectedmotif
EphysObj.Video.EachEventFR=EachEventFR;
%%
% Add information about times without sensor activation
try
Motifs(SelectedMotif+1).Name='NoSensor';
NoSensorRange=([0 double(EphysObj.Video.SimpleNumbers.VideoLastTTLTimestamp/60)]) ;
large_column_vector = SpikeData; %Unit data
lower_bounds=[Motifs(SelectedMotif).Bout.StartTime];
upper_bounds=[Motifs(SelectedMotif).Bout.EndTime];
SpikesInSensor = large_column_vector(any(large_column_vector >= lower_bounds & large_column_vector <= upper_bounds, 2));
Motifs(SelectedMotif+1).TotalSpikes = setdiff(SpikeData,SpikesInSensor);Motifs(SelectedMotif+1).TotalSpikes=Motifs(SelectedMotif+1).TotalSpikes';
Motifs(SelectedMotif+1).TotalSpikesNumber=length(Motifs(SelectedMotif+1).TotalSpikes);
Motifs(SelectedMotif+1).SumOfLengthSeconds=60*(NoSensorRange(end)-NoSensorRange(1))-sum([Motifs(SelectedMotif).Bout.BoutLength]);
Motifs(SelectedMotif+1).FR=Motifs(SelectedMotif+1).TotalSpikesNumber/Motifs(SelectedMotif+1).SumOfLengthSeconds;
catch
end
EphysObj.Video.Motifs=Motifs;
for ANOVA=1:1
MeanOfFRs=[];STDofFRs=[];CountofFRs=[];Xvalues=[];Yvalues=[];
for i=1:size(EphysObj.Video.EachEventFR,2)-3
MeanOfFRs=[MeanOfFRs mean(EphysObj.Video.EachEventFR(:,i),'omitnan')] ;
STDofFRs=[STDofFRs std(EphysObj.Video.EachEventFR(:,i),'omitnan')] ; 
CountofFRs=[CountofFRs length(find(~isnan(EphysObj.Video.EachEventFR(:,i))))];
Xvalues= [Xvalues i*ones(1,length(find(~isnan(EphysObj.Video.EachEventFR(:,i)))))];
Yvalues=[Yvalues (EphysObj.Video.EachEventFR((~isnan(EphysObj.Video.EachEventFR(:,i))),i))'];
end
figure
b=bar(MeanOfFRs, 'grouped','FaceColor','flat');
hold on
errorbar(b.XEndPoints',MeanOfFRs,STDofFRs,'k','LineWidth', 1 ,'linestyle','none','HandleVisibility','off');
if EphysObj.Variables.FoodBR
set(gca,'xticklabel',{'Empty','Food','Rearing','Walking','Running','Stopping','RightTurn','LeftTurn'});
else 
set(gca,'xticklabel',{'Food','Empty','Rearing','Walking','Running','Stopping','RightTurn','LeftTurn'});
end
xtickangle(90)
ylabel(['Firing rate'],'FontSize', 12)
hold on
scatter(Xvalues,Yvalues,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[0 0 0])
figure
[EphysObj.Video.p,tbl,EphysObj.Video.stats] = anova1(EphysObj.Video.EachEventFR(:,[1:8]));
% [EphysObj.Video.p,tbl,EphysObj.Video.stats] = kruskalwallis(EphysObj.Video.EachEventFR(:,[1:8]));
results = multcompare(EphysObj.Video.stats);
MotifAverage=mean(EphysObj.Video.EachEventFR,'omitnan');
EphysObj.Video.OneWayANOVA = array2table(results,"VariableNames",["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
writematrix(EphysObj.Video.EachEventFR,[EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\',EphysObj.Variables.FolderName,'\',EphysObj.Variables.UnitName,'_OneWayANOVA.xls'],'sheet','Each Event FR');
writetable(EphysObj.Video.OneWayANOVA,[EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\',EphysObj.Variables.FolderName,'\',EphysObj.Variables.UnitName,'_OneWayANOVA.xls'],'sheet','One-Way ANOVA');
writematrix(MotifAverage,[EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\',EphysObj.Variables.FolderName,'\',EphysObj.Variables.UnitName,'_OneWayANOVA.xls'],'sheet','Motif Averages');
saveas(b,[EphysObj.Variables.ComputerDir,'\',EphysObj.Variables.Date,'\',EphysObj.Variables.MouseName,'\',EphysObj.Variables.FolderName,'\',EphysObj.Variables.UnitName,'_Bar.jpg']);
writematrix(EphysObj.Video.EachEventFR,['D:\SummaryData\',EphysObj.Variables.UnitName,'_OneWayANOVA.xls'],'sheet','EachEventFR');
writematrix(MotifAverage,['D:\SummaryData\',EphysObj.Variables.UnitName,'_OneWayANOVA.xls'],'sheet','Motif Averages');
writetable(EphysObj.Video.OneWayANOVA,['D:\SummaryData\',EphysObj.Variables.UnitName,'_OneWayANOVA.xls'],'sheet','OneWayANOVA');
saveas(b,['D:\SummaryData\',EphysObj.Variables.UnitName,'_Bar.jpg']);
end
 for Plot=1:1
     if DisplayPlot
% plot times on a pie graph      
% PieX=extractfield(EphysObj.Video.Motifs,'TotalFrameTime'); 
% PieX=PieX([1 2 3 6 4 5 7 8]);
% PieLabels={'Close2TL';'Close2BR';'Rearing';'Stopping';'Walking';'Running';'RightTurn';'LeftTurn';'other'};
% pie(PieX,PieLabels)       
% % plot bar graph        
% figure;bar([1:(length(EphysObj.Video.Motifs)-1)],[EphysObj.Video.Motifs.MeanFR4EachBout]);hold on
% er = errorbar([1:(length(EphysObj.Video.Motifs)-1)],[EphysObj.Video.Motifs.MeanFR4EachBout],[EphysObj.Video.Motifs.SEM],[EphysObj.Video.Motifs.SEM]);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';  
% set(gca,'xticklabel',extractfield(EphysObj.Video.Motifs  ,'Name'))
% ylabel('Unit FR (Hz)')
% title('Average FR in each bout')
% % plot general average FR
% figure;bar([1:length(EphysObj.Video.Motifs)],[EphysObj.Video.Motifs.FR]);hold on
% set(gca,'xticklabel',extractfield(EphysObj.Video.Motifs  ,'Name'))
% ylabel('Unit FR (Hz)')
% title('Average FR overall')
 end
 end % for plot
end %for GetFR
%% 
% for GetCrossCorrelation=1:1
%     % correlation of FR and sensor rate
% [Sensorc,Sensorlags] = xcorr(EphysObj.Video.DataNPerFrame,EphysObj.Video.SensorNPerFrame);
%     % correlation of FR and velocity of mouse
% [Velocityc,Velocitylags] = xcorr(EphysObj.Video.DataNPerFrame,EphysObj.Video.BodyDistanceBetweenFrames);
% end
end % end function   
                
        function [ArrayForScatter,CollectFrequency]=BoutAnalysis(EphysObj)
%% Get infirmation about bout start and end
Raster=EphysObj.Video.Motifs;
MotifNames=extractfield(Raster,'Name');
% Get the cell timestamp
CellTimeStamp=(double(EphysObj.timestamps(EphysObj.cells==EphysObj.Variables.UnitNumber)-double(EphysObj.raw_data.TimeStamp(1))))/1000000; %timestamp in microseconds
%% get these parameters for later plots
BinEdges4Hz = linspace(-EphysObj.Variables.WindowSizeSec,EphysObj.Variables.WindowSizeSec,2*EphysObj.Variables.WindowSizeSec+1);
BinCenters = (BinEdges4Hz(1:end-1) + BinEdges4Hz(2:end)) / 2;
NumberOfBinsHz=2*EphysObj.Variables.WindowSizeSec;
if EphysObj.Variables.RasterAll
    EphysObj.Variables.Motif2Analyse=[1:8 10];
elseif EphysObj.Variables.FoodBR
    EphysObj.Variables.Motif2Analyse=2;
else
    EphysObj.Variables.Motif2Analyse=1;
end
for Motif=EphysObj.Variables.Motif2Analyse%
        RasterFigure=figure;
for m=1:2
    if m==2
        EphysObj.Variables.PlotOnset=false;
    end
subplot(3, 2, m);
TimesStartSec=extractfield(Raster(Motif).Bout,'StartTime')*60;  
TimesEndSec=extractfield(Raster(Motif).Bout,'EndTime')*60; 
BoutLengthSec=extractfield(Raster(Motif).Bout,'BoutLength')*60; 
% take out bouts that are smaller than Minimal bout length
TimesStartSec=TimesStartSec(BoutLengthSec>EphysObj.Variables.BoutLengthSecondsLimit);
TimesEndSec=TimesEndSec(BoutLengthSec>EphysObj.Variables.BoutLengthSecondsLimit); 
BoutLengthSec=BoutLengthSec(BoutLengthSec>EphysObj.Variables.BoutLengthSecondsLimit); 
% sort by motif length
if EphysObj.Variables.SortByMotifLength
[~,BoutLengthSecOrder]=sort(BoutLengthSec,'descend');
else
  BoutLengthSecOrder=[1:length(TimesStartSec)]  ;
end
TimesStartSec=TimesStartSec(BoutLengthSecOrder);
TimesEndSec=TimesEndSec(BoutLengthSecOrder);
if EphysObj.Variables.PlotOnset
Time2Use=TimesStartSec;
else
Time2Use=TimesEndSec;
end
% collect cell activity by Start stamps
ArrayForScatter=nan(length(Time2Use),length(CellTimeStamp));
for BoutNumber=BoutLengthSecOrder
Minimal2Collect=Time2Use(BoutNumber)-double(EphysObj.Variables.WindowSizeSec);
Maximal2Collect=Time2Use(BoutNumber)+double(EphysObj.Variables.WindowSizeSec);
ActivityInRange=CellTimeStamp(Minimal2Collect<CellTimeStamp & CellTimeStamp<Maximal2Collect);
ActivityInRangeZeroed=ActivityInRange-Time2Use(BoutNumber);
ArrayForScatter(BoutNumber,1:length(ActivityInRangeZeroed))=ActivityInRangeZeroed;
% plot all the timestamps bined to 1 sec
scatter(ActivityInRangeZeroed,BoutNumber*ones(length(ActivityInRange),1),10,'|','MarkerFaceColor','k','MarkerEdgeColor','k',...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);hold on
% Add bout start/end times
if EphysObj.Variables.PlotOnset
scatter(TimesEndSec(BoutNumber)-TimesStartSec(BoutNumber),BoutNumber,'r','|');hold on
else
scatter(TimesStartSec(BoutNumber)-TimesEndSec(BoutNumber),BoutNumber,'r','|');hold on
end
end; clear Minimal2Collect Maximal2Collect ActivityInRange ActivityInRangeZeroed
xlabel('Time from event onset (Sec)'); ylabel('Event'); xline(0);
if m==2
    title('offset');
else
    title('onset');
end
xlim([-EphysObj.Variables.WindowSizeSec EphysObj.Variables.WindowSizeSec])
%% Make a histogram of each event
subplot(3, 2, m+4);
% show each individule trace
for i=BoutLengthSecOrder 
[N,~] = histcounts(ArrayForScatter(i,:),NumberOfBinsHz); 
% plot(BinCenters,N,"Color",[0, 0, 0, 0.2]); hold on ; 
CollectFrequency(i,:)=N;
clear N
end 
xlim([-EphysObj.Variables.WindowSizeSec EphysObj.Variables.WindowSizeSec])
% Add the average trace
[Sum,~] = histcounts(ArrayForScatter,NumberOfBinsHz); hold on
% AverageFrequency=double(Sum/size(ArrayForScatter,1));
if EphysObj.Variables.PlotZscore
 CollectFrequency=zscore(CollectFrequency,1,'all');
end
data_Frequency=mean(CollectFrequency,1);
data_SEM = std(CollectFrequency,1)/sqrt(size(CollectFrequency,1));       % SEM Across Columns
plot(BinCenters,data_Frequency,'Color','k','LineWidth', 2);hold on
plot(BinCenters,data_Frequency+data_SEM, 'k--', 'LineWidth', 0.5);hold on
plot(BinCenters,data_Frequency-data_SEM, 'k--', 'LineWidth', 0.5);
xlabel('Time from event onset (Sec)'); ylabel('Firing rate (hz)'); xline(0);  
%% Make a heatmap of each event
subplot(3, 2, m+2);
hHM=heatmap(CollectFrequency,'Colormap',jet,'CellLabelColor','none'); %'ColorLimits',[0 20]
hHM.NodeChildren(3).YDir='normal';  
sgtitle(char(MotifNames(Motif))) 
end
%% save the figure
saveas(RasterFigure,['C:\Users\netas\OneDrive\Documents\Obeisty paper\EXPERIMENTS\Single unit\UnitImages\',EphysObj.Variables.Date,'_',EphysObj.Variables.MouseName,'_Tetrode_',num2str(EphysObj.Variables.TetrodeNumber),'_Unit_',num2str(EphysObj.Variables.UnitNumber),'_File_0',num2str(EphysObj.Variables.FileNumber),char(MotifNames(Motif)),'.jpg'])
end
%% Piezo bouts



%% make an array for all the FR in bouts
Data=nan(9,60);
% Get the FR of each event that is grater than minimal bout time
figure
for MotifNumber=[1:8 10]
    if MotifNumber<10
    Data(MotifNumber,1:length(EphysObj.Video.Motifs(MotifNumber).Bout))=extractfield(EphysObj.Video.Motifs(MotifNumber).Bout,'FRinBoat')';
    else
    Data(MotifNumber-1,1:length(EphysObj.Video.Motifs(MotifNumber).Bout))=extractfield(EphysObj.Video.Motifs(MotifNumber).Bout,'FRinBoat')';
    end
    EphysObj.Video.Motifs.AverageHeatmapData=mean(Data,2,'omitnan');
h = heatmap(Data,'Colormap',jet, 'CellLabelColor','none','MissingDataColor', [1,1,1]);
xlabel('Event #'); ylabel('Motif');
h.YDisplayLabels=MotifNames([1:8 10]);
end
% save the average dat to obj 
EphysObj.Video.AverageHeatmapData=mean(Data,2,'omitnan');
EphysObj.Video.AverageHeatmapDataZ=mean(zscore(Data,1,'omitnan'),2,'omitnan');
EphysObj.Video.AverageHeatmapDataZ=zscore(EphysObj.Video.AverageHeatmapData);
    end

end
end 
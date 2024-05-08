function [TrialArray] = AcutFeedingVideoAnalyse(DLCFileFolder,FileName,DisplayPlot,Threshold) 
MouseName=FileName(1:3);
%% other variables
timestamp_interval = 5;%Sec - this is the interval of the TTL signal from the bioobserve to NLX
MinimalBoutTime=0.3;
%% start code
% collect data
% read the dlc (csv) file
FullFileName=[DLCFileFolder,'\',FileName];
Info=readtable(FullFileName);
InfoArray=table2array(Info);
% read the video file
VideoFileName=[FullFileName(1:(strfind(FullFileName,'DLC')-1)),'.avi'];
try
VideoFileNameLabeled=[FileName(1:end-4),'_labeled.mp4'];
catch
VideoFileNameLabeled=VideoFileName;
end
VideoInfo=VideoReader(VideoFileName);
FrameTimes=linspace(0,VideoInfo.Duration,VideoInfo.NumFrames); % X values for frame basd analysis
FrameEdgesSec=linspace(0,FrameTimes(end)*2-mean([FrameTimes(end),FrameTimes(end-1)]),VideoInfo.NumFrames+1);
FrameEdgesMin=FrameEdgesSec/60;
%% work on the DeepLabCut output
% Take out data that has low accuracy
Temp=InfoArray(:,2:3);Temp(InfoArray(:,4)<0.85,1:2)=nan;Struct.Nose=Temp;clear Temp 
Temp=InfoArray(:,5:6);Temp(InfoArray(:,7)<0.85,1:2)=nan;Struct.Head=Temp;clear Temp
Temp=InfoArray(:,8:9);Temp(InfoArray(:,10)<0.85,1:2)=nan;Struct.Body=Temp;clear Temp
Temp=InfoArray(:,11:12);Temp(InfoArray(:,13)<0.85,1:2)=nan;Struct.TailBase=Temp;clear Temp
Temp=InfoArray(:,14:15);Temp(InfoArray(:,16)<0.85,1:2)=nan;Struct.Tail=Temp;clear Temp
Temp=InfoArray(:,17:18);Temp(InfoArray(:,19)<0.85,1:2)=nan;Struct.Jelly=Temp;clear Temp
Temp=InfoArray(:,20:21);Temp(InfoArray(:,22)<0.85,1:2)=nan;Struct.Empty=Temp;clear Temp
% collect data by frames
NewArray=nan(length(InfoArray),14);
NewArray(:,1:2)=Struct.Nose;
NewArray(:,3:4)=Struct.Head;
NewArray(:,5:6)=Struct.Body;
NewArray(:,7:8)=Struct.TailBase;
NewArray(:,9:10)=Struct.Tail;
NewArray(:,11:12)=Struct.Jelly;
NewArray(:,13:14)=Struct.Empty;
% define the arena
MaxX=max([NewArray(:,[1 3 5 7 9 11 13])]);MaxX=max(MaxX);
MinX=min([NewArray(:,[1 3 5 7 9 11 13])]);MinX=min(MinX);
MaxY=max(NewArray(:,[2 4 6 8 10 12 14]));MaxY=max(MaxY);
MinY=min(NewArray(:,[2 4 6 8 10 12 14]));MinY=min(MinY);
ArenaLength=MaxX-MinX;
ArenaLengthFactor=25/ArenaLength;
ArenaWidth=MaxY-MinY;
ArenaWidthFactor=25/ArenaWidth;
ArenaFactor=mean([ArenaLengthFactor,ArenaWidthFactor]);
ArenaAll(1,1:2)=[MinX MaxY];Arena.BottomLeft=[MinX MaxY];
ArenaAll(3,1:2)=[MinX MinY];Arena.TopLeft=[MinX MinY];
ArenaAll(2,1:2)=[MaxX MaxY];Arena.BottomRight=[MaxX MaxY];
ArenaAll(4,1:2)=[MaxX MinY];Arena.TopRight=[MaxX MinY];

%% calculate all the distances between the markers. 
NewArrayX=NewArray(:,[1:2:13]); % get x values of array
NewArrayY=NewArray(:,[2:2:14]); % get Y values of array
%% calculate distances between frames
fields = fieldnames(Struct);
for FieldNumber=1:length(fields)
DistancesArray(FieldNumber).DistanceBetweenFrames=...
sqrt((([NewArrayX(2:end,FieldNumber);0]-NewArrayX(:,FieldNumber)).^2)...
+ (([NewArrayY(2:end,FieldNumber);0]-NewArrayY(:,FieldNumber)).^2));
% take out unreasnable distances
NewArrayX(DistancesArray(FieldNumber).DistanceBetweenFrames*ArenaFactor*15>70,FieldNumber)=nan;
NewArrayY(DistancesArray(FieldNumber).DistanceBetweenFrames*ArenaFactor*15>70,FieldNumber)=nan;
DistancesArray(FieldNumber).DistanceBetweenFrames(DistancesArray(FieldNumber).DistanceBetweenFrames*ArenaFactor*15>70)=nan;
end
% update the NewArray variable to exclude weird values
NewArray(:,[1:2:13])=NewArrayX; % get x values of array
NewArray(:,[2:2:14])=NewArrayY; % get Y values of array
%% calculate distances between points
for FieldNumber=1:length(fields)
DistancesArray(FieldNumber).FieldName=fields(FieldNumber);
DistancesArray(FieldNumber).Distances=sqrt((NewArrayX-NewArrayX(:,FieldNumber)).^2 + (NewArrayY-NewArrayY(:,FieldNumber)).^2);
DistancesArray(FieldNumber).Distances=DistancesArray(FieldNumber).Distances(:,1:length(fields));
end
%%
for Plot=1:1
% plot each body part over time
    figure
    Xframes=[1:1:VideoInfo.NumFrames];
    for BodyPart=1:length(DistancesArray)
plot(Xframes,DistancesArray(BodyPart).DistanceBetweenFrames);hold on
    end   
xlabel('Frame Number')
ylabel('Distance between parts')  
%% Plot the general track
% plot each body part in the arena
     figure
    for BodyPart=[2 3 4] % 1:length(DistancesArray)
plot(NewArrayX(:,BodyPart),-1*NewArrayY(:,BodyPart));hold on
    end   
xlim([MinX MaxX]); ylim([-MaxY -MinY]);
% legend(char(fields(1)),char(fields(2)),char(fields(3)),char(fields(4)),char(fields(5)),char(fields(6)),char(fields(7)))
    if DisplayPlot
end % if
end %for plot
%% figure out the implant position relative to the plates
Times=FrameTimes/60;
% collect the times where there is threshold crossing
AllFrames=[1:VideoInfo.NumFrames];
% rearing is distance of head (2) and body (3) 
RearingFrames=find([DistancesArray(2).Distances(:,3)]<(Threshold/2)); % small distance of body and head indicates rearing
%% calculate angele between points
for GetAngles=1:1
XYBody = [NewArray(:,5) NewArray(:,6)]; % body
XYHead = [NewArray(:,3) NewArray(:,4)]; % head
XYTailBase = [NewArray(:,7) NewArray(:,8)]; % tail base
VectorHB=(XYHead-XYBody); VectorBT=(XYBody-XYTailBase);
for Temp=1:length(VectorHB)
BodyHeadAngle(Temp)= (180/pi)*atan2(VectorHB(Temp,1)*VectorBT(Temp,2)-VectorBT(Temp,1)*VectorHB(Temp,2),dot(VectorHB(Temp,:),VectorBT(Temp,:)));
end
BodyHeadAngle=BodyHeadAngle';
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
% get a sample of optional angles
Angles2Test=linspace(-180,180,13);
BodyHeadAngleNORearing=BodyHeadAngle;
BodyHeadAngleNORearing(RearingFrames)=nan;
for AngleRange=1:length(Angles2Test)-1
GroupAngle(AngleRange).Frame=(find(Angles2Test(AngleRange)<BodyHeadAngleNORearing&BodyHeadAngleNORearing<Angles2Test(AngleRange+1)));
GroupAngle(AngleRange).Values=(BodyHeadAngleNORearing(GroupAngle(AngleRange).Frame));
end
end %GetAngles
LeftTurnFrames=setdiff(GroupAngle(9).Frame,RearingFrames);%90
RightTurnFrames=setdiff(GroupAngle(4).Frame,RearingFrames);%-90
%% find frames of plates approaching
FoodPlate=[nanmedian(NewArrayX(:,6)),nanmedian(NewArrayY(:,6))];
FoodPlate = ginput(1);
FoodPlate(2)=FoodPlate(2)*-1;
%
EmptyPlate=[nanmedian(NewArrayX(:,7)),nanmedian(NewArrayY(:,7))];
% deal with the case it can't find the empty plate
% if isnan(EmptyPlate)
    try
EmptyPlate = ginput(1);
EmptyPlate(2)=EmptyPlate(2)*-1;
    catch
EmptyPlate = input("What is xy value? enter in the format [X Y] ");        
    end
% end
for FieldNumber=1:length(fields)
DistancesArray(FieldNumber).Distance2Empty=sqrt((EmptyPlate(1)-NewArrayX(:,FieldNumber)).^2 + (EmptyPlate(2)-NewArrayY(:,FieldNumber)).^2);
DistancesArray(FieldNumber).Distance2Food=sqrt((FoodPlate(1)-NewArrayX(:,FieldNumber)).^2 + (FoodPlate(2)-NewArrayY(:,FieldNumber)).^2);
end
hold on; scatter(round(EmptyPlate(1)),-1*round(EmptyPlate(2)),'k');
hold on; scatter(round(FoodPlate(1)),-1*round(FoodPlate(2)),'k');
figure
% plot distance of head and food
if DisplayPlot
figure
plot(Xframes,DistancesArray(2).Distance2Food);hold on
plot(Xframes,DistancesArray(2).Distance2Empty );hold on 
title ('distance to food/empty')
end
Close2FoodFrames=find([DistancesArray(2).Distance2Food(:,1)]<Threshold); % small distance indicates proximity to plate
Close2FoodFrames=setdiff(Close2FoodFrames,[RearingFrames;LeftTurnFrames;RightTurnFrames]); 
Close2EmptyFrames=find([DistancesArray(2).Distance2Empty(:,1)]<Threshold); % small distance indicates proximity to plate
Close2EmptyFrames=setdiff(Close2EmptyFrames,[RearingFrames;LeftTurnFrames;RightTurnFrames]); 
%% determine velocity between frames
for GetVelocity=1:1
BodyDistanceBetweenFrames=DistancesArray(3).DistanceBetweenFrames  ;
BodyDistanceBetweenFramesNorm=BodyDistanceBetweenFrames*ArenaFactor*15; % cm/sec  %15 is the sampling rate
% define the frames
BodyDistanceBetweenFramesMoving=BodyDistanceBetweenFramesNorm;
ExlcludeFrames=[RearingFrames;RightTurnFrames;LeftTurnFrames;Close2EmptyFrames;Close2FoodFrames];
BodyDistanceBetweenFramesMoving(ExlcludeFrames)=nan;
StoppingFrames=find(BodyDistanceBetweenFramesMoving<1); % threshold for low velocity
WalkingFrames=find(BodyDistanceBetweenFramesMoving>=1&BodyDistanceBetweenFramesMoving<10);
SpeedingFrames=find(BodyDistanceBetweenFramesNorm>=10&BodyDistanceBetweenFramesMoving<30);% threshold for high velocity
RunningFrames=find(BodyDistanceBetweenFramesNorm>=30&BodyDistanceBetweenFramesMoving<70);% threshold for high velocity
%%
RestOfFrames=setdiff(AllFrames,[Close2EmptyFrames;Close2FoodFrames;RearingFrames;StoppingFrames;WalkingFrames;SpeedingFrames;RunningFrames;RightTurnFrames;LeftTurnFrames]);
end
%% collect all frames of each motif
Motifes(1).Name={'Rearing'};
Motifes(1).Frames=RearingFrames;
Motifes(2).Name={'Stopping'};
Motifes(2).Frames=StoppingFrames;
Motifes(3).Name={'Walking'};
Motifes(3).Frames=WalkingFrames;
Motifes(4).Name={'Speeding'};
Motifes(4).Frames=SpeedingFrames;
Motifes(5).Name={'Running'};
Motifes(5).Frames=RunningFrames;
Motifes(6).Name={'RightTurn'};
Motifes(6).Frames=RightTurnFrames;
Motifes(7).Name={'LeftTurn'};
Motifes(7).Frames=LeftTurnFrames;
Motifes(8).Name={'Close2FoodFrames'};
Motifes(8).Frames=Close2FoodFrames;
Motifes(9).Name={'Close2EmptyFrames'};
Motifes(9).Frames=Close2EmptyFrames;
Motifes(10).Name={'Other'};
Motifes(10).Frames=RestOfFrames;
% collect times
for MotifNumber=1:length(Motifes)
Motifes(MotifNumber).Times=Times(Motifes(MotifNumber).Frames);
end
 for Plot=1:1
     DisplayMorePlot=false; 
% plot the activity on one timeline
%% get the relevant times from the frames
RearingTimes=Times(RearingFrames); 
StoppingTimes=Times(StoppingFrames);
WalkingTimes=Times(WalkingFrames);
SpeedingTimes=Times(SpeedingFrames);
RunningTimes=Times(RunningFrames);
RightTurnTimes=Times(RightTurnFrames);
LeftTurnTimes=Times(LeftTurnFrames);
Close2FoodTimes=Times(Close2FoodFrames);
Close2EmptyTimes=Times(Close2EmptyFrames);
RestOfFramesTimes=Times(RestOfFrames);
%plot
BehaviorFigure=figure;
scatter(StoppingTimes,zeros(1,length(StoppingTimes)),'|','MarkerEdgeColor',[0.231 0.741 0.871]);hold on
scatter(WalkingTimes,zeros(1,length(WalkingTimes)),'|','MarkerEdgeColor',[0.98 0.592 0]);hold on
scatter(SpeedingTimes,zeros(1,length(SpeedingTimes)),'|','MarkerEdgeColor',[0.45 0.652 0.1]);hold on
scatter(RunningTimes,zeros(1,length(RunningTimes)),'|','MarkerEdgeColor',[0.369 0.408 0.8]);hold on
scatter(RightTurnTimes,zeros(1,length(RightTurnTimes)),'|','MarkerEdgeColor',[0.094 0.871 0.224]);hold on
scatter(LeftTurnTimes,zeros(1,length(LeftTurnTimes)),'|','MarkerEdgeColor',[0.859 0.631 0.173]);hold on
scatter(RearingTimes,zeros(1,length(RearingTimes)),'|', 'MarkerEdgeColor',[0.447 0.569 0.475]);hold on
scatter(Close2EmptyTimes,zeros(1,length(Close2EmptyTimes)),'|', 'MarkerEdgeColor',[ 1 0.922 0]);hold on
scatter(Close2FoodTimes,zeros(1,length(Close2FoodTimes)),'|','MarkerEdgeColor',[0.957 0 1]);hold on
plot(Times,(BodyDistanceBetweenFramesNorm),'-k');hold on
% plot(Times,(smoothdata(BodyDistanceBetweenFramesNorm,'gaussian',100)));hold on
ylim([0 60])
% xlim([0 1])
legend('Stopping','Walking','Speeding','Running'...
    ,'RightTurn','LeftTurn','Rearing','Empty','Food','Velocity')
xlabel('Time (minutes)'); 
for Trial=[15:15:45]
xline([Trial],'--g');
end
if DisplayMorePlot
% collect random frames from each motif
for MotifNumber=1:length(Motifes)-1
    try
    Motifes(MotifNumber).SelectedFrames=[randsample(Motifes(MotifNumber).Frames,5)];
    catch
    try
    Motifes(MotifNumber).SelectedFrames=[randsample(Motifes(MotifNumber).Frames,4)];
    catch
    try
    Motifes(MotifNumber).SelectedFrames=[randsample(Motifes(MotifNumber).Frames,3)];
    catch
    try
    Motifes(MotifNumber).SelectedFrames=[randsample(Motifes(MotifNumber).Frames,2)];
    catch
    end
    end
    end
    end
end
% show the labeled frames
try
VideoInfo=VideoReader([DLCFileFolder,'\',VideoFileNameLabeled]);
catch
VideoInfo=VideoReader([VideoFileName]);
end
figure;count=1;
for MotifNumber=1:length(Motifes)-1
subplot(3,3,count)
imshow(read(VideoInfo,Motifes(MotifNumber).SelectedFrames(1)));
title(Motifes(MotifNumber).Name)
count=count+1;
end
% show the velocity
figure;
plot(Times,(BodyDistanceBetweenFramesNorm));hold on
title ('Velocity')
ylim([0 70])
xlim([0 2])
%% Plot a selected frame with markers
% plot the markers
count=1; figure; MarkerSize=2;
for MotifNumber=1:length(Motifes)-1
 for Frame2Plot=Motifes(MotifNumber).SelectedFrames'
subplot(length(Motifes)-1,5,count)
scatter(Struct.Nose(Frame2Plot,1),-1*Struct.Nose(Frame2Plot,2),MarkerSize);hold on
scatter(Struct.Jelly(Frame2Plot,1),-1*Struct.Jelly(Frame2Plot,2),MarkerSize);hold on
scatter(Struct.Empty(Frame2Plot,1),-1*Struct.Empty(Frame2Plot,2),MarkerSize);hold on
scatter(Struct.Head(Frame2Plot,1),-1*Struct.Head(Frame2Plot,2),MarkerSize);hold on
scatter(Struct.Body(Frame2Plot,1),-1*Struct.Body(Frame2Plot,2),MarkerSize);hold on
scatter(Struct.TailBase(Frame2Plot,1),-1*Struct.TailBase(Frame2Plot,2),MarkerSize);hold on
scatter(Struct.Tail(Frame2Plot,1),-1*Struct.Tail(Frame2Plot,2),MarkerSize);hold on
xlim([MinX MaxX]); ylim([-MaxY -MinY]);
title(num2str(Frame2Plot))
count=count+1;
 end
end
 legend('Nose','Jelly','Empty','Head','Body','TailBase','Tail')
% show the corresponding frames
VideoInfo=VideoReader([VideoFileName]);
try
VideoInfo=VideoReader([VideoFileName]);
figure; count=1;
for MotifNumber=1:length(Motifes)-1
for Frame2Plot=Motifes(MotifNumber).SelectedFrames'
subplot(length(Motifes)-1,5,count)
imshow(read(VideoInfo,Frame2Plot))
title(num2str(Frame2Plot))
count=count+1;
end
end
catch
VideoInfo=VideoReader([VideoFileName]);
figure; count=1;
for MotifNumber=1:length(Motifes)-1
for Frame2Plot=Motifes(MotifNumber).SelectedFrames'
subplot(length(Motifes)-1,5,count)
imshow(read(VideoInfo,Frame2Plot))
title(num2str(Frame2Plot))
count=count+1;
end
end
end
end %if
end % for plot
%% get the bouts for each motif
MotifTimes=60*Motifes(MotifNumber).Times; % time in seconds
%% work on each part of the timeline
Markers=linspace(1,VideoInfo.NumFrames,5);
for MotifNumber=1:length(Motifes)
TrialArray(MotifNumber).Name=Motifes(MotifNumber).Name;
TrialArray(MotifNumber).PreFeeding=Motifes(MotifNumber).Frames(Motifes(MotifNumber).Frames>Markers(1)&Motifes(MotifNumber).Frames<Markers(2));
TrialArray(MotifNumber).OFF1=Motifes(MotifNumber).Frames(Motifes(MotifNumber).Frames>Markers(2)&Motifes(MotifNumber).Frames<Markers(3));
TrialArray(MotifNumber).ON=Motifes(MotifNumber).Frames(Motifes(MotifNumber).Frames>Markers(3)&Motifes(MotifNumber).Frames<Markers(4));
TrialArray(MotifNumber).OFF2=Motifes(MotifNumber).Frames(Motifes(MotifNumber).Frames>Markers(4)&Motifes(MotifNumber).Frames<Markers(5));
TrialArray(MotifNumber).PreFeedingTimes=0.0667*length(TrialArray(MotifNumber).PreFeeding);
TrialArray(MotifNumber).OFF1Times=0.0667*length(TrialArray(MotifNumber).OFF1);
TrialArray(MotifNumber).ONTimes=0.0667*length(TrialArray(MotifNumber).ON);
TrialArray(MotifNumber).OFF2Times=0.0667*length(TrialArray(MotifNumber).OFF2);
end
UnitData.TrialArray=TrialArray;
UnitData.Motifes=Motifes;
UnitData.Arena=Arena;
UnitData.Motifes=Motifes;
UnitData.FileName=FileName;
UnitData.Struct=Struct;
UnitData.DistancesArray=DistancesArray;
UnitData.MouseName=MouseName;
save(['C:\Users\netas\Documents\Lammel lab\Obesity paper\EXPERIMENTS\NTS-OE\DLC\',MouseName,'.mat'],'UnitData')

%% plot the time spent in each motif
for Plot=1:1
%     if DisplayPlot
%         figure
% PieX=extractfield(Motifes,'TotalFrameTime'); 
% PieLabels={'Rearing';'Stopping';'Walking';'Speeding';'Running';'RightTurn';'LeftTurn'};
% pie(PieX,PieLabels) 
%     end
end % for plot
end % end function   

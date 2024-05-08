function[DistanceFromFood,DistanceFromEmpty,NewArray] = BehaviorTracking(FileName)
FileName='20200221_SUBLAT13-3_Jelly-ExposureDLC_resnet50_Jelly feeding SUBLATMay17shuffle1_100000_filtered.csv';
Info=readtable(FileName);
InfoArray=table2array(Info);
% read the video file
VideoInfo=VideoReader('20200221_SUBLAT13-3_Jelly-Exposure.avi');
% Take out data that has low accuracy
Temp=InfoArray(:,2:3);Temp(InfoArray(:,4)<0.85,1:2)=nan;Struct.Nose=Temp;clear Temp
Temp=InfoArray(:,5:6);Temp(InfoArray(:,7)<0.85,1:2)=nan;Struct.Implant=Temp;clear Temp
Temp=InfoArray(:,8:9);Temp(InfoArray(:,10)<0.85,1:2)=nan;Struct.Body=Temp;clear Temp
Temp=InfoArray(:,11:12);Temp(InfoArray(:,13)<0.85,1:2)=nan;Struct.Tail=Temp;clear Temp
Temp=InfoArray(:,14:15);Temp(InfoArray(:,16)<0.85,1:2)=nan;Struct.Empty=Temp;clear Temp
Temp=InfoArray(:,17:18);Temp(InfoArray(:,19)<0.85,1:2)=nan;Struct.Food=Temp;clear Temp
% collect data by frames
NewArray=nan(length(InfoArray),12);
NewArray(:,1:2)=Struct.Nose;
NewArray(:,3:4)=Struct.Implant;
NewArray(:,5:6)=Struct.Body;
NewArray(:,7:8)=Struct.Tail;
NewArray(:,9:10)=Struct.Empty;
NewArray(:,11:12)=Struct.Food;
%% define the arena
MaxX=max(NewArray(:,[1 3 5 7 9 11]));MaxX=max(MaxX);
MaxY=max(NewArray(:,[2 4 6 8 10 12]));MaxY=max(MaxY);
ArenaAll(1,1:2)=[0 MaxY];Arena.TopLeft=[0 MaxY];
ArenaAll(2,1:2)=[MaxX MaxY];Arena.TopRight=[MaxX MaxY];
ArenaAll(3,1:2)=[0 0];Arena.BottomLeft=[0 0];
ArenaAll(4,1:2)=[MaxX 0];Arena.BottomRight=[MaxX 0];
% define Food location
CenterFoodX=nanmedian(Struct.Food(:,1));
CenterFoodY=nanmedian(Struct.Food(:,2));
CenterEmptyX=nanmedian(Struct.Empty(:,1)); %%CenterEmptyX=MaxX-nanmedian(Struct.Food(:,1));
CenterEmptyY=nanmedian(Struct.Empty(:,2)); %%CenterEmptyY=MaxY-nanmedian(Struct.Food(:,2));
% calculate the distance of head from Food and empty in each frame
BodyPartToMeasure=Struct.Body;
figure;
DistanceFromFood = sqrt((BodyPartToMeasure(:,1)-CenterFoodX).^2 + (BodyPartToMeasure(:,1)-CenterFoodY).^2);
DistanceFromEmpty = sqrt((BodyPartToMeasure(:,1)-CenterEmptyX).^2 + (BodyPartToMeasure(:,1)-CenterEmptyY).^2);
plot(-1*DistanceFromFood);hold on; plot(-1*DistanceFromEmpty)
%%
%% Plot all the data
figure 
plot(Struct.Nose(:,1),-1*Struct.Nose(:,2));hold on
plot(Struct.Implant(:,1),-1*Struct.Implant(:,2));hold on
plot(Struct.Body(:,1),-1*Struct.Body(:,2));hold on
plot(Struct.Tail(:,1),-1*Struct.Tail(:,2));hold on
plot(Struct.Empty(:,1),-1*Struct.Empty(:,2));hold on
plot(Struct.Food(:,1),-1*Struct.Food(:,2));hold on
% xlim([0 MaxX]); ylim([0 -MaxY]);
legend('Nose','Implant','Body','Tail','Empty','Food')
%% Plot a selected frame
count=1;
figure
 for Frame=1:750:length(Struct.Nose)
subplot(5,5,count)
scatter(Struct.Nose(Frame,1),-1*Struct.Nose(Frame,2));hold on
scatter(Struct.Implant(Frame,1),-1*Struct.Implant(Frame,2));hold on
scatter(Struct.Body(Frame,1),-1*Struct.Body(Frame,2));hold on
scatter(Struct.Tail(Frame,1),-1*Struct.Tail(Frame,2));hold on
scatter(Struct.Empty(Frame,1),-1*Struct.Empty(Frame,2));hold on
scatter(Struct.Food(Frame,1),-1*Struct.Food(Frame,2));hold on
scatter(nanmedian(Struct.Food(:,1)),-1*nanmedian(Struct.Food(:,2)));hold on
scatter(CenterEmptyX,-1*CenterEmptyY);hold on
% scatter(ArenaAll(:,1),-1*ArenaAll(:,2))
% limit to Arena size
xlim([0 MaxX]); ylim([-1*MaxY 0]);
count=count+1;
 end
 legend('Nose','Implant','Body','Tail','Empty','Food','CenterFood','CenterEmpty','Arena')
 % show the corresponding frames
count=1;figure
for Frame=[1:750:length(Struct.Nose)]
 subplot(5,5,count)
    imshow(read(VideoInfo,Frame))
    count=count+1;
end
  end
% Icropped = imcrop(imshow(read(VideoInfo,Frame)))
% subplot(5,5,count)
% imshow(Icropped)

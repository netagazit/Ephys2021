% load files
TypeOfCell='Tagged';
TypeOfTest='kruskalwallis';%kruskalwallis
DirectoryName=['D:\SummaryMay2024\',TypeOfCell,'\',TypeOfTest];
Directory=dir(DirectoryName); Directory=extractfield(Directory,'name')';Directory=Directory(3:end);
CollecPvaluMatrix=zeros(8,8);
NotSignificant=0;
Significnant=0;
for i=1:length(Directory)
FileName= [DirectoryName,'\',Directory{i, 1}]  ;
load (FileName);
if Obj2Save.kruskalwallis_p<0.05
CollecPvaluMatrix=CollecPvaluMatrix+Obj2Save.PvalueMatrixSignificant;
Significnant=Significnant+1;
else 
NotSignificant=NotSignificant+1;
end
clear Obj2Save
end 
for i=1:8
CollecPvaluMatrix(i,i)=nan;
end
MotifsNames={'Empty','Food','Rearing','Walking','Running','Stopping','RightTurn','LeftTurn'};
CollecPvaluMatrix=round((double(CollecPvaluMatrix)*100/double(Significnant)));
h=heatmap(CollecPvaluMatrix, 'MissingDataColor', [1,1,1]);%'CellLabelColor','none'
% xlabel('Motif'); ylabel('Event #');
h.XDisplayLabels=MotifsNames;
h.YDisplayLabels=MotifsNames;
title([TypeOfCell,' ',TypeOfTest,' % cells ','Significantly responsive ',num2str(double(Significnant/(Significnant+NotSignificant)))]);
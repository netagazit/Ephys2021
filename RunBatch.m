% load files
DirectoryName='D:\SummaryMay2024';
Directory=dir(DirectoryName); Directory=extractfield(Directory,'name')';Directory=Directory(3:end);
CollecPvaluMatrix=zeros(8,8);
for i=1:length(Directory)
FileName= [DirectoryName,'\',Directory{i, 1}]  ;
load (FileName);
CollecPvaluMatrix=CollecPvaluMatrix+Obj2Save.PvalueMatrixSignificant  ;
clear Obj2Save
end; 
h=heatmap(CollecPvaluMatrix,'CellLabelColor','none');
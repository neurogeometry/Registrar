function FeatureExtraction(StackList_csv_pth)

StackList = table2cell(readtable(StackList_csv_pth,'Delimiter',','));
[PathStr,FolderName]=fileparts(StackList_csv_pth);
DataFolder=[PathStr,'/Results-',FolderName];
Folder = strcat(PathStr,'/Results-',FolderName);
newFolder = strcat(PathStr,'/Results-',FolderName,'/tmp/');
mkdir(newFolder);

% For Windows/Linux
% overlaps_csv = cell2mat(table2cell(readtable([DataFolder,'/overlaps.csv'],'Delimiter',',')));
% StackPositions_pixels_csv = cell2mat(table2cell(readtable([DataFolder,'/StackPositions_pixels.csv'],'Delimiter',',')));
% StackSizes_pixels_csv = cell2mat(table2cell(readtable([DataFolder,'/StackSizes_pixels.csv'],'Delimiter',',')));

% For Mac/Unix
overlaps_csv = cell2mat(table2cell(readtable(strcat(Folder,'/overlaps.csv'),'Delimiter',',')));
StackPositions_pixels_csv = cell2mat(table2cell(readtable(strcat(Folder,'/StackPositions_pixels.csv'),'Delimiter',',')));
StackSizes_pixels_csv = cell2mat(table2cell(readtable(strcat(Folder,'/StackSizes_pixels.csv'),'Delimiter',',')));

for i=1:size(StackList,1)
%     For Windows/Linux
%     tifFile = [PathStr,'\',char(StackList(i,1))];

%   For Mac/Unix
    tifFile = [PathStr,'/',char(StackList(i,1))];
    overlap_ind = overlaps_csv(i,find(overlaps_csv(i,:))); 
    % Run Feature Extraction for one Stack
    IM_Original=ImportStack(char(tifFile(3)),StackSizes_pixels_csv(i,:));
    FeatureExtractionFunc(IM_Original,i,Folder,StackPositions_pixels_csv(overlap_ind,:),StackSizes_pixels_csv(overlap_ind,:));
end
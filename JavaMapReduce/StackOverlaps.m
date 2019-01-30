function StackOverlaps(StackList_csv_pth)

StackList = table2cell(readtable(StackList_csv_pth,'Delimiter',','));
[PathStr,FolderName]=fileparts(StackList_csv_pth);
DataFolder=[PathStr,'/Results-',FolderName];
% index to remove invalid files
errIndxs = [];
StackSizes_pixels = zeros(size(StackList,1),3);
for i = 1:size(StackList,1)
    try
        [filepath,~,ext] = fileparts(char(StackList(i,1)));
        if size(ext,2)>1
            InfoImage=imfinfo(char(StackList(i,1)));
        else
            allfiles = dir(filepath);
            [~,~,ext] = fileparts(allfiles(1).name);
        end
    catch
        StackList{i} = [PathStr,'\',StackList{i}];
    end
    
    try
        [filepath,~,ext] = fileparts(char(StackList(i,1)));
        if size(ext,2)>1
            InfoImage=imfinfo(char(StackList(i,1)));
            StackSizes_pixels(i,1) = InfoImage.Height;
            StackSizes_pixels(i,2) = InfoImage.Width;
            StackSizes_pixels(i,3) = size(InfoImage,1);
        else
            allfiles = dir(filepath);
            NumFiles=length(allfiles);
            imgIdx = 0;
            for j = 1:NumFiles
                [~,~,ext] = fileparts(allfiles(j).name);
                if strcmpi(ext,'.tif') || strcmpi(ext,'.jp2') || strcmpi(ext,'.png') || strcmpi(ext,'.jpeg')
                    InfoImage=imfinfo(char([allfiles(j).folder,'\',allfiles(j).name]));
                    StackSizes_pixels(i,1) = InfoImage.Height;
                    StackSizes_pixels(i,2) = InfoImage.Width;
                    imgIdx = imgIdx +1;
                end
            end
            StackSizes_pixels(i,3) = imgIdx;
        end
    catch
        disp(['Can not read the file: ',char(StackList(i,1))]);
        errIndxs = [errIndxs, i];
    end
end
StackList (errIndxs,:) = [];
StackSizes_pixels (errIndxs,:) = [];
if size(StackList,2) == 7
    StackSizes_pixels = cell2mat(StackList(:,5:7));
end

% Get positions from the list
StackPositions_pixels = cell2mat(StackList(:,2:4)); %xlsread(StackPositions_pixels_csv_pth);
StackPositions_pixels = StackPositions_pixels(:,[2,1,3]);

All_overlaps = FindOverlaps(StackPositions_pixels,StackSizes_pixels);
overlaps = [];
for i=1:size(StackList,1)
    % Find overlap indexes
    overlap_ind = [i,find(All_overlaps(i,:)),find(All_overlaps(:,i))'];
    overlaps(i,1:size(overlap_ind,2))=overlap_ind;
end

dlmwrite([DataFolder,'\overlaps.csv'],overlaps);
dlmwrite([DataFolder,'\StackPositions_pixels.csv'],StackPositions_pixels);
dlmwrite([DataFolder,'\StackSizes_pixels.csv'],StackSizes_pixels);

function Stack = JPG2Mat(Path,StackSizes_pixels)
% close all
% Stack_id=1;
% % StackList_csv_pth='C:\Armen\Publications\Paper35 (Registration Seyed)\MicroscopeFiles\Results-Neocortical2_StackList\JPG\Neocortical2_StackList_JPG.csv';
% StackList_csv_pth='C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\Results-MouseLight_StackList\JPG\MouseLight_StackList_JPG.csv';
% StackList = table2cell(readtable(StackList_csv_pth,'Delimiter',','));
% 
% StackInfo = StackList(Stack_id,:);

FileList = dir(Path);
FileList = {FileList(3:end).name};

Stack=zeros(StackSizes_pixels,'uint8');
count=0;

for i = 1:size(FileList,2)

    temp=imread([Path,FileList{i}],'jpg');
    n_planes=size(temp,2)/StackSizes_pixels(2);
    Stack(:,:,count+1:count+n_planes) = reshape(temp,StackSizes_pixels(1),StackSizes_pixels(2),n_planes);
    count=count+n_planes;

%     figure,imshow(max(Stack,[],3));




%     fileID = fopen([StackInfo{1},FileList{i}]);
%     data = fread(fileID);
%     fclose(fileID);
%     assert (isa(uint8(data), 'uint8'));
%     jImg = javax.imageio.ImageIO.read(java.io.ByteArrayInputStream(data));
% %     h = jImg.getHeight;
%     w = jImg.getWidth;
%     temp = reshape(jImg.getData.getDataStorage, [StackInfo{5},w]);
%     
% 
%     n_planes=size(temp,2)/StackInfo{6};
%     Stack(:,:,count+1:count+n_planes) = reshape(temp,StackInfo{5},StackInfo{6},n_planes);
%     count=count+n_planes;

%     figure,imshow(max(Stack,[],3));
end

% figure,imshow(max(Stack,[],3));

% for i = 1:size(FileList,2)
%     fileID = fopen([StackInfo{1},FileList{i}],'rb');
%     data = fread(fileID, Inf, '*uint8');
%     fclose(fileID);
%     jImg = javax.imageio.ImageIO.read(java.io.ByteArrayInputStream(data));
% %     h = jImg.getHeight;
%     w = jImg.getWidth;
%     temp = reshape(jImg.getData.getDataStorage, [StackInfo{5},w]);
%     n_planes=size(temp,2)/StackInfo{6};
%     Stack(:,:,count+1:count+n_planes) = reshape(temp,StackInfo{5},StackInfo{6},n_planes);
%     count=count+n_planes;
% end


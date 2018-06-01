function Stack=JPG2Mat()

Stack_id=1;
StackList_csv_pth='C:\Armen\Publications\Paper35 (Registration Seyed)\MicroscopeFiles\Results-Neocortical2_StackList\JPG\Neocortical2_StackList_JPG.csv';
StackList = table2cell(readtable(StackList_csv_pth,'Delimiter',','));

StackInfo = StackList(Stack_id,:);

FileList = dir(StackInfo{1});
FileList = {FileList(3:end).name};

IM=zeros(StackInfo{5},StackInfo{6},StackInfo{7},'uint8');
count=0;
for i = 1:size(FileList,1)
    temp=imread([StackInfo{1},FileList{i}]);
    n_planes=size(temp,2)/StackInfo{6};
    IM(:,:,count+1:count+n_planes) = reshape(temp,StackInfo{5},StackInfo{6},n_planes);
    count=count+n_planes;
end
figure,imshow(max(IM,[],3));


function Stack2JPG()

 StackList_csv_pth='C:\Armen\Publications\Paper35 (Registration Seyed)\MicroscopeFiles\Neocortical2_StackList.csv';
%StackList_csv_pth='C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\MouseLight_StackList.csv';

StackList = table2cell(readtable(StackList_csv_pth,'Delimiter',','));

[FilePath,FileName,~]=fileparts(StackList_csv_pth);
SaveFolder = [FilePath,'\Results-',FileName];
mkdir([SaveFolder,'\JPG'])
M=cell(size(StackList,1),7);
M(:,2:4)=StackList(:,2:4);
for i = 1:size(StackList,1)
    %     try
    [filepath,filename,ext] = fileparts(char(StackList(i,1)));
    temp=strsplit(filepath,'\');
    StackName=temp{end};
    if size(ext,2)>1
        InfoImage=imfinfo(char(StackList(i,1)));
        StackSizes_pixels(1) = InfoImage.Height;
        StackSizes_pixels(2) = InfoImage.Width;
        StackSizes_pixels(3) = size(InfoImage,1);
    else
        allfiles = dir(filepath);
        NumFiles=length(allfiles);
        imgIdx = 0;
        for j = 1:NumFiles
            [~,~,ext] = fileparts(allfiles(j).name);
            if strcmpi(ext,'.tif') || strcmpi(ext,'.jp2') || strcmpi(ext,'.png') || strcmpi(ext,'.jpeg')
                InfoImage=imfinfo(char([allfiles(j).folder,'\',allfiles(j).name]));
                StackSizes_pixels(1) = InfoImage.Height;
                StackSizes_pixels(2) = InfoImage.Width;
                imgIdx = imgIdx +1;
            end
        end
        StackSizes_pixels(3) = imgIdx;
        filepath = [filepath,'\'];
    end
    MaxIntensityValue = InfoImage(1).MaxSampleValue;
    M{i,1}=[SaveFolder,'\JPG\',StackName,'\'];
    M{i,5}=StackSizes_pixels(1);
    M{i,6}=StackSizes_pixels(2);
    M{i,7}=StackSizes_pixels(3);
    mkdir([SaveFolder,'\JPG\',StackName]);
            IM_Original=ImportStack(filepath,StackSizes_pixels);
    % IM_Original=ImportStack([filepath,'\',filename,'.tif'],StackSizes_pixels);
    IM_Original = uint8(double(IM_Original)./double(MaxIntensityValue)*255);
%     figure,imshow(max(IM_Original,[],3));
    N_planes=fix(8192/StackSizes_pixels(2));%65535
    Nfiles=ceil(StackSizes_pixels(3)/N_planes);
    for j=1:Nfiles-1
        temp=reshape(IM_Original(:,:,N_planes*(j-1)+1:N_planes*j),StackSizes_pixels(1),StackSizes_pixels(2)*N_planes);
        imwrite(temp,[[SaveFolder,'\JPG\',StackName,'\'],num2str(j),'.jpg'],'jpg');
    end
    res=StackSizes_pixels(3)-N_planes*(Nfiles-1);
    temp=reshape(IM_Original(:,:,end-res+1:end),StackSizes_pixels(1),StackSizes_pixels(2)*res);
    imwrite(temp,[SaveFolder,'\JPG\',StackName,'\',num2str(Nfiles),'.jpg'],'jpg');
    
    %     catch
    %         disp(['Can not read the file: ',char(StackList(i,1))]);
    %     end
end

writetable(cell2table(M),[SaveFolder,'\JPG\',FileName,'_JPG','.csv'],'WriteVariableNames',false)


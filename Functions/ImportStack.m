function FinalImage=ImportStack(File,StackSizes_pixels)
warning('off','all');



[filepath,~,ext] = fileparts(char(File));
if size(ext,2)>1
    TifLink = Tiff(File, 'r');
    InfoImage=imfinfo(File);
    FinalImage=zeros(StackSizes_pixels,['uint',num2str(InfoImage(1).BitsPerSample)]);
    for i=1:StackSizes_pixels(3)
        TifLink.setDirectory(i);
        FinalImage(:,:,i)=TifLink.read();
    end
    TifLink.close();
else
    allfiles = dir(filepath);
    InfoImage=imfinfo([File,'\',allfiles(3).name]);
    FinalImage=zeros(StackSizes_pixels,['uint',num2str(InfoImage(1).BitsPerSample)]);
    for i = 1:StackSizes_pixels(3)
        [~,~,ext] = fileparts(allfiles(i).name);
        if strcmp(ext,'.tif') || strcmp(ext,'.jp2') || strcmp(ext,'.png') || strcmp(ext,'.jpeg')
            FinalImage(:,:,i) = imread(char([allfiles(i).folder,'\',allfiles(i).name]));
        end
    end
end

end


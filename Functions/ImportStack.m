function FinalImage=ImportStack(Path,StackSizes_pixels)
% warning('off','all');
% InfoImage=imfinfo(FileTif);
% mImage=InfoImage(1).Width;
% nImage=InfoImage(1).Height;
% NumberImages=length(InfoImage);
% FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
% TifLink = Tiff(FileTif, 'r');
% for i=1:NumberImages
%    TifLink.setDirectory(i);
%    FinalImage(:,:,i)=TifLink.read();
% end
% TifLink.close();

warning('off','all');



[filepath,~,ext] = fileparts(char(Path));
if size(ext,2)>1
    TifLink = Tiff(Path, 'r');
    InfoImage=imfinfo(Path);
    FinalImage=zeros(StackSizes_pixels,['uint',num2str(InfoImage(1).BitsPerSample)]);
    for i=1:StackSizes_pixels(3)
        TifLink.setDirectory(i);
        FinalImage(:,:,i)=TifLink.read();
    end
    TifLink.close();
else
    allfiles = dir(filepath);
    InfoImage=imfinfo([Path,'/',allfiles(3).name]);
    if strcmp(InfoImage.Format,'jpg') && StackSizes_pixels(3) > 1
%         FinalImage=zeros(StackSizes_pixels,['uint',num2str(InfoImage(1).BitDepth)]);
        FinalImage = JPG2Mat(Path,StackSizes_pixels);
    else
        FinalImage=zeros(StackSizes_pixels,['uint',num2str(InfoImage(1).BitsPerSample)]);
        for i = 1:StackSizes_pixels(3)
            [~,~,ext] = fileparts(allfiles(i).name);
            if strcmp(ext,'.tif') || strcmp(ext,'.jp2') || strcmp(ext,'.png') || strcmp(ext,'.jpeg') || strcmp(ext,'.jpg')
                FinalImage(:,:,i) = imread(char([allfiles(i).folder,'/',allfiles(i).name]));
            end
        end
    end
    
end

end


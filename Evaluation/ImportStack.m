function FinalImage=ImportStack(File)
warning('off','all');

[filepath,~,ext] = fileparts(char(File));
if size(ext,2)>1
    InfoImage=imfinfo(File);
    mImage=InfoImage(1).Width;
    nImage=InfoImage(1).Height;
    NumberImages=length(InfoImage);
 %    FinalImage=zeros(nImage,mImage,NumberImages,'uint8');
    FinalImage=zeros(nImage,mImage,NumberImages,['uint',num2str(InfoImage(1).BitsPerSample)]);
    TifLink = Tiff(File, 'r');
    for i=1:NumberImages
        TifLink.setDirectory(i);
        FinalImage(:,:,i)=TifLink.read();
    end
    TifLink.close();
else
%     FinalImage = [];
%     FinalImage=zeros(nImage,mImage,NumberImages,['uint',num2str(InfoImage(1).BitsPerSample)]);
    allfiles = dir(filepath);
    NumFiles=length(allfiles);
    j = 1;
    for i = 1:NumFiles
        [~,~,ext] = fileparts(allfiles(i).name);
        if strcmp(ext,'.tif') || strcmp(ext,'.jp2') || strcmp(ext,'.png') || strcmp(ext,'.jpeg')
            FinalImage(:,:,j) = imread(char([allfiles(i).folder,'\',allfiles(i).name]));
            j = j + 1;
        end
    end
    
end





% clc
% tic
% InfoImage=imfinfo(FileTif);
% mImage=InfoImage(1).Width;
% nImage=InfoImage(1).Height;
% NumberImages=length(InfoImage);
% FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
%
% for i=1:NumberImages
%     FinalImage(:,:,i)=imread(FileTif,'Index',i,'Info',InfoImage);
% end
% first = toc

% tic
% InfoImage=imfinfo(FileTif);
% mImage=InfoImage(1).Width;
% nImage=InfoImage(1).Height;
% NumberImages=length(InfoImage);
% FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
% FileID = tifflib('open',FileTif,'r');
% rps = tifflib('getField',FileID,Tiff.TagID.RowsPerStrip);
% stripNum = 1;
% for i=1:NumberImages-1
%    tifflib('setDirectory',FileID,i);
%    % Go through each strip of data.
%    rps = min(rps,nImage);
%    for r = 1:rps:nImage
%       row_inds = r:min(nImage,r+rps-1);
% %        stripNum = tifflib('computeStrip',FileID,r,i);
%
%       FinalImage(row_inds,:,i) = tifflib('readEncodedStrip',FileID,int8(stripNum));
%    end
% end
% tifflib('close',FileID);
end


function stitched = StitichAll(Image1,Image2,Image3,Image4,Dataset,P2,P3,P4)
% close all;

% Dataset = 'Neuromuscular1';

% Image1 = 106;
% Image2 = 107;
% Image3 = 114;
% Image4 = 115;

if strcmp(Dataset,'Mine')
    load('../data/StackData.mat');
elseif strcmp(Dataset, 'DIADEM1')
    load('../data/StackData_DIADEM.mat');
elseif strcmp(Dataset, 'DIADEM2')
    load('../data/StackData_DIADEM2.mat');
elseif strcmp(Dataset, 'Svoboda')
    load('../data/StackData_Svoboda.mat');
elseif strcmp(Dataset, 'Neuromuscular1')
    load('../data/StackData_Neuromuscular1.mat');
elseif strcmp(Dataset, 'Holtmaat_2_1')
    load('../data/StackData_Holtmaat_2_1.mat');
elseif strcmp(Dataset, 'Visual')
    load('../data/StackData_Visual.mat');
end

if strcmp(Dataset, 'Holtmaat_2_1')|| strcmp(Dataset, 'Visual')|| strcmp(Dataset, 'Neuromuscular1')
%     Source_Stack_Folder = [char(StackList(SourceID,2)),'/'];
%     Source_Stack_File = [char(StackList(SourceID,1)),'.tif'];
    IM1=ImportStack([[char(StackList(Image1,2)),'/'],[char(StackList(Image1,1)),'.tif']]);
    IM2=ImportStack([[char(StackList(Image2,2)),'/'],[char(StackList(Image2,1)),'.tif']]);
     IM3=ImportStack([[char(StackList(Image3,2)),'/'],[char(StackList(Image3,1)),'.tif']]);
     IM4=ImportStack([[char(StackList(Image4,2)),'/'],[char(StackList(Image4,1)),'.tif']]);
else
%     tifFile_1 = dir([char(StackList(SourceID,2)) '/*.tif' ]);
%     Source_Stack_Folder = [tifFile_1(1).folder,'/'];
%     Source_Stack_File = tifFile_1(1).name;
%     IM_Source=ImportStack([Source_Stack_Folder,Source_Stack_File]);
IM1=ImportStack([[char(StackList(Image1,2)),'/'],[char(StackList(Image1,1)),'-ngc.0.tif']]);
    IM2=ImportStack([[char(StackList(Image2,2)),'/'],[char(StackList(Image2,1)),'-ngc.0.tif']]);
     IM3=ImportStack([[char(StackList(Image3,2)),'/'],[char(StackList(Image3,1)),'-ngc.0.tif']]);
     IM4=ImportStack([[char(StackList(Image4,2)),'/'],[char(StackList(Image4,1)),'-ngc.0.tif']]);
end



%     figure,imshow(max(IM1,[],3),[0 max(IM1(:))]);
%     figure,imshow(max(IM2,[],3),[0 max(IM2(:))]);
%      figure,imshow(max(IM3,[],3),[0 max(IM3(:))]);
%      figure,imshow(max(IM4,[],3),[0 max(IM4(:))]);

[Y1,X1,Z1]=size(IM1);
[Y2,X2,Z2]=size(IM2);
[Y3,X3,Z3]=size(IM3);
[Y4,X4,Z4]=size(IM4);
% ---------------------------- 1-2
    DX = P2(1)-X1;%
    DY = P2(2);%
    if DX < 0
        IM_Target_tmp = padarray(IM2,[0 X1+DX],0,'pre');
        IM_Source_tmp = padarray(IM1,[0 X1+DX],0,'post');
    else
        IM_Target_tmp = padarray(IM2,[0 X1+DX],0,'post');
        IM_Source_tmp = padarray(IM1,[0 X1+DX],0,'pre');
    end
    
    if DY > 0
        IM_Target_tmp = padarray(IM_Target_tmp,[abs(DY) 0],0,'post');
        IM_Source_tmp = padarray(IM_Source_tmp,[abs(DY) 0],0,'pre');
    else
        IM_Target_tmp = padarray(IM_Target_tmp,[abs(DY) 0],0,'pre');
        IM_Source_tmp = padarray(IM_Source_tmp,[abs(DY) 0],0,'post');
    end
    sourceMax = max(IM_Source_tmp,[],3);
    targetMax = max(IM_Target_tmp,[],3);
    %     figure(17),imshow(max(IM_Source_tmp,[],3),[0 M_source]) ;
    %     figure(18),imshow(max(IM_Target_tmp,[],3),[0 M_source]) ;
    IM_OUT1 = [];
%     IM_OUT1 = zeros(Y1+abs(DY),X1+DX+X2);
    IM_OUT1(:,:,:) = max(IM_Target_tmp(:,:,:),IM_Source_tmp(:,:,:));
%      figure(),imshow(IM_OUT1,[0 max(IM_OUT1(:))]) ;
%      figure(),imshow(max(IM_OUT1,[],3),[0 max(IM_OUT1(:))]);
    
% ----------------------------- 2 -3
    DX = P3(1)-X3;
    %
    DY = P3(2)+1;%
    if DX < 0
        IM_Target_tmp1 = padarray(IM4,[0 X1+DX],0,'pre');
        IM_Source_tmp1 = padarray(IM3,[0 X1+DX],0,'post');
    else
        IM_Target_tmp1 = padarray(IM4,[0 X1+DX],0,'post');
        IM_Source_tmp1 = padarray(IM3,[0 X1+DX],0,'pre');
    end
    
    if DY > 0
        IM_Target_tmp1 = padarray(IM_Target_tmp1,[abs(DY) 0],0,'post');
        IM_Source_tmp1 = padarray(IM_Source_tmp1,[abs(DY) 0],0,'pre');
    else
        IM_Target_tmp1 = padarray(IM_Target_tmp1,[abs(DY) 0],0,'pre');
        IM_Source_tmp1 = padarray(IM_Source_tmp1,[abs(DY) 0],0,'post');
    end
    sourceMax1 = max(IM_Source_tmp1,[],3);
    targetMax1 = max(IM_Target_tmp1,[],3);
    %     figure(17),imshow(max(IM_Source_tmp,[],3),[0 M_source]) ;
    %     figure(18),imshow(max(IM_Target_tmp,[],3),[0 M_source]) ;
    IM_OUT2 = [];
%     IM_OUT2 = zeros(Y3+abs(DY),X3+DX+X4);
    IM_OUT2(:,:,:) = max(IM_Target_tmp1(:,:,:),IM_Source_tmp1(:,:,:));
%      figure(),imshow(IM_OUT2,[0 max(IM_OUT2(:))]) ;
%  figure(),imshow(max(IM_OUT2,[],3),[0 max(IM_OUT2(:))]);
% ----------------------------- all
    

[Y1,X1,Z1]=size(IM_OUT1);
[Y2,X2,Z2]=size(IM_OUT2);
IM_OUT2 = imresize(IM_OUT2,[size(IM_OUT1,1) size(IM_OUT1,2)],'bicubic'); % or 'linear', etc.
    DX = P4(1)-1;
    DY = P4(2)-Y2+1;
    
%     IM_OUT = zeros(Y1-abs(DY)+Y2,X1+abs(DX),Z1);
    IM_Target_tmp = [];
    IM_Source_tmp = [];
    if DX > 0
        IM_Target_tmp = padarray(IM_OUT2,[0 abs(DX) ],0,'pre');
        IM_Source_tmp = padarray(IM_OUT1,[0 abs(DX) ],0,'post');
    else
        IM_Target_tmp = padarray(IM_OUT2,[0 abs(DX) ],0,'post');
        IM_Source_tmp = padarray(IM_OUT1,[0 abs(DX) ],0,'pre');
    end
    
    if DY < 0
        IM_Target_tmp = padarray(IM_Target_tmp,[Y1-abs(DY) 0],0,'pre');
        IM_Source_tmp = padarray(IM_Source_tmp,[Y1-abs(DY) 0],0,'post');
    else
        IM_Target_tmp = padarray(IM_Target_tmp,[Y1-abs(DY) 0],0,'post');
        IM_Source_tmp = padarray(IM_Source_tmp,[Y1-abs(DY) 0],0,'pre');
    end
    
%     figure(17),imshow(max(IM_Source_tmp,[],3),[0 M_source]) ;
%     figure(18),imshow(max(IM_Target_tmp,[],3),[0 M_source]) ;
    stitched = [];
%     IM_OUT = zeros(Y1+abs(DY),X1+DX+X2);
    %IM_Target_tmp = IM_Source_tmp(IM_Source_tmp>0);
    stitched(:,:,:) = max(IM_Source_tmp(:,:,:),IM_Target_tmp(:,:,:));  
%     Error_diff = sum(sum(sum(abs(IM_Source_tmp - IM_Target_tmp))))
%      figure(),imshow(max(IM_OUT,[],3),[0 max(IM_OUT(:))]) ;
%       figure(),imshow(max(stitched,[],3),[0 max(stitched(:))]);






end
function [IM_Source, IM_Target, IM_source_max, IM_target_max, M_source] = showFeatures1(StackList,SourceID,TargetID,Source_seed,Target_seed)

    
    Source_Stack_File = [StackList{SourceID,2},'\',StackList{SourceID,1},'.tif'];
    
    Target_Stack_File = [StackList{TargetID,2},'\',StackList{TargetID,1},'.tif'];

%     IM_Source1=ImportStack('E:\Datasets\TimeLaps\Images\DL083C001G.tif');
%     IM_Source2=ImportStack('E:\Datasets\TimeLaps\Images\DL083D001G.tif');
%     IM_Source3=ImportStack('E:\Datasets\TimeLaps\Images\DL083E001G.tif');
%     IM_Source4=ImportStack('E:\Datasets\TimeLaps\Images\DL083F001G.tif');
%     IM_Source6=ImportStack('E:\Datasets\TimeLaps\Images\DL083G001G.tif');
%     IM_Source7=ImportStack('E:\Datasets\TimeLaps\Images\DL083H001G.tif');
%     f1 = figure;imshowpair(max(IM_Source1,[],3),max(IM_Source2,[],3),'falsecolor','Scaling','joint');caxis([0 500])
%     f2 = imshowpair(max(f1.CData,[],3),max(IM_Source3,[],3),'falsecolor','Scaling','joint');caxis([0 500])
%     figure;imshowpair(max(f2.CData,[],3),max(IM_Source4,[],3),'falsecolor','Scaling','joint');caxis([0 500])

    IM_Source=ImportStack(char(Source_Stack_File));

    IM_Target=ImportStack(char(Target_Stack_File));


[X1,Y1,Z1]=size(IM_Source);
[X2,Y2,Z2]=size(IM_Target);
%     IM_Source1 = IM_Source./max(IM_Target(:))+eps;
% if ~strcmp(Dataset, 'Svoboda')
%     if strcmp(Dataset, 'Holtmaat_2_1')
%         figure,imshow(max(IM_Source,[],3),[0 500]);
%     else
        figure,imshow(max(IM_Source,[],3),[0 max(IM_Source(:))]);
%     end
% else
%     figure,imshow(max(im2double(IM_Source),[],3),[0 1]);
% end

hold on
plot(Source_seed(:,1),Source_seed(:,2),'r*');

% if ~strcmp(Dataset, 'Svoboda')
%     if strcmp(Dataset, 'Holtmaat_2_1')
%         figure,imshow(max(IM_Target,[],3),[0 500]);
%     else
        figure,imshow(max(IM_Target,[],3),[0 max(IM_Target(:))]);
%     end
%     
% else
%     figure,imshow(max(im2double(IM_Target),[],3),[0 1]);
% end

hold on
plot(Target_seed(:,1),Target_seed(:,2),'r*');
IM_source_max=max(IM_Source,[],3);
IM_target_max=max(IM_Target,[],3);
M_source=max(IM_Source(:));
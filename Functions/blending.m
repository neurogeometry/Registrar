function [Tile3D_org,Tile3D,stop] = blending(StackPositions_pixels_original,StackSizes_pixels,StackList,stackID,L,b,DataFolder)
% ============================== About ====================================
% -------------------------------------------------------------------------
%
% Purpose: Main Registration function to call other functions
% Input:
%   k                  : 1*1, The folder address of stack
%   Mesh                    : 1*1, The file name of stack
%   StackList_csv_pth          : 1*1, The address to save the features
%   filterValue          : 1*1, The address to save the features
%   TransformationValue          : 1*1, The address to save the features
%   Seq_Par          : 1*1, The address to save the features
%   Par_workers          : 1*1, The address to save the features
%   blendingSID          : 1*1, The address to save the features
%
% -------------------------------------------------------------------------
% Author: Seyed Mostafa Mousavi Kahaki, Armen Stepanyants
% Northeastern University, USA
% kahaki@neu.edu, a.stepanyants.neu.edu
% =========================================================================
% -------------------------------------------------------------------------
addpath('../Functions');
parameters;
% paramsBPremove_pad = 0;
% axes_xz = findobj(Registrar,'Tag', 'v');
options.overwrite = 1;
tb11 = findobj(Registrar,'Tag', 'pushbutton10');
set(tb11,'userdata',0);
stop = 0;

% allfiles = dir(char(StackList(stackID,1)));
% InfoImage=imfinfo([allfiles(end).folder,'\', allfiles(end).name]);
% cl = ['uint',num2str(InfoImage.BitDepth)];

disp(['Reading ',char(StackList(stackID,1))]);
[~,~,ext] = fileparts(char(StackList(stackID,1)));
if strcmp(ext,'.JPG') || strcmp(ext,'.png')
    IM_Source = imread(char(StackList(stackID,1)));
else
    %     IM_Source = ImportStack(char(StackList(stackID,1)));
    if paramsREuseHDF5
        IM_Source = hdf5read([DataFolder,'\tmp\temp_',num2str(stackID),'.h5'], '/dataset1');
    else
        IM_Source = ImportStack(char(StackList(stackID,1)),StackSizes_pixels(stackID,:));
    end
end
% IM_Source(:,1:3,:) = 60000;
% IM_Source(:,end-2:end,:) = 60000;
% IM_Source(1:3,:,:) = 60000;
% IM_Source(end-2:end,:,:) = 60000;

StackPositions_pixels_original = round(StackPositions_pixels_original);
Neighbors = findStackNeighbors(stackID,StackSizes_pixels,StackPositions_pixels_original);

Tile3D_org=zeros(size(IM_Source)+2*params.BP.extSize,class(IM_Source));
Tile3D_org(params.BP.extSize(1)+1:end-params.BP.extSize(1),params.BP.extSize(2)+1:end-params.BP.extSize(2),params.BP.extSize(3)+1:end-params.BP.extSize(3))=IM_Source;

if size(L,2) > 1
    [IM_Source,StackPosition_reg_source]=Perform_Linear_Transform(IM_Source,StackPositions_pixels_original(stackID,:),L(:,(3*stackID)-2:3*stackID),b(:,stackID));
else
    StackPosition_reg_source = StackPositions_pixels_original(stackID,:)+b(:,stackID)'; %comming stack position
end

Tile3D=zeros(size(IM_Source)+2*params.BP.extSize,class(IM_Source));
Tile3D(params.BP.extSize(1)+1:end-params.BP.extSize(1),params.BP.extSize(2)+1:end-params.BP.extSize(2),params.BP.extSize(3)+1:end-params.BP.extSize(3))=IM_Source;

StackPositions_neighb=StackPositions_pixels_original(Neighbors,:)-ones(length(Neighbors),1)*StackPositions_pixels_original(stackID,:)+params.BP.extSize+1;

for i = 1:size(Neighbors,1)
    tb11 = findobj(Registrar,'Tag', 'pushbutton10');
    if get(tb11,'userdata') || stop% stop condition
        disp(num2str(tb11.UserData));
        disp('Process Stopped');
        stop = 1;
        break;
    end
    SourceID = Neighbors(i);
    disp(['Reading ',char(StackList(SourceID,1))]);
    [~,~,ext] = fileparts(char(StackList(SourceID,1)));
    if strcmp(ext,'.JPG') || strcmp(ext,'.png')
        IM_Source = imread(char(StackList(SourceID,1)));
    else
        %         IM_Source = ImportStack(char(StackList(SourceID,1)));
        if paramsREuseHDF5
            IM_Source = hdf5read([DataFolder,'\tmp\temp_',num2str(SourceID),'.h5'], '/dataset1');
        else
            IM_Source = ImportStack(char(StackList(SourceID,1)),StackSizes_pixels(SourceID,:));
        end
    end

    Start=max([1,1,1;StackPositions_neighb(i,:)]);
    End=min(size(Tile3D_org),StackPositions_neighb(i,:)+size(IM_Source)-1);
    Start1=max([1,1,1;1-StackPositions_neighb(i,:)+1]);
    End1=min(size(IM_Source),size(Tile3D_org)-StackPositions_neighb(i,:)+1);
    Tile3D_org(Start(1):End(1),Start(2):End(2),Start(3):End(3))=max(Tile3D_org(Start(1):End(1),Start(2):End(2),Start(3):End(3)),...
        IM_Source(Start1(1):End1(1),Start1(2):End1(2),Start1(3):End1(3)));
    
    if size(L,2) > 1
        [IM_Source,StackPosition_registered]=Perform_Linear_Transform(IM_Source,StackPositions_pixels_original(SourceID,:),L(:,(3*SourceID)-2:3*SourceID),b(:,SourceID));
    else
        StackPosition_registered = StackPositions_pixels_original(SourceID,:)+b(:,SourceID)'; %comming stack position
    end
    
    StackPosition_registered=StackPosition_registered-StackPosition_reg_source+params.BP.extSize+1;
    
    Start=max([1,1,1;StackPosition_registered]);
    End=min(size(Tile3D),StackPosition_registered+size(IM_Source)-1);
    Start1=max([1,1,1;1-StackPosition_registered+1]);
    End1=min(size(IM_Source),size(Tile3D)-StackPosition_registered+1);
    Tile3D(Start(1):End(1),Start(2):End(2),Start(3):End(3))=max(Tile3D(Start(1):End(1),Start(2):End(2),Start(3):End(3)),...
        IM_Source(Start1(1):End1(1),Start1(2):End1(2),Start1(3):End1(3)));
    
    %               figure();imshow(max(Tile3D,[],3),[0 max(Tile3D(:))]);
    %      figure();imshow(max(IM_Source_Pad_translate,[],3),[0 max(IM_Source_Pad_translate(:))]);
end


if params.BP.saveImages
    saveastiff(Tile3D,[DataFolder,'\stacks_registered.tif'],options);
    saveastiff(Tile3D_org,[DataFolder,'\stacks_before_register.tif'],options);
end

% for Fig 3-A with below parameters    
%Parameters to control the size of the visualization tile size
% params.BP.extSize = [2024 2024 0];
% figure,imshow(max(Tile3D_org(2025:5005,1160:3047,:),[],3) ,[0 max(Tile3D_org(:))]);
end
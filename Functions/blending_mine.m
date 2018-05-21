function [Tile3D_org,Tile3D,stop] = blending_mine(StackPositions_pixels_original,StackSizes_pixels,StackList,stackID,L,b,Seq_Par,Par_workers,DataFolder)
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
paramsBPextSize = [200 200 100];
paramsBPremove_pad = 0;
StackPositions_pixels_original = round(StackPositions_pixels_original);

axes_xz = findobj(NCT_Registration,'Tag', 'v');

tb11 = findobj(NCT_Registration,'Tag', 'pushbutton10');
set(tb11,'userdata',0);
stop = 0;

% if axes_xz.Value == 1
% StackPositions_pixels_org(:,1) = max(StackPositions_pixels_original(:,1))-StackPositions_pixels_original(:,1);
% StackPositions_pixels_org(:,2) = max(StackPositions_pixels_original(:,2))-StackPositions_pixels_original(:,2);
% StackPositions_pixels_org(:,3) = max(StackPositions_pixels_original(:,3))-StackPositions_pixels_original(:,3);
% else
% StackPositions_pixels_org = StackPositions_pixels_original;
% end
StackPositions_pixels_org = StackPositions_pixels_original;
StackPositions_Registered = StackPositions_pixels_org + b';
StackPositions_Registered = round(StackPositions_Registered);

% StackPositions_pixels(:,1) = max(StackPositions_Registered(:,1))-StackPositions_Registered(:,1);
% StackPositions_pixels(:,2) = max(StackPositions_Registered(:,2))-StackPositions_Registered(:,2);
% StackPositions_pixels(:,3) = max(StackPositions_Registered(:,3))-StackPositions_Registered(:,3);

Neighbors = findStackNeighbors(stackID,StackSizes_pixels,StackPositions_Registered);

% tifFile_1 = dir([char(StackList(stackID,2)) '/*.tif' ]);
% Source_Stack_Folder = [tifFile_1(1).folder,'/'];
% Source_Stack_File = tifFile_1(1).name;
% disp(['Reading ',Source_Stack_File]);

% ------------ TEST
% [filepath,name,ext]= fileparts(char(StackList(stackID,1)));
% IM_Sourcemain = ImportStackJ(filepath,{['\',name,ext]},1,1,1,'max');
% IM_Sourcemain = ImportStackJ('E:\Tiling\2014-11-26\00\00735',{['\','00735-ngc.0.tif']},1,1,1,'max');
% -------------------

disp(['Reading ',char(StackList(stackID,1))]);
[filepath,name,ext] = fileparts(char(StackList(stackID,1)));
if strcmp(ext,'.JPG') || strcmp(ext,'.JPEG')
    IM_Sourcemain = imread(char(StackList(stackID,1)));
else
    IM_Sourcemain=ImportStack(char(StackList(stackID,1)));
end

if size(L,2) > 1
    tform = affine3d([[L(:,(3*stackID)-2:3*stackID)';b(:,stackID)'],[0;0;0;1]]);
    IM_Sourcemain = imwarp(IM_Sourcemain,tform);
    %        figure();imshow(max(IM_Sourcemain,[],3),[0 max(IM_Sourcemain(:))]);
end
% IM_Sourcemain = ImportStack(char(StackList(stackID,1)));
% IM_Sourcemain = ImportStack([Source_Stack_Folder,Source_Stack_File]);
%     figure(1);imshow(max(IM_Sourcemain,[],3),[0 max(IM_Sourcemain(:))]);

Tile3D = padarray(IM_Sourcemain,paramsBPextSize);
Tile3D_org = Tile3D;
%     figure(1);imshow(max(Tile3D_org,[],3),[0 max(Tile3D_org(:))]);
% tileSize = size(IM_Pad);
S = StackPositions_Registered(stackID,:); % main stack position
S_org = StackPositions_pixels_org(stackID,:);
% TilePosition =  S  - (extSize/2); % Tile Position
% Tile3D = IM_Pad;
% load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList\MatchedPoints_Translation.mat');


for i = 1:size(Neighbors,1)
    tb11 = findobj(NCT_Registration,'Tag', 'pushbutton10');
    %         disp(num2str(tb11.UserData));
    if get(tb11,'userdata') || stop% stop condition
        disp(num2str(tb11.UserData));
        disp('Process Stopped');
        stop = 1;
        %                 listboxItems{v}  = 'Process Stopped ';
        %                 v = v + 1;
        %                 tb = findobj(NCT_Registration,'Tag', 'listbox1');
        %                 set(tb, 'String', listboxItems);drawnow
        %                 tb.Value = v-1;drawnow
        break;
    end
    SourceID = Neighbors(i);
    %     tifFile_1 = dir([char(StackList(SourceID,2)) '/*.tif' ]);
    %     Source_Stack_Folder = [tifFile_1(1).folder,'/'];
    %     Source_Stack_File = tifFile_1(1).name;
    %     disp(['Reading ',Source_Stack_File]);
    %     IM_Source = ImportStack([Source_Stack_Folder,Source_Stack_File]);
    disp(['Reading ',char(StackList(SourceID,1))]);
    [filepath,name,ext] = fileparts(char(StackList(SourceID,1)));
    if strcmp(ext,'.JPG') || strcmp(ext,'.png')
        IM_Source = imread(char(StackList(SourceID,1)));
    else
        IM_Source = ImportStack(char(StackList(SourceID,1)));
    end
    
    if size(L,2) > 1
        tform = affine3d([[L(:,(3*SourceID)-2:3*SourceID)';b(:,SourceID)'],[0;0;0;1]]);
        IM_Source = imwarp(IM_Source,tform);
        %        figure();imshow(max(IM_Source,[],3),[0 max(IM_Source(:))]);
    end
    StackPosition = StackPositions_Registered(SourceID,:); %comming stack position
    StackPosition_org = StackPositions_pixels_org(SourceID,:); %comming stack position
    
    IM_Source_Pad = padarray(IM_Source,paramsBPextSize);
    %      figure();imshow(max(IM_Source_Pad,[],3),[0 max(IM_Source_Pad(:))]);
    translate = (StackPosition -S  ) ;
    translate_org = (StackPosition_org - S_org) ;
    
    %     translate = median(Matched{stackID,SourceID}(:,4:6)- Matched{stackID,SourceID}(:,1:3),1);
    %     translate (1) = translate(1) -1;
    %      translate (2) = translate(2) +1;
    %      translate (3) = translate(3) +1;
    IM_Source_Pad_translate = imtranslate(IM_Source_Pad,translate([2,1,3]));
    %            figure();imshow(max(IM_Source_Pad_translate,[],3),[0 max(IM_Source_Pad_translate(:))]);
    IM_Source_Pad_translate_org = imtranslate(IM_Source_Pad,translate_org([2,1,3]));
    
    if size(Tile3D,3) ~= size(IM_Source_Pad_translate,3) || size(Tile3D,2) ~= size(IM_Source_Pad_translate,2) || size(Tile3D,1) ~= size(IM_Source_Pad_translate,1)
        IM_Source_Pad_translate = imresize3(IM_Source_Pad_translate,[size(Tile3D,1) size(Tile3D,2) size(Tile3D,3)]) ;
    end
    if size(Tile3D_org,3) ~= size(IM_Source_Pad_translate_org,3) || size(Tile3D_org,2) ~= size(IM_Source_Pad_translate_org,2) || size(Tile3D_org,1) ~= size(IM_Source_Pad_translate_org,1)
        IM_Source_Pad_translate_org = imresize3(IM_Source_Pad_translate_org,[size(Tile3D_org,1) size(Tile3D_org,2) size(Tile3D_org,3)]) ;
    end
    
    Tile3D_org = max(Tile3D_org,IM_Source_Pad_translate_org);
    
    %     figure();imshow(max(Tile3D,[],3),[0 max(Tile3D(:))]);
    
    if paramsBPremove_pad && axes_xz.Value == 1
        [counts,x] = imhist(Tile3D,16);
        %     figure;stem(x,counts);
        [val, loc] = min(counts);
        
        thr = x(loc+2) + 20;
        Tile3D(find(Tile3D< thr))=0;
        IM_Source_Pad_translate(find(IM_Source_Pad_translate< thr))=0;
        IM_Source_Pad_translate(IM_Source_Pad_translate<thr) = 0;
        Tile3D(IM_Source_Pad_translate>0) =  IM_Source_Pad_translate(IM_Source_Pad_translate>0);
    else
        Tile3D = max(Tile3D,IM_Source_Pad_translate);
    end
    %               figure();imshow(max(Tile3D,[],3),[0 max(Tile3D(:))]);
    %      figure();imshow(max(IM_Source_Pad_translate,[],3),[0 max(IM_Source_Pad_translate(:))]);
end

if params.BP.saveImages
    saveastiff(Tile3D,[DataFolder,'\stacks_registered.tif']);
    saveastiff(Tile3D_org,[DataFolder,'\stacks_before_register.tif']);
end

end
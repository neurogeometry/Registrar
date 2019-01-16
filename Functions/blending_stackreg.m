function [Tile3D_org,Tile3D_reg] = blending_stackreg(StackPositions_pixels_original,StackSizes_pixels,StackList,stackID,L,b,DataFolder)
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


for i = 1:size(StackList,1)
    IM = imread(char(StackList(i,1)));
    Tile3D_org(:,:,i) = IM; 
    [IMs{i},StackPositions_Registered(:,i)]=Perform_Linear_Transform(IM,[1,1,1],[],b(:,i));
    Sizes(i,:)=[size(IMs{i}),1];
end
StackPositions_Registered=round(StackPositions_Registered');
Start=min(StackPositions_Registered);
End=max(StackPositions_Registered+Sizes-1);
Shift=1-Start;
Tile3D_reg=zeros(End-Shift,class(IM));
for i = 1:size(StackList,1)
    Tile3D_reg(StackPositions_Registered(i,1)+Shift(1):StackPositions_Registered(i,1)+Shift(1)+Sizes(i,1)-1,...
        StackPositions_Registered(i,2)+Shift(2):StackPositions_Registered(i,2)+Shift(2)+Sizes(i,2)-1,i+Shift(3))=IMs{i};
end




% Tile3D = zeros(StackSizes_pixels(1,1),StackSizes_pixels(1,2),size(StackSizes_pixels,1));
% Tile3D_org  = zeros(StackSizes_pixels(1,1),StackSizes_pixels(1,2),size(StackSizes_pixels,1));
% StackPositions_Registered_C = StackPositions_Registered(:,[2,1]);
% for i = 1:size(StackPositions_Registered,1)
% %     tb11 = findobj(NCT_Registration,'Tag', 'pushbutton10');
%     %         disp(num2str(tb11.UserData));
% %     if get(tb11,'userdata') || stop% stop condition
% %         disp(num2str(tb11.UserData));
% %         disp('Process Stopped');
% %         stop = 1;
% %         %                 listboxItems{v}  = 'Process Stopped ';
% %         %                 v = v + 1;
% %         %                 tb = findobj(NCT_Registration,'Tag', 'listbox1');
% %         %                 set(tb, 'String', listboxItems);drawnow
% %         %                 tb.Value = v-1;drawnow
% %         break;
% %     end
%     IM = imread(char(StackList(i,1)));
%     Tile3D_org(:,:,i) = IM;
%     IM = imtranslate(IM,StackPositions_Registered_C(i,:));
%     Tile3D(:,:,i) = IM;
% end


end
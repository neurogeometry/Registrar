close all;
clear all;
clc;
addpath ('NeuronTracerV20');
addpath ('../functions');




Displacement = [864.919440349651,0,0];



k = 1;
load('../data/StackData_Neuromuscular1.mat');
StackPositions_pixels_4Stacks(1,:) = StackPositions_pixels(106,:);
StackPositions_pixels_4Stacks(2,:) = StackPositions_pixels(107,:);
StackPositions_pixels_4Stacks(3,:) = StackPositions_pixels(114,:);
StackPositions_pixels_4Stacks(4,:) = StackPositions_pixels(115,:);
StackList_4Stacks(1,:) = StackList(106,:);
StackList_4Stacks(2,:) = StackList(107,:);
StackList_4Stacks(3,:) = StackList(114,:);
StackList_4Stacks(4,:) = StackList(115,:);
StackList = StackList_4Stacks;
showImage = 0;
showTraces = 0;
ppm=3;


StackPositions_pixels = StackPositions_pixels_4Stacks(:,[2,1,3]);

sourceID = 1;
targetID = 2;


Source_StackPositions = StackPositions_pixels(sourceID,:);
Target_StackPositions = StackPositions_pixels(targetID,:);

temp=StackList{sourceID,1};
sourcePath = ['../data/evaluation/DIADEM_Neuro/106-107/',temp,'_opt.swc'];%00760_Opt.swc';%
temp=StackList{targetID,1};
targetPath = ['../data/evaluation/DIADEM_Neuro/106-107/',temp,'_opt.swc'];%00759_Opt.swc';%



[AM_Source,r_Source,R_Source]=swc2AM(sourcePath);
[AM_Target,r_Target,R_Target]=swc2AM(targetPath);

[AM_Source,r_Source,R_Source] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
[AM_Target,r_Target,R_Target] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);



% Rigid
    FeaturePositions = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-Neuromuscular_StackList_4Stacks\MatchedPoints_Rigid.mat');
    second = FeaturePositions.Matched{sourceID,targetID}(:,1:3);
    first = FeaturePositions.Matched{sourceID,targetID}(:,4:6);
    
    [L_rigid,b_rigid]=Optimal_Rigid_Transform(first',second');
    r_Source_rigid = r_Source;%L*r_Source'+b'*ones(1,size(r_Source',2));
    r_Target_rigid = L_rigid*r_Target' + b_rigid*ones(1,size(r_Target',2));
%     r_Source_rigid = L_rigid*r_Source(:,[2,1,3])'+ b_rigid*ones(1,size(r_Source(:,[2,1,3])',2));
%     r_Target_rigid = r_Target(:,[2,1,3]);
    figure
        PlotAM(AM_Source,r_Source_rigid,'r')
        axis equal
        PlotAM(AM_Target,r_Target_rigid','g')
        axis equal
    [Distances_rigid,AllDistances_rigid] = TraceDistance(AM_Source, r_Source_rigid, AM_Target, r_Target_rigid',0);
    Rigid = mean(Distances_rigid)







    FeaturePositions = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-Neuromuscular_StackList_4Stacks\StackPositions_Registered.mat');
    
    NEW_StackPositions_pixels = FeaturePositions.StackPositions_Registered;
    r_Source_Global = r_Source + [NEW_StackPositions_pixels(sourceID,1),NEW_StackPositions_pixels(sourceID,2),NEW_StackPositions_pixels(sourceID,3)]-1;
    r_Target_Global = r_Target + [NEW_StackPositions_pixels(targetID,1),NEW_StackPositions_pixels(targetID,2),NEW_StackPositions_pixels(targetID,3)]-1;
    
%     Source_Stack_File = [char(StackList(sourceID,2)),'\',char(StackList(sourceID,1)),'.tif'];
%     Target_Stack_File = [char(StackList(targetID,2)),'\',char(StackList(targetID,1)),'.tif'];
%     IM_Source=ImportStack(char(Source_Stack_File));
%     IM_Target=ImportStack(char(Target_Stack_File));
%     IM_source_max=max(IM_Source,[],3);
%     IM_target_max=max(IM_Target,[],3);
%     
%     figure;imshow(IM_source_max,[0 max(IM_source_max(:))]);
%        hold on; PlotAM(AM_Source,r_Source,'r')
%     figure;imshow(IM_target_max,[0 max(IM_target_max(:))]);
%        hold on; PlotAM(AM_Target,r_Target,'r')
       
%     figure
%         PlotAM(AM_Source,r_Source_Global,'r')
%         axis equal
%         PlotAM(AM_Target,r_Target_Global,'g')
%         axis equal
    [Distances_Translation,AllDistances_Translation] = TraceDistance(AM_Source, r_Source_Global, AM_Target, r_Target_Global,1);
    mean(Distances_Translation)

    

     FeaturePositions = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-Neuromuscular_StackList_4Stacks\MatchedPoints_Affine.mat');
   % FeaturePositions = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList\MatchedPoints_Affine_ok.mat');
    
       
%        figure;
%        PlotAM(AM_Source,r_Source,'r')
%         axis equal
%         PlotAM(AM_Target,r_Target,'g')
%         axis equal
    [~,L,b]=Global_Optimal_Affine(StackPositions_pixels,FeaturePositions.Matched); %Matched_Global
    
    r_Source_Affine = L(:,(sourceID-1)*3+1:sourceID*3)*(r_Source+ones(size(r_Source,1),1)*Source_StackPositions)'+b(:,sourceID)*ones(1,size(r_Source',2));   
    r_Target_Affine = L(:,(targetID-1)*3+1:targetID*3)*(r_Target+ones(size(r_Target,1),1)*Target_StackPositions)'+b(:,targetID)*ones(1,size(r_Target',2));
    
%     figure;imshow(IM_source_max,[0 max(IM_source_max(:))]);
%        hold on; PlotAM(AM_Source,r_Source_Affine','r')
%     figure;imshow(IM_target_max,[0 max(IM_target_max(:))]);
%        hold on; PlotAM(AM_Target,r_Target_Affine','r')
%     figure
%         PlotAM(AM_Source,r_Source_Affine','r')
%         axis equal
%         PlotAM(AM_Target,r_Target_Affine','g')
%         axis equal
    
    % r_Target_Affine = L(:,(sourceID-1)*3+1:sourceID*3)*(r_Target(:,[2,1,3])+ones(size(r_Target,1),1)*Source_StackPositions)'+b(:,sourceID)*ones(1,size(r_Target(:,[2,1,3])',2));
   %  r_Source_Affine = L(:,(targetID-1)*3+1:targetID*3)*(r_Source(:,[2,1,3])+ones(size(r_Source,1),1)*Target_StackPositions)'+b(:,targetID)*ones(1,size(r_Source(:,[2,1,3])',2));

[Distances_Affine,AllDistances_Affine] = TraceDistance(AM_Source, r_Source_Affine', AM_Target, r_Target_Affine',0);
mean(Distances_Affine)


dx_GT = 902.8;
        dy_GT = 2.28;
        dz_GT = 4;
        
        
        
        % Fiji time = 270930ms - 68515ms --  Globally optimal stitching of tiled 3D microscopic image acquisitions
        dx_Fiji = 902;
        dy_Fiji = 2;
        dz_Fiji = 4;
        
        dx_XUV = 903;
        dy_XUV = -2;
        dz_XUV = 4; %  Z Ignored Z = 16.9375
        
        dx_Tara = 902;
            dy_Tara = 2;
            dz_Tara = 17;
        


OLD_StackPositions_pixels=StackPositions_pixels; % global stack positions in pixels


% N_L = 10;
% r_G_Nonrigid = r_G(:,[2,1,3])+Target_StackPositions;%L*r_G'+b'*ones(1,size(r_G',2));
% r_A_Nonrigid = r_A(1:554,[2,1,3])+Source_StackPositions;
% transformed = Perform_Nonrigid_Transform(r_A_Nonrigid',XYZlmn,N_L,Min,Max);
% [Distances_Nonrigid,AllDistances_Nonrigid] = TraceDistance( AM_G, r_G_Nonrigid,AM_A, transformed',showTraces);
%  mean(Distances_Nonrigid)


% N_L = 10;
% r_G_Nonrigid = r_G(:,[2,1,3])+Target_StackPositions;%L*r_G'+b'*ones(1,size(r_G',2));
% r_A_Nonrigid = r_A(:,[2,1,3])+Source_StackPositions;
% transformed = Perform_Nonrigid_Transform(r_G_Nonrigid',XYZlmn,N_L,Min,Max);




r_Source_old = r_Source + [OLD_StackPositions_pixels(sourceID,1),OLD_StackPositions_pixels(sourceID,2),OLD_StackPositions_pixels(sourceID,3)]-1;
r_Target_old = r_Target + [OLD_StackPositions_pixels(targetID,1),OLD_StackPositions_pixels(targetID,2),OLD_StackPositions_pixels(targetID,3)]-1;


r_Source_Fiji = r_Source;
r_Target_Fiji = r_Target +[dy_Fiji,dx_Fiji,dz_Fiji];

r_Source_Tara = r_Source;%
r_Target_Tara = r_Target +[dy_Tara,dx_Tara,dz_Tara];

r_Source_XUV = r_Source;
r_Target_XUV = r_Target +[dy_XUV,dx_XUV,dz_XUV];

r_Source_GT = r_Source;
r_Target_GT = r_Target +[dy_GT,dx_GT,dz_GT];





[Distances_Fiji,AllDistances_Fiji] = TraceDistance(AM_Source, r_Source_Fiji, AM_Target, r_Target_Fiji,showTraces);
[Distances_Tara,AllDistances_Tara] = TraceDistance(AM_Source, r_Source_Tara, AM_Target, r_Target_Tara,showTraces);
[Distances_XUV,AllDistances_XUV] = TraceDistance(AM_Source, r_Source_XUV, AM_Target, r_Target_XUV,showTraces);
[Distances_GT,AllDistances_GT] = TraceDistance(AM_Source, r_Source_GT, AM_Target, r_Target_GT,showTraces);
[Distances_old,AllDistances_old] = TraceDistance(AM_Source, r_Source_old, AM_Target, r_Target_old,showTraces);

M_Translation = mean(Distances_Translation)
M_Affine = mean(Distances_Affine)
M_Rigid = mean(Distances_rigid)
%     M_NonRigid = mean(Distances_Nonrigid)

M_Fiji = mean(Distances_Fiji)
M_Tera = mean(Distances_Tara)
M_XUV = mean(Distances_XUV)
M_GT = mean(Distances_GT)
M_Old = mean(Distances_old)

M_NonRigid = 0;
All = [M_Old,M_XUV,M_Tera,M_Fiji,M_Translation,M_Rigid,M_Affine,M_NonRigid]




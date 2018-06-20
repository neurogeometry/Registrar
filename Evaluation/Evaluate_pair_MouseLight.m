close all;
clear all;
clc;
addpath ('NeuronTracerV20');
addpath ('../functions');
addpath('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\Functions\');



Displacement = [864.919440349651,0,0];

load('../data/MatchedPoints');

k = 1;
load('../data/StackData.mat');
showImage = 0;
showTraces = 0;
ppm=3;

sourceID = 3;%
targetID = 4;%

temp=StackPositions_pixels;
StackPositions_pixels(:,1) = max(temp(:,2))-temp(:,2);
StackPositions_pixels(:,2) = max(temp(:,1))-temp(:,1);
StackPositions_pixels(:,3) = max(temp(:,3))-temp(:,3); % Not Sure - Need to check

Source_StackPositions = StackPositions_pixels(sourceID,:);
Target_StackPositions = StackPositions_pixels(targetID,:);

temp=StackList{sourceID,1};
sourcePath = ['../data/evaluation/3-759-760/',temp(1:find(temp=='-')-1),'_opt.swc'];%00760_Opt.swc';%
temp=StackList{targetID,1};
targetPath = ['../data/evaluation/3-759-760/',temp(1:find(temp=='-')-1),'_opt.swc'];%00759_Opt.swc';%

% sourcePath = '../data/evaluation/6-735-736/00736_opt.swc';%00760_Opt.swc';%
% targetPath = '../data/evaluation/6-735-736/00735_opt.swc';%00759_Opt.swc';%

[AM_Source,r_Source,R_Source]=swc2AM(sourcePath);
[AM_Target,r_Target,R_Target]=swc2AM(targetPath);

[AM_Source,r_Source,R_Source] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
[AM_Target,r_Target,R_Target] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);
    
   

FeaturePositions = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList\StackPositions_Registered.mat');
    
    NEW_StackPositions_pixels = FeaturePositions.StackPositions_Registered;
    r_Source_Global = r_Source + [NEW_StackPositions_pixels(sourceID,1),NEW_StackPositions_pixels(sourceID,2),NEW_StackPositions_pixels(sourceID,3)]-1;
    r_Target_Global = r_Target + [NEW_StackPositions_pixels(targetID,1),NEW_StackPositions_pixels(targetID,2),NEW_StackPositions_pixels(targetID,3)]-1;
    
    [Distances_Translation,AllDistances_Translation] = TraceDistance(AM_Source, r_Source_Global, AM_Target, r_Target_Global,0);
    Translation = mean(Distances_Translation)


     FeaturePositions = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList\MatchedPoints_Affine.mat');
   % FeaturePositions = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList\MatchedPoints_Affine_ok.mat');

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
Affine = mean(Distances_Affine)


 % Rigid
    FeaturePositions = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList\MatchedPoints_Rigid.mat');
    second = FeaturePositions.Matched{sourceID,targetID}(:,1:3);
    first = FeaturePositions.Matched{sourceID,targetID}(:,4:6);
    
    [L_rigid,b_rigid]=Optimal_Rigid_Transform(first',second');
    r_Source_rigid = r_Source;%L*r_Source'+b'*ones(1,size(r_Source',2));
    r_Target_rigid = L_rigid*r_Target' + b_rigid*ones(1,size(r_Target',2));
%     r_Source_rigid = L_rigid*r_Source(:,[2,1,3])'+ b_rigid*ones(1,size(r_Source(:,[2,1,3])',2));
%     r_Target_rigid = r_Target(:,[2,1,3]);
%     figure
%         PlotAM(AM_Source,r_Source_rigid,'r')
%         axis equal
%         PlotAM(AM_Target,r_Target_rigid','g')
%         axis equal
    [Distances_rigid,AllDistances_rigid] = TraceDistance(AM_Source, r_Source_rigid, AM_Target, r_Target_rigid',0);
    Rigid = mean(Distances_rigid)








dx_Fiji = -870;%pairwise
dy_Fiji = 9;
dz_Fiji = 1;


dx_Tara = -857;
dy_Tara = 9;
dz_Tara = 1;


dx_XUV = -432*2; %Pairwise
dy_XUV = 5*2;
dz_XUV = 1*2;

dx_GT = -865;
dy_GT = 9.8;
dz_GT = 1;


OLD_StackPositions_pixels=StackPositions_pixels; % global stack positions in pixels


% OLD_StackPositions_pixels(:,1) = max(OLD_StackPositions_pixels(:,1))-OLD_StackPositions_pixels(:,1);
% OLD_StackPositions_pixels(:,2) = max(OLD_StackPositions_pixels(:,2))-OLD_StackPositions_pixels(:,2);
% OLD_StackPositions_pixels(:,3) = max(OLD_StackPositions_pixels(:,2))-OLD_StackPositions_pixels(:,3);




% N_L = 10;
% r_Source_Nonrigid = r_Source(:,[2,1,3])+Target_StackPositions;%L*r_Source'+b'*ones(1,size(r_Source',2));
% r_Target_Nonrigid = r_Target(1:554,[2,1,3])+Source_StackPositions;
% transformed = Perform_Nonrigid_Transform(r_Target_Nonrigid',XYZlmn,N_L,Min,Max);
% [Distances_Nonrigid,AllDistances_Nonrigid] = TraceDistance( AM_Source, r_Source_Nonrigid,AM_Target, transformed',showTraces);
%  mean(Distances_Nonrigid)

%     N_L = 10;
% r_Source_Nonrigid = r_Source(:,[2,1,3])+Target_StackPositions;%L*r_Source'+b'*ones(1,size(r_Source',2));
% r_Target_Nonrigid = r_Target(:,[2,1,3])+Source_StackPositions;
% transformed = Perform_Nonrigid_Transform(r_Source_Nonrigid',XYZlmn,N_L,Min,Max);


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
% All = [All_Mean_Proposed_Translation,All_Mean_Proposed_Rigid,All_Mean_Proposed_Affine,All_Mean_Proposed_NonRigid,All_Mean_Proposed,All_Mean_Fiji,All_Mean_Tera,All_Mean_XUV,All_Mean_GT,All_Mean_Old,All_Mean_Affine,All_Mean_Projective]


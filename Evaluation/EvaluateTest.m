close all;
clear all;
clc;
addpath ('NeuronTracerV20');
addpath ('../functions');




Displacement = [864.919440349651,0,0];

load('../data/MatchedPoints');

k = 1;
load('../data/StackData.mat');
showImage = 0;
showTraces = 0;
ppm=3;

sourceID = 2;%
targetID = 1;%


StackPositions_pixels(:,1) = max(StackPositions_pixels(:,1))-StackPositions_pixels(:,1);
StackPositions_pixels(:,2) = max(StackPositions_pixels(:,2))-StackPositions_pixels(:,2);
StackPositions_pixels(:,3) = max(StackPositions_pixels(:,2))-StackPositions_pixels(:,3);

Source_StackPositions = StackPositions_pixels(sourceID,:);
Target_StackPositions = StackPositions_pixels(targetID,:);
sourcePath = '../data/evaluation/6-735-736/00736_opt.swc';%00760_Opt.swc';%
targetPath = '../data/evaluation/6-735-736/00735_opt.swc';%00759_Opt.swc';%

[AM_Source,r_Source,R_Source]=swc2AM(sourcePath);
[AM_Target,r_Target,R_Target]=swc2AM(targetPath);

[AM_Source,r_Source,R_Source] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
[AM_Target,r_Target,R_Target] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);
doglobal = 0;
if ~doglobal
    FeaturePositions = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList\MatchedPoints_Affine.mat');
    Matched = FeaturePositions.Matched;
    NewMatched=Matched;

%     for i = 1:size(Matched,1)
%         for j=1:size(Matched,2)
%             if i<=j && ~isempty(Matched{i,j})
%                 NewMatched{i,j}=[Matched{i,j}(:,4:6),Matched{i,j}(:,1:3)];
%             end
%         end
%     end
    
    
    
    [X_aligned,L,b]=Global_Optimal_Affine(StackPositions_pixels,NewMatched); %Matched_Global
    
  %  r_Target_Affine = L(:,(targetID-1)*3+1:targetID*3)*(r_Target(:,[2,1,3])+ones(size(r_Target,1),1)*Source_StackPositions)'+b(:,targetID)*ones(1,size(r_Target(:,[2,1,3])',2));
  %  r_Source_Affine = L(:,(sourceID-1)*3+1:sourceID*3)*(r_Source(:,[2,1,3])+ones(size(r_Source,1),1)*Target_StackPositions)'+b(:,sourceID)*ones(1,size(r_Source(:,[2,1,3])',2));
    r_Target_Affine = L(:,(sourceID-1)*3+1:sourceID*3)*(r_Target(:,[2,1,3])+ones(size(r_Target,1),1)*Source_StackPositions)'+b(:,sourceID)*ones(1,size(r_Target(:,[2,1,3])',2));
     r_Source_Affine = L(:,(targetID-1)*3+1:targetID*3)*(r_Source(:,[2,1,3])+ones(size(r_Source,1),1)*Target_StackPositions)'+b(:,targetID)*ones(1,size(r_Source(:,[2,1,3])',2));
else
    FeaturePositions = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList\MatchedPoints_Affine_Global.mat');
    [X_aligned,L,b]=Global_Optimal_Affine(StackPositions_pixels,FeaturePositions.Matched_Global); %Matched_Global
     r_Target_Affine = L(:,(sourceID-1)*3+1:sourceID*3)*r_Target(:,[2,1,3])'+b(:,sourceID)*ones(1,size(r_Target(:,[2,1,3])',2));
    r_Source_Affine = L(:,(targetID-1)*3+1:targetID*3)*r_Source(:,[2,1,3])'+b(:,targetID)*ones(1,size(r_Source(:,[2,1,3])',2));
end

[Distances_Affine,AllDistances_Affine] = TraceDistance(AM_Source, r_Source_Affine', AM_Target, r_Target_Affine',showTraces);
mean(Distances_Affine)


b_rigid = [880.83598632946336692839, -10.49563173218984957202, -0.27151260280197675456];
L_rigid = [0.99991204079478157585, -0.01244598255129452345, -0.00458346942262180014;0.01244988197625845720, 0.99992215835216591291, 0.00082321025477287142;0.00457286697734222199, -0.00088020149920634450, 0.99998915700767887493];



% 3  00760 00759
P1_P2_DX = 864.97;% 864.756756756757;
P1_P2_DY = -9.34;%-9.37837837837838;%-10
P1_P2_DZ = 1.12;%1.05405405405405;%-10

dx_Fiji = 870;%pairwise
dy_Fiji = -9;
dz_Fiji = 1;


dx_Tara = 857;
dy_Tara = -9;
dz_Tara = 1;


dx_XUV = 432*2; %Pairwise
dy_XUV = 5*2;
dz_XUV = 1*2;

dx_GT = 865;
dy_GT = -9.8;
dz_GT = 1;


OLD_StackPositions_pixels=StackPositions_pixels; % global stack positions in pixels


OLD_StackPositions_pixels(:,1) = max(OLD_StackPositions_pixels(:,1))-OLD_StackPositions_pixels(:,1);
OLD_StackPositions_pixels(:,2) = max(OLD_StackPositions_pixels(:,2))-OLD_StackPositions_pixels(:,2);
OLD_StackPositions_pixels(:,3) = max(OLD_StackPositions_pixels(:,2))-OLD_StackPositions_pixels(:,3);





r_Source_rigid = r_Source(:,[2,1,3]);%L*r_Source'+b'*ones(1,size(r_Source',2));
r_Target_rigid = L_rigid*r_Target(:,[2,1,3])' + b_rigid'*ones(1,size(r_Target(:,[2,1,3])',2));


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


r_Source_old = r_Source + [OLD_StackPositions_pixels(sourceID,2),OLD_StackPositions_pixels(sourceID,1),OLD_StackPositions_pixels(sourceID,3)]-1;
r_Target_old = r_Target + [OLD_StackPositions_pixels(targetID,2),OLD_StackPositions_pixels(targetID,1),OLD_StackPositions_pixels(targetID,3)]-1;


r_Source_Matching = r_Source;
r_Target_Matching = r_Target +[P1_P2_DY,P1_P2_DX,P1_P2_DZ];

r_Source_Fiji = r_Source;
r_Target_Fiji = r_Target +[dy_Fiji,dx_Fiji,dz_Fiji];

r_Source_Tara = r_Source;%
r_Target_Tara = r_Target +[dy_Tara,dx_Tara,dz_Tara];

r_Source_XUV = r_Source;
r_Target_XUV = r_Target +[dy_XUV,dx_XUV,dz_XUV];

r_Source_GT = r_Source;
r_Target_GT = r_Target +[dy_GT,dx_GT,dz_GT];



[Distances_rigid,AllDistances_rigid] = TraceDistance(AM_Source, r_Source_rigid, AM_Target, r_Target_rigid',showTraces);
[Distances_Proposed, AllDistances_Proposed] = TraceDistance(AM_Source, r_Source_Matching, AM_Target, r_Target_Matching,showTraces);
[Distances_Fiji,AllDistances_Fiji] = TraceDistance(AM_Source, r_Source_Fiji, AM_Target, r_Target_Fiji,showTraces);
[Distances_Tara,AllDistances_Tara] = TraceDistance(AM_Source, r_Source_Tara, AM_Target, r_Target_Tara,showTraces);
[Distances_XUV,AllDistances_XUV] = TraceDistance(AM_Source, r_Source_XUV, AM_Target, r_Target_XUV,showTraces);
[Distances_GT,AllDistances_GT] = TraceDistance(AM_Source, r_Source_GT, AM_Target, r_Target_GT,showTraces);
[Distances_old,AllDistances_old] = TraceDistance(AM_Source, r_Source_old, AM_Target, r_Target_old,showTraces);

M_Affine = mean(Distances_Affine)
M_Rigid = mean(Distances_rigid)
%     M_NonRigid = mean(Distances_Nonrigid)
M_Translation = mean(Distances_Proposed)
M_Fiji = mean(Distances_Fiji)
M_Tera = mean(Distances_Tara)
M_XUV = mean(Distances_XUV)
M_GT = mean(Distances_GT)
M_Old = mean(Distances_old)

M_NonRigid = 0;
All = [M_Affine,M_Rigid,M_NonRigid,M_Translation,M_Fiji,M_Tera,M_XUV,M_Old]
% All = [All_Mean_Proposed_Translation,All_Mean_Proposed_Rigid,All_Mean_Proposed_Affine,All_Mean_Proposed_NonRigid,All_Mean_Proposed,All_Mean_Fiji,All_Mean_Tera,All_Mean_XUV,All_Mean_GT,All_Mean_Old,All_Mean_Affine,All_Mean_Projective]


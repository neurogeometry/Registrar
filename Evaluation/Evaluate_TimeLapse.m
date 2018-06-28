% clear all;
% clc
% close all;
ppm = 2;
affine = 1;
addpath('NeuronTracerV20');
% GTpath = 'E:\Datasets\TimeLaps\';
% GTpath = E:\Shih-Luen\Lab\Projects\RegistrationEvaluation\TimeLaps\';
% Functionspath = 'E:\Shih-Luen\Lab\Projects\registrar\Functions\';
% resultpath = 'C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-TimeLapse_Holtmaat_StackList\';


GTpath = '..\..\RegistrationEvaluation\TimeLaps\';
Functionspath = '..\';
resultpath = '..\..\RegistrationEvaluation\Results-TimeLapse_Holtmaat_StackList\';

addpath([Functionspath,'Functions']);
% load('../data/StackData_TimeLapse_Holtmaat.mat');
load([GTpath,'Matches\DL083-001-Matches.mat']);
showTranceonImage = 0;

T_Names = {'B','C','D','E','F','G','H','I','J','K','L','M','N'};

% Show Traces
% sourceID = 2;
% targetID = 3;
tic
for sourceID = 1:12
%     disp(sourceID)
    targetID = sourceID + 1;
    FeaturePositions_NR = load([resultpath,'MatchedPoints_Non-Rigid.mat']);
    
    % B-Spline-NonRigid
    Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
    Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
    Minimum = min(Global_Matched_Source,[],2);
    Maximum = max(Global_Matched_Source,[],2);
    nxyz = [512;512;156]; %image size
%     Nxyz = ceil((Maximum-Minimum)./nxyz');
    [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine);
    result{sourceID}.T.L = L;
    result{sourceID}.T.b = b;
    result{sourceID}.T.Cxyz = Cxyz;
    result{sourceID}.T.Nxyz = Nxyz;
    result{sourceID}.T.nxyz = nxyz;
    result{sourceID}.T.Grid_start = Grid_start;
    result{sourceID}.T.affine = affine;


fname_First = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
fname_First={fname_First.name}';
fname_Second = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
fname_Second={fname_Second.name}';
AllDistances = zeros(size(fname_First,1),6);
pixelSize = [0.26 0.26 0.8];

        %Use Boutons
        Boutons = struct2cell(Matches);
        B1 = Boutons{sourceID,1};
        SourcePoints = B1.r2;
        TargetPoints = B1.r1; 
%         D_NonRigid_voxel = mean(mean((SourcePoints_NR'-TargetPoints).^2,1).^0.5);
%         D_NonRigid_voxel = SourcePoints_NR'-TargetPoints;
[SourcePoints_NR,~]=Perform_Bspline_Transform(SourcePoints',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
    result{sourceID}.Bouton.r1 = SourcePoints;
    result{sourceID}.Bouton.r1_NR = SourcePoints_NR';
    result{sourceID}.Bouton.r2 = TargetPoints;

for i=1:size(fname_First,1)
        % Using Trace
        sourcePath = [GTpath,'Matches\Traces\',fname_Second{i}];%'E:\Datasets\TimeLaps\Matches\Traces\DL083B001-A001.swc';
        targetPath = [GTpath,'Matches\Traces\',fname_First{i}];%'E:\Datasets\TimeLaps\Matches\Traces\DL083C001-A001.swc';
        %         sourcePath = 'E:\Datasets\TimeLaps\Matches\Traces\DL083B001-A001.swc';
        %         targetPath = 'E:\Datasets\TimeLaps\Matches\Traces\DL083C001-A001.swc';
        
        [AM_Source,r_Source,R_Source]=swc2AM(sourcePath);
        [AM_Target,r_Target,R_Target]=swc2AM(targetPath);
        [AM_Source,r_Source,~] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
        [AM_Target,r_Target,~] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);
        SourcePoints = r_Source;
        TargetPoints = r_Target;

    
    
    
%     sourceID = 1;
%     targetID = 2;
    
    % Source_StackPositions = StackPositions_pixels(sourceID,:);
    % Target_StackPositions = StackPositions_pixels(targetID,:);
    
%     if showTranceonImage
%         Source_Stack_File = [char(StackList(sourceID,2)),'\',char(StackList(sourceID,1)),'.tif'];
%         Target_Stack_File = [char(StackList(targetID,2)),'\',char(StackList(targetID,1)),'.tif'];
%         IM_Source=ImportStack(char(Source_Stack_File));
%         IM_Target=ImportStack(char(Target_Stack_File));
%         IM_source_max=max(IM_Source,[],3);
%         IM_target_max=max(IM_Target,[],3);
%         if useTrace
%             figure;imshow(IM_source_max,[0 500]);
%             hold on; PlotAM(AM_Source,r_Source,'r')
%             figure;imshow(IM_target_max,[0 500]);
%             hold on; PlotAM(AM_Target,r_Target,'g')
%         else
%             figure;imshow(IM_source_max,[0 500]);
%             hold on; plot(SourcePoints(:,2),SourcePoints(:,1),'*r')
%             figure;imshow(IM_target_max,[0 500]);
%             hold on; plot(TargetPoints(:,2),TargetPoints(:,1),'*r')
%         end
%     end
    %         figure
    %         PlotAM(AM_Source,SourcePoints,'r')
    %         axis equal
    %         PlotAM(AM_Target,TargetPoints,'g')
    %         axis equal
%     if useTrace
%         [Dis_Before_um,Dis_Before_voxel] = TraceDistance(AM_Source, SourcePoints, AM_Target, TargetPoints,pixelSize,0);
%         D_Before_um = mean(Dis_Before_um);
%         D_Before_voxel = mean(Dis_Before_voxel);
%     else
% %         D_Before_voxel = mean(mean((SourcePoints-TargetPoints).^2,1).^0.5);
%         D_Before_voxel = SourcePoints-TargetPoints;
% 
%     end
   
%     [SourcePoints_NR,~]=Perform_Bspline_Transform(SourcePoints',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);

%         [Dis_NonRigid_um,Dis_NonRigid_voxel] = TraceDistance(AM_Source, SourcePoints_NR', AM_Target, TargetPoints,pixelSize,0);
%         D_NonRigid_um = mean(Dis_NonRigid_um);
%         D_NonRigid_voxel = mean(Dis_NonRigid_voxel);
    [SourcePoints_NR,~]=Perform_Bspline_Transform(SourcePoints',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
    result{sourceID}.Trace.r1{i} = SourcePoints;
    result{sourceID}.Trace.r1_NR{i} = SourcePoints_NR';
    result{sourceID}.Trace.r2{i} = TargetPoints;
    result{sourceID}.Trace.AM1{i} = AM_Source;
    result{sourceID}.Trace.AM2{i} = AM_Target;
    result{sourceID}.Trace.R1{i} = R_Source;
    result{sourceID}.Trace.R2{i} = R_Target;

%     if ~useTrace
% 
%     end
%     disp('-------------------------------------');
%     
%     D_Before_voxel
%     D_NonRigid_voxel
%     
%     if useTrace
%     AllDistances_um(i,1) = D_Before_um;
%     AllDistances_um(i,6) = D_NonRigid_um;
%     AllDistances_voxel(i,1) = D_Before_voxel;
%     AllDistances_voxel(i,6) = D_NonRigid_voxel;
%     else
%       AllDistances_bouton = D_Before_voxel;
%     AllDistances_bouton_after = D_NonRigid_voxel;
%     end
%     if ~useTrace
%         break;
%     end
end
end

% save(['E:\Shih-Luen\Lab\Projects\RegistrationEvaluation\test',T_Names{sourceID},'and',T_Names{targetID},'.mat'],'A')
save(['..\..\RegistrationEvaluation\result_T.mat'],'result','pixelSize')
toc
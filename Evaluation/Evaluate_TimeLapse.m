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


GTpath = 'E:\Shih-Luen\Lab\Projects\RegistrationEvaluation\TimeLaps\';
Functionspath = 'E:\Shih-Luen\Lab\Projects\registrar\';
resultpath = 'E:\Shih-Luen\Lab\Projects\RegistrationEvaluation\Results-TimeLapse_Holtmaat_StackList\';

addpath([Functionspath,'Functions']);
% load('../data/StackData_TimeLapse_Holtmaat.mat');
load([GTpath,'Matches\DL083-001-Matches.mat']);
useTrace = 1; % 0 = use boutons 1 = use trace
showTranceonImage = 0;

T_Names = {'B','C','D','E','F','G','H','I','J','K','L','M','N'};

% Show Traces
% sourceID = 2;
% targetID = 3;
for sourceID = 1:100
%     disp(sourceID)
    targetID = sourceID + 1;
% for useTrace = 0:1

fname_First = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
fname_First={fname_First.name}';
fname_Second = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
fname_Second={fname_Second.name}';
AllDistances = zeros(size(fname_First,1),6);
pixelSize = [0.26 0.26 0.8];
%     Source_Stack_File = [char(StackList(sourceID,2)),'\',char(StackList(sourceID,1)),'.tif'];
%         Target_Stack_File = [char(StackList(targetID,2)),'\',char(StackList(targetID,1)),'.tif'];
%         IM_Source=ImportStack(char(Source_Stack_File));
%         IM_Target=ImportStack(char(Target_Stack_File));
%         IM_source_max=max(IM_Source,[],3);
%         IM_target_max=max(IM_Target,[],3);
%
%          figure(1);imshow(IM_source_max,[0 500]);
%
%             figure(2);imshow(IM_target_max,[0 500]);
%
% for i=1:size(fname_B,1)
%     sourcePath = ['E:\Datasets\TimeLaps\Matches\Traces\',fname_B{i}];%'E:\Datasets\TimeLaps\Matches\Traces\DL083B001-A001.swc';
%     targetPath = ['E:\Datasets\TimeLaps\Matches\Traces\',fname_C{i}];%'E:\Datasets\TimeLaps\Matches\Traces\DL083C001-A001.swc';
%     [AM_Source,r_Source,R_Source]=swc2AM(sourcePath);
%         [AM_Target,r_Target,R_Target]=swc2AM(targetPath);
% %         [AM_Source,r_Source,~] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
% %         [AM_Target,r_Target,~] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);
%         SourcePoints = r_Source;
%         TargetPoints = r_Target;
%         Tracecolor = rand(1,3);
%     figure(1);hold on; PlotAM(AM_Source,r_Source, Tracecolor)
%     figure(2);hold on; PlotAM(AM_Target,r_Target, Tracecolor)
% end



for i=1:size(fname_First,1)
    if useTrace
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
    else
        %Use Boutons
        Boutons = struct2cell(Matches);
        B1 = Boutons{sourceID,1};
        SourcePoints = B1.r2;
        TargetPoints = B1.r1; 
    end
    
    
    
%     sourceID = 1;
%     targetID = 2;
    
    % Source_StackPositions = StackPositions_pixels(sourceID,:);
    % Target_StackPositions = StackPositions_pixels(targetID,:);
    
    if showTranceonImage
        Source_Stack_File = [char(StackList(sourceID,2)),'\',char(StackList(sourceID,1)),'.tif'];
        Target_Stack_File = [char(StackList(targetID,2)),'\',char(StackList(targetID,1)),'.tif'];
        IM_Source=ImportStack(char(Source_Stack_File));
        IM_Target=ImportStack(char(Target_Stack_File));
        IM_source_max=max(IM_Source,[],3);
        IM_target_max=max(IM_Target,[],3);
        if useTrace
            figure;imshow(IM_source_max,[0 500]);
            hold on; PlotAM(AM_Source,r_Source,'r')
            figure;imshow(IM_target_max,[0 500]);
            hold on; PlotAM(AM_Target,r_Target,'g')
        else
            figure;imshow(IM_source_max,[0 500]);
            hold on; plot(SourcePoints(:,2),SourcePoints(:,1),'*r')
            figure;imshow(IM_target_max,[0 500]);
            hold on; plot(TargetPoints(:,2),TargetPoints(:,1),'*r')
        end
    end
    %         figure
    %         PlotAM(AM_Source,SourcePoints,'r')
    %         axis equal
    %         PlotAM(AM_Target,TargetPoints,'g')
    %         axis equal
    if useTrace
        [Dis_Before_um,Dis_Before_voxel] = TraceDistance(AM_Source, SourcePoints, AM_Target, TargetPoints,pixelSize,0);
        D_Before_um = mean(Dis_Before_um);
        D_Before_voxel = mean(Dis_Before_voxel);
    else
%         D_Before_voxel = mean(mean((SourcePoints-TargetPoints).^2,1).^0.5);
        D_Before_voxel = SourcePoints-TargetPoints;

    end
    
    %     Dag_temp=(bsxfun(@minus,SourcePoints(:,1),TargetPoints(:,1)').^2+bsxfun(@minus,SourcePoints(:,2),TargetPoints(:,2)').^2+bsxfun(@minus,SourcePoints(:,3),TargetPoints(:,3)').^2).^0.5;
    %     [AM,~]=Hungarian_fast(Dag_temp);
    %     [D_Before,~] = TraceDistance(AM, SourcePoints, AM, TargetPoints,0)
    
    
    %     FeaturePositions_T = load([resultpath,'MatchedPoints_Translation.mat']);
    %     FeaturePositions_R = load([resultpath,'MatchedPoints_Rigid.mat']);
    %     FeaturePositions_A = load([resultpath,'MatchedPoints_Affine.mat']);
    FeaturePositions_NR = load([resultpath,'MatchedPoints_Non-Rigid.mat']);
    
    
    
    %     % Trnslation
    %     Global_Matched_Source = FeaturePositions_T.Matched{sourceID,targetID}(:,1:3)';
    %     Global_Matched_Target = FeaturePositions_T.Matched{sourceID,targetID}(:,4:6)';
    %     b=Optimal_Translation_Transform(Global_Matched_Source,Global_Matched_Target);
    %
    %     SourcePoints_TranslationPairs = SourcePoints+b';
    %     TargetPoints_TranslationPairs = TargetPoints;
    %
    %     if useTrace
    %         [Distances_Translation,~] = TraceDistance(AM_Source, SourcePoints_TranslationPairs, AM_Target, TargetPoints_TranslationPairs,pixelSize,0);
    %         D_Translation = mean(Distances_Translation)
    %     else
    %         D_Translation = mean(mean((SourcePoints_TranslationPairs-TargetPoints_TranslationPairs).^2,1).^0.5)
    %     end
    %
    %     % Rigid
    %     Global_Matched_Source = FeaturePositions_R.Matched{sourceID,targetID}(:,1:3)';
    %     Global_Matched_Target = FeaturePositions_R.Matched{sourceID,targetID}(:,4:6)';
    %     [R,b]=Optimal_Rigid_Transform(Global_Matched_Source,Global_Matched_Target);
    %
    %     SourcePoints_Rigid = (R*(SourcePoints)'+b)';
    %     TargetPoints_Rigid = TargetPoints;
    %     if useTrace
    %         [Distances_Rigid,~] = TraceDistance(AM_Source, SourcePoints_Rigid, AM_Target, TargetPoints_Rigid,pixelSize,0);
    %         D_Rigid = mean(Distances_Rigid)
    %     else
    %         D_Rigid = mean(mean((SourcePoints_Rigid-TargetPoints_Rigid).^2,1).^0.5)
    %     end
    %
    %     % Affine
    %     Global_Matched_Source = FeaturePositions_A.Matched{sourceID,targetID}(:,1:3)';
    %     Global_Matched_Target = FeaturePositions_A.Matched{sourceID,targetID}(:,4:6)';
    %     [L,b]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target);
    %
    %     SourcePoints_Affine = (L*(SourcePoints)'+b)';
    %     TargetPoints_Affine = TargetPoints;
    %     if useTrace
    %         [Distances_Affine,~] = TraceDistance(AM_Source, SourcePoints_Affine, AM_Target, TargetPoints_Affine,pixelSize,0);
    %         D_Affine = mean(Distances_Affine)
    %     else
    %         D_Affine = mean(mean((SourcePoints_Affine-TargetPoints_Affine).^2,1).^0.5)
    %     end
    %
    %
    %
    %     % AffineNonRigid
    %     Global_Matched_Source = FeaturePositions_A.Matched{sourceID,targetID}(:,1:3)';
    %     Global_Matched_Target = FeaturePositions_A.Matched{sourceID,targetID}(:,4:6)';
    %     [L,b]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target);
    %
    %     SourcePoints_Affine = (L*(SourcePoints)'+b)';
    %     TargetPoints_Affine = TargetPoints;
    %
    %     N_L=3;
    %     temp = (L*Global_Matched_Source+b);
    %     Min=min([min(temp,[],2),min(Global_Matched_Target,[],2)],[],2);
    %     Max=max([max(temp,[],2),max(Global_Matched_Target,[],2)],[],2);
    %     [~,XYZlmn,N_L,Min,Max]=Optimal_Nonrigid_Transform(temp,Global_Matched_Target,N_L,Min,Max);
    %     SourcePoints_Affine=Perform_Nonrigid_Transform(SourcePoints_Affine',XYZlmn,N_L,Min,Max)';
    %     if useTrace
    %         [Distances_NonRigid,~] = TraceDistance(AM_Source, SourcePoints_Affine, AM_Target, TargetPoints_Affine,pixelSize,0);
    %         D_NonRigidAffine = mean(Distances_NonRigid)
    %     else
    %         D_NonRigidAffine = mean(mean((SourcePoints_Affine-TargetPoints_Affine).^2,1).^0.5)
    %     end
    
    
    
    %     % NonRigid
    %     Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
    %     Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
    %     N_L=3;
    %     Min=min([min(Global_Matched_Source,[],2),min(Global_Matched_Target,[],2)],[],2);
    %     Max=max([max(Global_Matched_Source,[],2),max(Global_Matched_Target,[],2)],[],2);
    %     [~,XYZlmn,N_L,Min,Max]=Optimal_Nonrigid_Transform(Global_Matched_Source,Global_Matched_Target,N_L,Min,Max);
    %     SourcePoints_NR=Perform_Nonrigid_Transform(SourcePoints',XYZlmn,N_L,Min,Max);
    %     if useTrace
    %         [Distances_NonRigid,~] = TraceDistance(AM_Source, SourcePoints_NR', AM_Target, TargetPoints,pixelSize,0);
    %         D_NonRigid = mean(Distances_NonRigid);
    %     else
    %         D_NonRigid = mean(mean((SourcePoints_NR-Global_Matched_Target).^2,1).^0.5)
    %     end
    %     disp('-------------------------------------');
    
    % B-Spline-NonRigid
    Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
    Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
%     N_L=3;
%     Min=min([min(Global_Matched_Source,[],2),min(Global_Matched_Target,[],2)],[],2);
%     Max=max([max(Global_Matched_Source,[],2),max(Global_Matched_Target,[],2)],[],2);
%     [~,XYZlmn,N_L,Min,Max]=Optimal_Nonrigid_Transform(Global_Matched_Source,Global_Matched_Target,N_L,Min,Max);
%     SourcePoints_NR=Perform_Nonrigid_Transform(SourcePoints',XYZlmn,N_L,Min,Max);
    Minimum = min(Global_Matched_Source,[],2);
    Maximum = max(Global_Matched_Source,[],2);
    nxyz = [512;512;156]; %image size
%     Nxyz = ceil((Maximum-Minimum)./nxyz');
    [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine);
    [SourcePoints_NR,~]=Perform_Bspline_Transform(SourcePoints',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
    if useTrace
        [Dis_NonRigid_um,Dis_NonRigid_voxel] = TraceDistance(AM_Source, SourcePoints_NR', AM_Target, TargetPoints,pixelSize,0);
        D_NonRigid_um = mean(Dis_NonRigid_um);
        D_NonRigid_voxel = mean(Dis_NonRigid_voxel);
    else
%         D_NonRigid_voxel = mean(mean((SourcePoints_NR'-TargetPoints).^2,1).^0.5);
        D_NonRigid_voxel = SourcePoints_NR'-TargetPoints;

    end
    disp('-------------------------------------');
    
    D_Before_voxel
    D_NonRigid_voxel
    
%     D_Before_um
%     D_NonRigid_um
    
    if useTrace
    AllDistances_um(i,1) = D_Before_um;
%     AllDistances_um(i,2) = D_Translation_um;
%     AllDistances_um(i,3) = D_Rigid_um;
%     AllDistances_um(i,4) = D_Affine_um;
%     AllDistances_um(i,5) = D_NonRigidAffine_um;
    AllDistances_um(i,6) = D_NonRigid_um;
    AllDistances_voxel(i,1) = D_Before_voxel;
%     AllDistances_um(i,2) = D_Translation_um;
%     AllDistances_um(i,3) = D_Rigid_um;
%     AllDistances_um(i,4) = D_Affine_um;
%     AllDistances_um(i,5) = D_NonRigidAffine_um;
    AllDistances_voxel(i,6) = D_NonRigid_voxel;
    else
      AllDistances_bouton = D_Before_voxel;
%     AllDistances_um(i,2) = D_Translation_um;
%     AllDistances_um(i,3) = D_Rigid_um;
%     AllDistances_um(i,4) = D_Affine_um;
%     AllDistances_um(i,5) = D_NonRigidAffine_um;
    AllDistances_bouton_after = D_NonRigid_voxel;
    end
    if ~useTrace
        break;
    end
end

end
save(['E:\Shih-Luen\Lab\Projects\RegistrationEvaluation\test',num2str(sourceID),'.mat'],'AllDistances_um','AllDistances_voxel','AllDistances_bouton_before')
% end
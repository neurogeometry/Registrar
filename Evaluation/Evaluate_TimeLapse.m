clear all;
clc
close all;
ppm = 2;
addpath('NeuronTracerV20');
addpath('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\Functions');
% load('../data/StackData_TimeLapse_Holtmaat.mat');
load('E:\Datasets\TimeLaps\Matches\DL083-001-Matches.mat');
useTrace = 1;
showTranceonImage = 0;

fname_First = dir('E:\Datasets\TimeLaps\Matches\Traces\DL083C001-A0*');
fname_First={fname_First.name}';
fname_Second = dir('E:\Datasets\TimeLaps\Matches\Traces\DL083D001-A0*');
fname_Second={fname_Second.name}';
AllDistances = zeros(size(fname_First,1),6);
pixelSize = [0.26 0.26 0.8];

% Show Traces
sourceID = 1;
targetID = 2;


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
        sourcePath = ['E:\Datasets\TimeLaps\Matches\Traces\',fname_First{i}];%'E:\Datasets\TimeLaps\Matches\Traces\DL083B001-A001.swc';
        targetPath = ['E:\Datasets\TimeLaps\Matches\Traces\',fname_Second{i}];%'E:\Datasets\TimeLaps\Matches\Traces\DL083C001-A001.swc';
        %         sourcePath = 'E:\Datasets\TimeLaps\Matches\Traces\DL083B001-A001.swc';
        %         targetPath = 'E:\Datasets\TimeLaps\Matches\Traces\DL083C001-A001.swc';
        
        [AM_Source,r_Source,R_Source]=swc2AM(sourcePath);
        [AM_Target,r_Target,R_Target]=swc2AM(targetPath);
        [AM_Source,r_Source,~] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
        [AM_Target,r_Target,~] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);
        SourcePoints = r_Source;
        TargetPoints = r_Target;
    else
        SourcePoints = Matches.BC.r1;
        TargetPoints = Matches.BC.r2;
    end
    
    
    
    sourceID = 1;
    targetID = 2;
    
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
            [Distances_Before,~] = TraceDistance(AM_Source, SourcePoints, AM_Target, TargetPoints,pixelSize,0);
            D_Before = mean(Distances_Before)
        else
            D_Before = mean(mean((SourcePoints-TargetPoints).^2,1).^0.5)
        end
    
    %     Dag_temp=(bsxfun(@minus,SourcePoints(:,1),TargetPoints(:,1)').^2+bsxfun(@minus,SourcePoints(:,2),TargetPoints(:,2)').^2+bsxfun(@minus,SourcePoints(:,3),TargetPoints(:,3)').^2).^0.5;
    %     [AM,~]=Hungarian_fast(Dag_temp);
    %     [D_Before,~] = TraceDistance(AM, SourcePoints, AM, TargetPoints,0)
    
    
    FeaturePositions_T = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-TimeLapse_Holtmaat_StackList\MatchedPoints_Translation.mat');
    FeaturePositions_R = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-TimeLapse_Holtmaat_StackList\MatchedPoints_Rigid.mat');
    FeaturePositions_A = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-TimeLapse_Holtmaat_StackList\MatchedPoints_Affine.mat');
    FeaturePositions_NR = load('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\MicroscopeFiles\Results-TimeLapse_Holtmaat_StackList\MatchedPoints_Non-Rigid.mat');
    
    
    
    % Trnslation
    Global_Matched_Source = FeaturePositions_T.Matched{sourceID,targetID}(:,1:3)';
    Global_Matched_Target = FeaturePositions_T.Matched{sourceID,targetID}(:,4:6)';
    b=Optimal_Translation_Transform(Global_Matched_Source,Global_Matched_Target);
    
    SourcePoints_TranslationPairs = SourcePoints+b';
    TargetPoints_TranslationPairs = TargetPoints;
    
    if useTrace
        [Distances_Translation,~] = TraceDistance(AM_Source, SourcePoints_TranslationPairs, AM_Target, TargetPoints_TranslationPairs,pixelSize,0);
        D_Translation = mean(Distances_Translation)
    else
        D_Translation = mean(mean((SourcePoints_TranslationPairs-TargetPoints_TranslationPairs).^2,1).^0.5)
    end
    
    % Rigid
    Global_Matched_Source = FeaturePositions_R.Matched{sourceID,targetID}(:,1:3)';
    Global_Matched_Target = FeaturePositions_R.Matched{sourceID,targetID}(:,4:6)';
    [R,b]=Optimal_Rigid_Transform(Global_Matched_Source,Global_Matched_Target);
    
    SourcePoints_Rigid = (R*(SourcePoints)'+b)';
    TargetPoints_Rigid = TargetPoints;
    if useTrace
        [Distances_Rigid,~] = TraceDistance(AM_Source, SourcePoints_Rigid, AM_Target, TargetPoints_Rigid,pixelSize,0);
        D_Rigid = mean(Distances_Rigid)
    else
        D_Rigid = mean(mean((SourcePoints_Rigid-TargetPoints_Rigid).^2,1).^0.5)
    end
    
    % Affine
    Global_Matched_Source = FeaturePositions_A.Matched{sourceID,targetID}(:,1:3)';
    Global_Matched_Target = FeaturePositions_A.Matched{sourceID,targetID}(:,4:6)';
    [L,b]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target);
    
    SourcePoints_Affine = (L*(SourcePoints)'+b)';
    TargetPoints_Affine = TargetPoints;
    if useTrace
        [Distances_Affine,~] = TraceDistance(AM_Source, SourcePoints_Affine, AM_Target, TargetPoints_Affine,pixelSize,0);
        D_Affine = mean(Distances_Affine)
    else
        D_Affine = mean(mean((SourcePoints_Affine-TargetPoints_Affine).^2,1).^0.5)
    end
    
    
    
    % AffineNonRigid
    Global_Matched_Source = FeaturePositions_A.Matched{sourceID,targetID}(:,1:3)';
    Global_Matched_Target = FeaturePositions_A.Matched{sourceID,targetID}(:,4:6)';
    [L,b]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target);
    
    SourcePoints_Affine = (L*(SourcePoints)'+b)';
    TargetPoints_Affine = TargetPoints;
    
    N_L=3;
    temp = (L*Global_Matched_Source+b);
    Min=min([min(temp,[],2),min(Global_Matched_Target,[],2)],[],2);
    Max=max([max(temp,[],2),max(Global_Matched_Target,[],2)],[],2);
    [~,XYZlmn,N_L,Min,Max]=Optimal_Nonrigid_Transform(temp,Global_Matched_Target,N_L,Min,Max);
    SourcePoints_Affine=Perform_Nonrigid_Transform(SourcePoints_Affine',XYZlmn,N_L,Min,Max)';
    if useTrace
        [Distances_NonRigid,~] = TraceDistance(AM_Source, SourcePoints_Affine, AM_Target, TargetPoints_Affine,pixelSize,0);
        D_NonRigidAffine = mean(Distances_NonRigid)
    else
        D_NonRigidAffine = mean(mean((SourcePoints_Affine-TargetPoints_Affine).^2,1).^0.5)
    end
    
    
    
    % NonRigid
    Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
    Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
    N_L=3;
    Min=min([min(Global_Matched_Source,[],2),min(Global_Matched_Target,[],2)],[],2);
    Max=max([max(Global_Matched_Source,[],2),max(Global_Matched_Target,[],2)],[],2);
    [~,XYZlmn,N_L,Min,Max]=Optimal_Nonrigid_Transform(Global_Matched_Source,Global_Matched_Target,N_L,Min,Max);
    SourcePoints_NR=Perform_Nonrigid_Transform(SourcePoints',XYZlmn,N_L,Min,Max);
    if useTrace
        [Distances_NonRigid,~] = TraceDistance(AM_Source, SourcePoints_NR', AM_Target, TargetPoints,pixelSize,0);
        D_NonRigid = mean(Distances_NonRigid);
    else
        D_NonRigid = mean(mean((SourcePoints_NR-Global_Matched_Target).^2,1).^0.5)
    end
    disp('-------------------------------------');
    
    AllDistances(i,1) = D_Before;
    AllDistances(i,2) = D_Translation;
    AllDistances(i,3) = D_Rigid;
    AllDistances(i,4) = D_Affine;
    AllDistances(i,5) = D_NonRigidAffine;
    AllDistances(i,6) = D_NonRigid;
end
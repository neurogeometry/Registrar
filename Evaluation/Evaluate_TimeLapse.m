% clear all;
% clc
% close all;
ppm = 2;
affine = 1;
addpath('NeuronTracerV20');
showDiff = 0;
csvPath = '..\..\RegistrationEvaluation\TimeLapse_Holtmaat_StackList.csv';
StackList = table2cell(readtable(csvPath,'Delimiter',','));

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
    
    
    
    if showDiff
        Source_Stack_File = char(StackList(sourceID,1));
        Target_Stack_File = char(StackList(targetID,1));
        IM_Source=ImportStack(char(Source_Stack_File));
        IM_Source = uint8(double(IM_Source)./double(max(IM_Source(:))).*255);
        IM_Target=ImportStack(char(Target_Stack_File));
        IM_Target = uint8(double(IM_Target)./double(max(IM_Target(:))).*255);
        IM_source_max=max(IM_Source,[],3);
        IM_target_max=max(IM_Target,[],3);
        
        
        
        %before
        figure(1)
        imshowpair(IM_source_max,IM_target_max,'Scaling','independent')
        
        %                     (1:500,1:500,1:100)
        
        [IM_Target_NR,StackPosition_prime,~]=Perform_Bspline_Transform(IM_Target,[1;1;1],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
        IM_Target_NR_max=max(IM_Target_NR,[],3);
        MIN=min([1;1],StackPosition_prime(1:2));
        MAX=max(size(IM_source_max)',size(IM_target_max)'+StackPosition_prime(1:2)-1);
        temp=zeros(MAX(1)-MIN(1)+1,MAX(2)-MIN(2)+1,'uint8');
        
        IM_Target_NR_max_P=temp;
        IM_Target_NR_max_P(StackPosition_prime(1)-MIN(1)+1:StackPosition_prime(1)-MIN(1)+size(IM_Target_NR_max,1),...
            StackPosition_prime(2)-MIN(2)+1:StackPosition_prime(2)-MIN(2)+size(IM_Target_NR_max,2))=IM_Target_NR_max;
        IM_source_max_P=temp;
        IM_source_max_P(2-MIN(1):1-MIN(1)+size(IM_source_max,1),...
            2-MIN(2):1-MIN(2)+size(IM_source_max,2))=IM_source_max;
        
        figure(5),imshowpair(IM_Target_NR_max_P,IM_source_max_P,'Scaling','independent')    
        
        %                     IM_source_maxR = imresize(IM_source_max,size(IM_Target_NR_max));
        %                     figure(3),imshowpair(IM_Target_NR_max,IM_source_maxR,'Scaling','independent')
        
%         [IM_Source_NR,~]=Perform_Bspline_Transform(IM_Source,[0,0,0],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
%         IM_Source_NR_max=max(IM_Source_NR,[],3);
%         figure(3)
%         imshowpair(IM_Source_NR_max,IM_target_max,'Scaling','independent')
%               
%         [optimizer, metric] = imregconfig('multimodal');
%         optimizer.InitialRadius = 0.009;
%         optimizer.Epsilon = 1.5e-4;
%         optimizer.GrowthFactor = 1.01;
%         optimizer.MaximumIterations = 300;
%         movingRegistered = imregister(IM_Source, IM_Target, 'affine', optimizer, metric);
%         figure(4)
%         imshowpair(IM_target_max,max(movingRegistered,[],3),'Scaling','independent')
    end
    
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





















% % Trnslation
%         Global_Matched_Source = FeaturePositions_T.Matched{SourceID,TargetID}(:,1:3)';
%         Global_Matched_Target = FeaturePositions_T.Matched{SourceID,TargetID}(:,4:6)';
%         b=Optimal_Translation_Transform(Global_Matched_Source,Global_Matched_Target);
%
%         SourcePoints_TranslationPairs = SourcePoints+b';
%         TargetPoints_TranslationPairs = TargetPoints;
%
%         if useTrace
%             [Distances_Translation,~] = TraceDistance(AM_Source, SourcePoints_TranslationPairs, AM_Target, TargetPoints_TranslationPairs,pixelSize,0);
%             D_Translation = mean(Distances_Translation)
%         else
%             D_Translation = mean(mean((SourcePoints_TranslationPairs-TargetPoints_TranslationPairs).^2,1).^0.5)
%         end
%
%         % Rigid
%         Global_Matched_Source = FeaturePositions_R.Matched{SourceID,TargetID}(:,1:3)';
%         Global_Matched_Target = FeaturePositions_R.Matched{SourceID,TargetID}(:,4:6)';
%         [R,b]=Optimal_Rigid_Transform(Global_Matched_Source,Global_Matched_Target);
%
%         SourcePoints_Rigid = (R*(SourcePoints)'+b)';
%         TargetPoints_Rigid = TargetPoints;
%         if useTrace
%             [Distances_Rigid,~] = TraceDistance(AM_Source, SourcePoints_Rigid, AM_Target, TargetPoints_Rigid,pixelSize,0);
%             D_Rigid = mean(Distances_Rigid)
%         else
%             D_Rigid = mean(mean((SourcePoints_Rigid-TargetPoints_Rigid).^2,1).^0.5)
%         end
%
%         % Affine
%         Global_Matched_Source = FeaturePositions_A.Matched{SourceID,TargetID}(:,1:3)';
%         Global_Matched_Target = FeaturePositions_A.Matched{SourceID,TargetID}(:,4:6)';
%         [L,b]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target);
%
%         SourcePoints_Affine = (L*(SourcePoints)'+b)';
%         TargetPoints_Affine = TargetPoints;
%         if useTrace
%             [Distances_Affine,~] = TraceDistance(AM_Source, SourcePoints_Affine, AM_Target, TargetPoints_Affine,pixelSize,0);
%             D_Affine = mean(Distances_Affine)
%         else
%             D_Affine = mean(mean((SourcePoints_Affine-TargetPoints_Affine).^2,1).^0.5)
%         end
%
%
%
%         % AffineNonRigid
%         Global_Matched_Source = FeaturePositions_A.Matched{SourceID,TargetID}(:,1:3)';
%         Global_Matched_Target = FeaturePositions_A.Matched{SourceID,TargetID}(:,4:6)';
%         [L,b]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target);
%
%         SourcePoints_Affine = (L*(SourcePoints)'+b)';
%         TargetPoints_Affine = TargetPoints;
%
%         N_L=3;
%         temp = (L*Global_Matched_Source+b);
%         Min=min([min(temp,[],2),min(Global_Matched_Target,[],2)],[],2);
%         Max=max([max(temp,[],2),max(Global_Matched_Target,[],2)],[],2);
%         [~,XYZlmn,N_L,Min,Max]=Optimal_Nonrigid_Transform(temp,Global_Matched_Target,N_L,Min,Max);
%         SourcePoints_Affine=Perform_Nonrigid_Transform(SourcePoints_Affine',XYZlmn,N_L,Min,Max)';
%         if useTrace
%             [Distances_NonRigid,~] = TraceDistance(AM_Source, SourcePoints_Affine, AM_Target, TargetPoints_Affine,pixelSize,0);
%             D_NonRigidAffine = mean(Distances_NonRigid)
%         else
%             D_NonRigidAffine = mean(mean((SourcePoints_Affine-TargetPoints_Affine).^2,1).^0.5)
%         end
%
%
%
%         % NonRigid
%         Global_Matched_Source = FeaturePositions_NR.Matched{SourceID,TargetID}(:,1:3)';
%         Global_Matched_Target = FeaturePositions_NR.Matched{SourceID,TargetID}(:,4:6)';
%         N_L=3;
%         Min=min([min(Global_Matched_Source,[],2),min(Global_Matched_Target,[],2)],[],2);
%         Max=max([max(Global_Matched_Source,[],2),max(Global_Matched_Target,[],2)],[],2);
%         [~,XYZlmn,N_L,Min,Max]=Optimal_Nonrigid_Transform(Global_Matched_Source,Global_Matched_Target,N_L,Min,Max);
%         SourcePoints_NR=Perform_Nonrigid_Transform(SourcePoints',XYZlmn,N_L,Min,Max);
%         if useTrace
%             [Distances_NonRigid,~] = TraceDistance(AM_Source, SourcePoints_NR', AM_Target, TargetPoints,pixelSize,0);
%             D_NonRigid = mean(Distances_NonRigid);
%         else
%             D_NonRigid = mean(mean((SourcePoints_NR-Global_Matched_Target).^2,1).^0.5)
%         end
%         disp('-------------------------------------');
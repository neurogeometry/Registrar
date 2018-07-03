% clear all;
% clc
% close all;
ppm = 1;
affine = 1;
addpath('NeuronTracerV20');
showDiff = 0;
csvPath = '..\..\RegistrationEvaluation\TimeLapse_Holtmaat_StackList.csv';
StackList = table2cell(readtable(csvPath,'Delimiter',','));

Affine_Fiji{1} =   [0.9936921821167053, 0.008007812145490931, 0.0315078570167579, -24.27614925202483;-0.0054197055191250315, 0.993164734920597, 9.180206475727926E-4, 0.9227623205139378;0.010330514866540269, 3.326481269751852E-4, 0.9955273091289829, 13.451674407899446];
Affine_Fiji{2} =   [1.0017279611119831, -0.008739748867177977, 7.598120178636258E-4, 7.201407511322642;0.009988781967444614, 0.9909520362960338, 0.016058198093810594, -0.5843885119905292;0.0019325126371294576, -0.0015279713437594648, 0.9697108879368073, -5.7824216586693495];
Affine_Fiji{3} =   [0.9844017709049991, -0.004106107397898372, 0.05665728500606049, 4.011870426481522;0.002314511903363372, 0.9976149176091117, 0.05721206380947885, 7.477534636215167;2.0194609914565842E-4, 0.0028578837113763713, 0.9947226573336245, 0.834333092281156];
Affine_Fiji{4} =   [0.9954810548175131, 9.28705942175697E-4, -0.02071482289753468, 0.04078623784920676;0.001797893952538314, 0.9956744162603566, -0.006329119904485814, -3.156202588085548;-0.001888915946940975, -9.880429642617553E-4, 0.993944675667947, 1.654007709967047];
Affine_Fiji{5} =   [1.0026400364341825, -0.01351132543444835, -0.0034545862571120356, 8.423251468402011; 0.012552151208417454, 0.9928482947217806, 0.04999547881086985, -3.512752634002786; -6.879998375250163E-4, -0.0018044605030352878, 1.00283949298545, 4.691198095246264];
Affine_Fiji{6} =   [1.00001500997013, 0.0033570535987453437, 1.655919355156854E-4, -5.623941605866036; 3.151418458432198E-4, 1.0003423716493687, 0.02932621523650758, -4.074664038821975; 0.00231113287979684, 5.64096714365081E-4, 0.9863876661587397, -4.489086279849488];
Affine_Fiji{7} =   [0.9971674331463407, -0.006928105953497653, 0.02093248610197307, 3.78219629909988; 0.0012900147841593476, 0.9968799828718048, -0.014788309514328968, 2.122652745631963; 0.0012393347579878963, 3.341368993613203E-5, 0.99157193992645, 4.480019454011426];
Affine_Fiji{8} =   [0.9994495252473322, 0.005602715955703927, 0.010142786923411196, -5.073125721516959; -0.003309083440242014, 0.9985107591039083, 0.03356060956302653, 3.5201961790928458; -0.003021802481539562, -0.00327379981313812, 0.9974418907957041, 2.3091148114248625];
Affine_Fiji{9} =   [0.9960913049019193, -0.002631490674806336, -0.04185196685120762, 3.2674121530160285; 2.4524281739611575E-4, 1.003278303407285, -0.03460646887457852, 3.3864888411433762; 0.0044193437006529585, 0.002611826289340785, 1.0225941005220958, 5.327124369787498];
Affine_Fiji{10} =  [0.9927585128466225, 0.02331069556632598, -0.03705097060629248, -21.60857613836826; -0.018675009003586276, 0.9987739155669406, -0.0480350149519263, 16.000861986977647; -0.006566536391762391, -3.921913255836018E-4, 0.9721557587333788, 1.2558925790011557];
Affine_Fiji{11} =  [0.9863087955784957, -0.013385167092749258, -0.03964070888656113, 14.191123332558618; -0.005094676029152205, 1.0019087385416316, 0.028051602718283688, -7.579518932245131; 0.003347758456653957, -0.0016588264544113145, 0.9960583050214211, -2.3971075801555344];
Affine_Fiji{12} =  [0.9861824146461196, -0.009038573199904351, 0.014602462685948447, 13.867741740447787; 0.033014489240465995, 1.0094896667987125, 0.07722290456304083, -15.21638104893514; -0.011171966867808941, 0.0020511317808233898, 1.019145093243156, 0.6197985081796844];
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
mu_Names = {'0','10'};
Dis_NonRigid_voxelall = [];
% Show Traces
% sourceID = 2;
% targetID = 3;
tic
for nummu = 1:size(mu_Names,2)
for sourceID = 1:12
    %     disp(sourceID)
    targetID = sourceID + 1;
    FeaturePositions_NR = load([resultpath,'MatchedPoints_Non-Rigid_mu',mu_Names{nummu},'.mat']);
    
    % B-Spline-NonRigid
    Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
    Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
    Minimum = min(Global_Matched_Source,[],2);
    Maximum = max(Global_Matched_Source,[],2);
    nxyz = [512;512;156]; %image size
    %     Nxyz = ceil((Maximum-Minimum)./nxyz');
    [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,str2double(mu_Names{nummu}));
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
    BoutonsBefore = mean(mean((TargetPoints-SourcePoints).^2,1).^0.5)
    
    %         D_NonRigid_voxel = mean(mean((SourcePoints_NR'-TargetPoints).^2,1).^0.5);
    %         D_NonRigid_voxel = SourcePoints_NR'-TargetPoints;
    [SourcePoints_NR,~]=Perform_Bspline_Transform(SourcePoints',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
    result{sourceID}.Bouton.r1 = SourcePoints;
    result{sourceID}.Bouton.r1_NR = SourcePoints_NR';
    result{sourceID}.Bouton.r2 = TargetPoints;
    %     mean(mean((TargetPoints-SourcePoints_NR').^2,1).^0.5)
    BoutonsNonrigid = mean(mean((TargetPoints-SourcePoints_NR').^2,1).^0.5)
    
    Affine_L = Affine_Fiji{sourceID}(:,1:3);
    Affine_b = Affine_Fiji{sourceID}(:,4);
    
    
    SourcePoints_Affine = ((SourcePoints(:,[2,1,3])));
    TargetPoints_Affine = (Affine_L*TargetPoints(:,[2,1,3])'+Affine_b)';%(Affine_L*(TargetPoints)'-Affine_b)';
    result{sourceID}.Bouton.r1_fiji_Affine = SourcePoints_Affine;
    result{sourceID}.Bouton.r2_fiji_Affine = TargetPoints_Affine;
    BoutonsFiji = mean(mean((SourcePoints_Affine-TargetPoints_Affine).^2,1).^0.5)
    %         Affine_b = -Affine_b_temp([2,1,3])
    %         b
    %     result{sourceID}.Bouton.r1_fiji_Affine = (Affine_L*(SourcePoints)'+Affine_b)';
    %         TargetPoints_Affine = TargetPoints;
    %         SourcePoints_Affine'-TargetPoints_Affine
    
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
%         [Distances_before,~] = TraceDistance(AM_Source, SourcePoints, AM_Target, TargetPoints,pixelSize,0);

        
        
        
        
        
        
        
        
        
        
        
        
        
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
        %         [Dis_NonRigid_um,Dis_NonRigid_voxel] = TraceDistance(AM_Source, SourcePoints_NR', AM_Target, TargetPoints,pixelSize,0);
        %         Dis_NonRigid_voxelall(i)=Dis_NonRigid_voxel
        %----------- Fiji
        
        %Affine
        
        SourcePoints_Affine = ((SourcePoints(:,[2,1,3])));
        TargetPoints_Affine = (Affine_L*TargetPoints(:,[2,1,3])'+Affine_b)';%(Affine_L*(TargetPoints)'-Affine_b)';
        %         [Affine_Dis_um,Affine_Dis_voxel] = TraceDistance(AM_Source, SourcePoints_Affine, AM_Target, TargetPoints_Affine,pixelSize,0);
        %         Affine_Dis_voxelall(i) = Affine_Dis_voxel
        result{sourceID}.Trace.r1_fiji_Affine{i} = SourcePoints_Affine;
        result{sourceID}.Trace.r2_fiji_Affine{i} = TargetPoints_Affine;
        %         result{sourceID}.Trace.r1_fiji_Affine{i} = (Affine_L*(SourcePoints)'+Affine_b)';
        %         TargetPoints_Affine = TargetPoints;
        
        % [Before_Dis_um,Before_Dis_voxel] = TraceDistance(AM_Source, SourcePoints, AM_Target, TargetPoints, pixelSize, 0);
        % Before_Dis_voxel
        
        %         SourcePoints_Affine = ((SourcePoints)+Affine_b');
        %         TargetPoints_Affine = (Affine_L*TargetPoints(:,[2,1,3])')';%(Affine_L*(TargetPoints)'-Affine_b)';
        %         TargetPoints_Affine = TargetPoints_Affine(:,[2,1,3]);
        %         [Affine_Dis_um,Affine_Dis_voxel] = TraceDistance(AM_Source, SourcePoints_Affine, AM_Target, TargetPoints_Affine,pixelSize,0);
        % Affine_Dis_voxel
        %          PlotAM(AM_Source,SourcePoints_Affine,'r')
        %                 axis equal
        %                 PlotAM(AM_Target,TargetPoints_Affine,'g')
        %                 axis equal
        
        
        
        
        % l=1
        
        
        %         if useTrace
        %             [Distances_Affine,~] = TraceDistance(AM_Source, SourcePoints_Affine', AM_Target, TargetPoints_Affine',pixelSize,0);
        %             D_Affine = mean(Distances_Affine)
        %         else
        %             D_Affine = mean(mean((SourcePoints_Affine-TargetPoints_Affine).^2,1).^0.5)
        %         end
        
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
save(['..\..\RegistrationEvaluation\result_NR_fiji_',mu_Names{nummu},'.mat'],'result','pixelSize')
end
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
% clear all;
% clc
% close all;
ppm = 2;
affine = 1;
addpath('NeuronTracerV20');
showDiff = 0;
csvPath = '..\..\RegistrationEvaluation\TimeLapse_Holtmaat_StackList.csv';
StackList = table2cell(readtable(csvPath,'Delimiter',','));

Affine_Fiji{1} =  [1.0069332813089047, -0.00854589834317538, -0.03260691097378371, 25.087531465827908; 
0.005935560672269434, 1.00641289904624, -7.341383252498224E-4, -0.8644961248943994; 
-0.010410748719979171, -4.036321915032448E-4, 1.0052884173593035, -13.772167004564864];
Affine_Fiji{2} =  [0.9983291029115744, 0.008773565768463169, -0.001373618163566681, -7.140981589639295; 
    -0.0099040317030493, 1.009027825065207, -0.015988862814655364, 0.4232911614558521; 
    -0.0019108435254093711, 0.0017134007740634552, 1.0315824267901006, 5.824324365351288];       
Affine_Fiji{3} =  [1.0157573748738304, 0.00415654503524179, -0.05758314461282188, -3.983125058716479;
    -0.0024844422778687383, 1.0023028206740348, -0.05773061104034549, -7.283894609986068;
    -2.4568539617978585E-4, -0.0029588706054060956, 1.0055313704308426, -0.7932550093767929];         
Affine_Fiji{4} =  [0.9954810548175131, 9.28705942175697E-4, -0.02071482289753468, 0.04078623784920676;
    0.001797893952538314, 0.9956744162603566, -0.006329119904485814, -3.156202588085548;
    -0.001888915946940975, -9.880429642617553E-4, 0.993944675667947, 1.654007709967047]; 
Affine_Fiji{5} =  [1.0026400364341825, -0.01351132543444835, -0.0034545862571120356, 8.423251468402011;
    0.012552151208417454, 0.9928482947217806, 0.04999547881086985, -3.512752634002786;
    -6.879998375250163E-4, -0.0018044605030352878, 1.00283949298545, -4.691198095246264]; 
Affine_Fiji{6} =  [1.0000110899235126, -0.0033838858919999076, -1.9264684199549187E-5, 5.6187768903928115;
    -4.5456432500028815E-4, 0.9997120138287525, -0.029755361036161165, 4.244980057741089;
    -0.002278505893924164, -5.063353878100585E-4, 1.0137185236019786, -4.606189806547832]; 
Affine_Fiji{7} =  [1.0026927764455296, 0.0069590405075884015, -0.02055518045961091, -3.722207026804554;
    -0.001322915753459289, 1.0031868854742882, 0.014634420090856137, -2.181350452959463;
    -0.0011290052809648163, 9.111217662877147E-5, 1.0085835186456233, -4.612587163461214];
 Affine_Fiji{8} =  [0.9994495252473322, 0.005602715955703927, 0.010142786923411196, -5.073125721516959;
     -0.003309083440242014, 0.9985107591039083, 0.03356060956302653, 3.5201961790928458;
     -0.003021802481539562, -0.00327379981313812, 0.9974418907957041, 2.3091148114248625];
 Affine_Fiji{9} =  [1.004072835731811, 0.0024854782393286395, 0.041581630657794966, -3.173496556765299;
     -2.804664633419829E-4, 0.9967128858467997, 0.03403311224415906, -3.3139990321841015;
     -0.0044491647758947706, -0.0027637141205848273, 0.9777999091126199, 5.345392011504245]; 
 Affine_Fiji{10} =  [1.0071443003151253, -0.02349817173098602, 0.0372157383038727, 22.195378410021682;
     0.019417956473857377, 1.0007603172016564, 0.050720146585283676, -15.661595364287251;
     0.006939088789655483, 2.0629444192779567E-4, 1.0288819188079352, 1.4087576181626957];
 Affine_Fiji{11} =  [1.0137174032775829, 0.013061027721580418, 0.038120924863384587, -13.802574567119578;
     0.005150729724184285, 0.9979231639239432, -0.02775927334378541, 7.657849702262158;
     -0.003404912468337795, 0.001478880608151754, 1.0035055919922646, -2.251670111535077]; 
 Affine_Fiji{12} =  [1.0139986568140658, 0.009573975132356805, -0.01446398825983286, -14.415653029635457;
     -0.03409538428580029, 0.9901566907991396, -0.07465939906395067, 15.642695011823143;
     0.010985289730242875, -0.0018223180411751985, 0.980679918183361, 0.48934445024409];
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
    
    
        Affine_L = Affine_Fiji{sourceID}(:,1:3);
        Affine_b = Affine_Fiji{sourceID}(:,4);
        SourcePoints_Affine = (Affine_L*(SourcePoints)'+Affine_b)';
        TargetPoints_Affine = TargetPoints;
        SourcePoints_Affine'-TargetPoints_Affine
    
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
        
        
        
        
        
        %----------- Fiji
        
        %Affine
        
        SourcePoints_Affine = (Affine_L*(SourcePoints)'+Affine_b)';
        TargetPoints_Affine = TargetPoints;
% %         [Distances_Affine,~] = TraceDistance(AM_Source, SourcePoints_Affine, AM_Target, TargetPoints_Affine,pixelSize,0);
%         if useTrace
%             [Distances_Affine,~] = TraceDistance(AM_Source, SourcePoints_Affine', AM_Target, TargetPoints_Affine',pixelSize,0);
%             D_Affine = mean(Distances_Affine)
%         else
%             D_Affine = mean(mean((SourcePoints_Affine-TargetPoints_Affine).^2,1).^0.5)
%         end
         








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
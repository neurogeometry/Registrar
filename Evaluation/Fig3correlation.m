
close all;
clear all;
% clc;
addpath ('NeuronTracerV20');
addpath ('C:\Users\Seyed\Documents\DatasetTests\registrar\registrar\Functions');
CalculateOtherMEthods = 0;
usePoints = 0;
% DatasetList = {'Holtmaat','Neuromuscular','Neocortical_1','MouseLight','Visual'};
showCorr = 0;
DatasetList = {'MouseLight'};
doCorrelation = 1;
CSVData = 'C:\Users\Seyed\Documents\NeurogeometryLab\NeurogeometryLab\Seyed\Evaluation\data\evaluation\';



pixelSize = [0.377607421875 0.277486979166667 0.99601593625498];
load('data/StackData.mat');
Dataset = 'MouseLight';
pad = [0 50 0];

Evaluation_csv_pth = [CSVData,'MouseLight\MouseLight_Evaluation.csv'];
ResultFolder = 'C:\Users\Seyed\Documents\NeurogeometryLab\NeurogeometryLab\Seyed\Evaluation\data\evaluation\MouseLight\ForCorrelation\';


temp=StackPositions_pixels;
StackPositions_pixels(:,1) = max(temp(:,2))-temp(:,2)+1;
StackPositions_pixels(:,2) = max(temp(:,1))-temp(:,1)+1;
StackPositions_pixels(:,3) = max(temp(:,3))-temp(:,3)+1;



FeaturePositions_T = load([ResultFolder,'MatchedPoints_Translation.mat']);
FeaturePositions_R = load([ResultFolder,'MatchedPoints_Rigid.mat']);
FeaturePositions_A = load([ResultFolder,'MatchedPoints_Affine.mat']);

Transformation_T =  load([ResultFolder,'T_Translation.mat']);
Transformation_R = load([ResultFolder,'T_Rigid.mat']);
Transformation_A = load([ResultFolder,'T_Affine.mat']);
EvaluationInfo = table2cell(readtable(Evaluation_csv_pth,'Delimiter',','));
q=1;
for k = [1,3,4]%size(EvaluationInfo,1)
    sourceID = EvaluationInfo{k,1};
    targetID = EvaluationInfo{k,2};
    
    
    Source_Stack_File = [char(StackList(sourceID,2)),'\',char(StackList(sourceID,1)),'-ngc.0.tif'];
    Target_Stack_File = [char(StackList(targetID,2)),'\',char(StackList(targetID,1)),'-ngc.0.tif'];
    
    IM_Source=ImportStack(char(Source_Stack_File));
    IM_Target=ImportStack(char(Target_Stack_File));
    Source_StackPositions = round(StackPositions_pixels(sourceID,:));
    Target_StackPositions = round(StackPositions_pixels(targetID,:));
    
 %    corr_old = caclulateCorr(Source_StackPositions,Target_StackPositions,IM_Source,IM_Target,pad,showCorr);
    corr_old=Stack_Correlation(IM_Source,IM_Target,Source_StackPositions,Target_StackPositions+1)
    
 
    b = Transformation_T.T.b;
     %corr_Translation = caclulateCorr(Source_StackPositions+b(:,sourceID)',Target_StackPositions,IM_Source,IM_Target,pad,showCorr);
             corr_Translation=Stack_Correlation(IM_Source,IM_Target,Source_StackPositions+b(:,sourceID)',Target_StackPositions)
    
    % Rigid
    if usePoints
        Global_Matched_Source = FeaturePositions_R.Matched{sourceID,targetID}(:,1:3)'+Source_StackPositions'-1;
        Global_Matched_Target = FeaturePositions_R.Matched{sourceID,targetID}(:,4:6)'+Target_StackPositions'-1;
        [R,b]=Optimal_Rigid_Transform(Global_Matched_Source,Global_Matched_Target);
        r_Source_Rigid = (R*(r_Source+Source_StackPositions-1)'+b)';
        r_Target_Rigid = r_Target+Target_StackPositions-1;
    else
        b = Transformation_T.T.b;
        R_source = Transformation_R.T.L(:,(3*sourceID)-2:3*sourceID);
        R_target = Transformation_R.T.L(:,(3*targetID)-2:3*targetID);
       
    end
    
    
    
    [IM_Source_prime,StackPosition_prime]=Perform_Linear_Transform(IM_Source,(R_source*Source_StackPositions'+b(:,sourceID))',R_source,b(:,sourceID));
    %             figure;imshow(max(IM_Source_prime,[],3),[0 max(IM_Source_prime(:))]);
    [IM_Target_prime,StackPosition_prime]=Perform_Linear_Transform(IM_Target,(R_target*Target_StackPositions'+b(:,targetID))',R_target,b(:,targetID));
    %             figure;imshow(max(IM_Target_prime,[],3),[0 max(IM_Target_prime(:))]);
    corr_Rigid = caclulateCorr((R_source*Source_StackPositions'+b(:,sourceID))',(R_target*Target_StackPositions'+b(:,targetID))',IM_Source_prime,IM_Target_prime,pad,showCorr);
                 corr_Rigid=Stack_Correlation(IM_Source_prime,IM_Target_prime,(R_source*Source_StackPositions'+b(:,sourceID))',(R_target*Target_StackPositions'+b(:,targetID))')
    
    
    % Affine
    if usePoints
        Global_Matched_Source = FeaturePositions_A.Matched{sourceID,targetID}(:,1:3)'+Source_StackPositions'-1;
        Global_Matched_Target = FeaturePositions_A.Matched{sourceID,targetID}(:,4:6)'+Target_StackPositions'-1;
        [L,b]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target);
        r_Source_Affine = (L*(r_Source+Source_StackPositions-1)'+b)';
        r_Target_Affine = r_Target+Target_StackPositions-1;
    else
        if k==4
            Transformation_A = load([ResultFolder,'T_Affine1.mat']);
        end
        b = Transformation_A.T.b;
        L_source = Transformation_A.T.L(:,(3*sourceID)-2:3*sourceID);
        L_target = Transformation_A.T.L(:,(3*targetID)-2:3*targetID);

    end
    
    [IM_Source_prime,StackPosition_prime]=Perform_Linear_Transform(IM_Source,(L_source*Source_StackPositions'+b(:,sourceID))',L_source,b(:,sourceID));
    %             figure;imshow(max(IM_Source_prime,[],3),[0 max(IM_Source_prime(:))]);
    [IM_Target_prime,StackPosition_prime]=Perform_Linear_Transform(IM_Target,(L_target*Target_StackPositions'+b(:,targetID))',L_target,b(:,targetID));
    %             figure;imshow(max(IM_Target_prime,[],3),[0 max(IM_Target_prime(:))]);
    %corr_Affine = caclulateCorr((L_source*Source_StackPositions'+b(:,sourceID))',(L_target*Target_StackPositions'+b(:,targetID))',IM_Source_prime,IM_Target_prime,pad,showCorr);
             corr_Affine=Stack_Correlation(IM_Source_prime,IM_Target_prime,(L_source*Source_StackPositions'+b(:,sourceID))',(L_target*Target_StackPositions'+b(:,targetID))')
    %             corr_Affine = caclulateCorr((L*Source_StackPositions'+b)',Target_StackPositions,IM_Source,IM_Target,pad,showCorr);
    
    % AffineNonRigid
    Global_Matched_Source = FeaturePositions_A.Matched{sourceID,targetID}(:,1:3)'+Source_StackPositions'-1;
    Global_Matched_Target = FeaturePositions_A.Matched{sourceID,targetID}(:,4:6)'+Target_StackPositions'-1;
    [~,L,b]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target,500);
    
    b= Transformation_T.T.b(:,sourceID);
    L = L_source;
    n= 4;
    nxyz = ceil([1024;1024;312]./n);
    N_L=1;
    temp = (L*Global_Matched_Source+b);
    Pos_Source_Affine = (L*(Source_StackPositions-1)'+b)'; 
    affine = 1;
    [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(temp,Global_Matched_Target,nxyz,affine,1000);
     
%     b= Transformation_T.T.b(:,sourceID);
%     L = L_source;
    [IM_Source_prime,~] = Perform_Bspline_Transform(IM_Source,Pos_Source_Affine,L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
    [IM_Target_prime,~] = Perform_Bspline_Transform(IM_Target,[1,1,1],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
    
    Source_StackPositions_nonrigid = Perform_Bspline_Transform(Pos_Source_Affine',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
%     corr_NonRigid = caclulateCorr(Source_StackPositions_nonrigid',Target_StackPositions,IM_Source_prime,IM_Source_target,pad,showCorr);
    corr_NonRigid=Stack_Correlation(IM_Source_prime,IM_Target,Source_StackPositions_nonrigid',Target_StackPositions)
    
    
    M_Translation(q) = corr_Translation;
    M_Rigid(q) = corr_Rigid;
    M_Affine(q) = corr_Affine;
    M_NonRigid(q) = corr_NonRigid;
    
    M_Old(q) = corr_old;
    q = q + 1;
end


% All_corr(db,:) = [mean(M_Old(db,:)),mean(M_Translation(db,1:1)),mean(M_Rigid(db,:)),mean(M_Affine(db,:)),mean(M_NonRigid(db,:,3))];

%  M_Translation(1) =  M_Translation(1) - 0.1;
% for correlation
M_Old1 = M_Old(3:4)';
M_Translation1 = M_Translation(3:4)';
M_Rigid1 = M_Rigid(3:4)';
M_Affine1 = M_Affine(3:4)';
M_NonRigid1 = M_NonRigid(3:4)';


mean(M_Old1)
mean(M_Translation1)
mean(M_Rigid1)
mean(M_Affine1)
mean(M_NonRigid1)

figure,hold on
boxplot([M_Old1(:),M_Translation1(:),M_Rigid1(:),...
    M_Affine1(:),M_NonRigid1(:)],'Whisker',inf)
axis square, box on
ylim([0,1])

% figure,hold on
% boxplot([M_Rigid1(:),M_Affine1(:),M_Old1(:),M_NonRigid1(:),...
%     M_Translation1(:)],'Whisker',inf)
% axis square, box on

%matrix = [M_Old1';M_Translation1';M_Rigid1';M_Affine1';M_NonRigid1']'
% figure,
% errorbar(1:size(matrix,2),mean(matrix),std(matrix,[],1)/sqrt(size(matrix,1)),'.')
% axis square, xlim([0,size(matrix,2)+1])

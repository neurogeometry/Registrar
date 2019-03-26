
close all;
clear all;
% clc;
addpath ('NeuronTracerV20');
addpath ('C:\Users\Seyed\Documents\DatasetTests\registrar\registrar\Functions');

showCorr = 0;
pad = [0 60 0];

csv_pth = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\CorrelationResults\';
files = dir(csv_pth);
% 
% k = 1;
% correllation = {};
% for i = 3:size(files,1)-1
%     correllation{k} = load([files(i).folder,'\',files(i).name,'\','corr.mat']);
%     k=k+1;
% end
% j = 1;
% for i=1:4
%     O(j) = correllation{i}.corr_old;
%     T(j) = correllation{i}.corr_Translation;
%     R(j) = correllation{i}.corr_Rigid;
%     A(j) = correllation{i}.corr_Affine;
%     NR(j) = correllation{i}.corr_NonRigid;
%     j = j +1
% end
% 
% boxplot([O(:),T(:),R(:),...
%     A(:),NR(:)],'Whisker',inf)
% axis square, box on

for i = 3:size(files,1)-1
    
    ResultFolder = [files(i).folder,'\',files(i).name,'\Results-',files(i).name,'\'];
    csvFile = [files(i).folder,'\',files(i).name,'\',files(i).name,'.csv'];
    PositionsCSV = table2cell(readtable(csvFile,'Delimiter',','));
    StackPositions_pixels = cell2mat(PositionsCSV(:,2:4));
    
    
    temp=StackPositions_pixels;
    StackPositions_pixels(:,1) = max(temp(:,2))-temp(:,2);
    StackPositions_pixels(:,2) = max(temp(:,1))-temp(:,1);
    StackPositions_pixels(:,3) = max(temp(:,3))-temp(:,3); % Not Sure - Need to check
    
    
    Transformation_T =  load([ResultFolder,'T_Translation.mat']);
    Transformation_R = load([ResultFolder,'T_Rigid.mat']);
    Transformation_A = load([ResultFolder,'T_Affine.mat']);
    
    FeaturePositions_T = load([ResultFolder,'MatchedPoints_Translation.mat']);
    FeaturePositions_R = load([ResultFolder,'MatchedPoints_Rigid.mat']);
    FeaturePositions_A = load([ResultFolder,'MatchedPoints_Affine.mat']);
    
    
    sourceID = 1;
    targetID = 2;
    
    
    Source_Stack_File = PositionsCSV{sourceID,1};
    Target_Stack_File = PositionsCSV{targetID,1};
    
    IM_Source=ImportStack(char(Source_Stack_File));
    IM_Target=ImportStack(char(Target_Stack_File));
    Source_StackPositions = round(StackPositions_pixels(sourceID,:));
    Target_StackPositions = round(StackPositions_pixels(targetID,:));
    
    % % --------------------------------------------------------------------------------------------------
    %before
    corr_old=Stack_Correlation(IM_Source,IM_Target,Source_StackPositions,Target_StackPositions)
    
    % --------------------------------------------------------------------------------------------------
    % Translation
    b = Transformation_T.T.b;
    corr_Translation=Stack_Correlation(IM_Source,IM_Target,Source_StackPositions+b(:,sourceID)',Target_StackPositions)
    
    
    
    % --------------------------------------------------------------------------------------------------
    % Rigid
    b = Transformation_R.T.b;
    R_source = Transformation_R.T.L(:,(3*sourceID)-2:3*sourceID);
    R_target = Transformation_R.T.L(:,(3*targetID)-2:3*targetID);
    
    [IM_Source_prime,SourceStackPosition_prime]=Perform_Linear_Transform(IM_Source,(R_source*Source_StackPositions'+b(:,sourceID))',R_source,b(:,sourceID));
    [IM_Target_prime,TargetStackPosition_prime]=Perform_Linear_Transform(IM_Target,(R_target*Target_StackPositions'+b(:,targetID))',R_target,b(:,targetID));
    
    corr_Rigid=Stack_Correlation(IM_Source_prime,IM_Target_prime,SourceStackPosition_prime,TargetStackPosition_prime)
    % corr_Rigid=Stack_Correlation(IM_Source_prime,IM_Target_prime,(R_source*Source_StackPositions'+b(:,sourceID))',(R_target*Target_StackPositions'+b(:,targetID))')
    
    % --------------------------------------------------------------------------------------------------
    % Affine
    
    b = Transformation_A.T.b;
    L_source = Transformation_A.T.L(:,(3*sourceID)-2:3*sourceID);
    L_target = Transformation_A.T.L(:,(3*targetID)-2:3*targetID);
    
    [IM_Source_prime,StackPosition_prime]=Perform_Linear_Transform(IM_Source,(L_source*Source_StackPositions'+b(:,sourceID))',L_source,b(:,sourceID));
    [IM_Target_prime,StackPosition_prime]=Perform_Linear_Transform(IM_Target,(L_target*Target_StackPositions'+b(:,targetID))',L_target,b(:,targetID));
    corr_Affine=Stack_Correlation(IM_Source_prime,IM_Target_prime,(L_source*Source_StackPositions'+b(:,sourceID))',(L_target*Target_StackPositions'+b(:,targetID))')
    
    % --------------------------------------------------------------------------------------------------
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
    [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(temp,Global_Matched_Target,nxyz,affine,1024);
    
    %     b= Transformation_T.T.b(:,sourceID);
    %     L = L_source;
    [IM_Source_prime,~] = Perform_Bspline_Transform(IM_Source,Pos_Source_Affine,L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
    [IM_Target_prime,~] = Perform_Bspline_Transform(IM_Target,[1,1,1],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
    
    Source_StackPositions_nonrigid = Perform_Bspline_Transform(Pos_Source_Affine',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
    corr_NonRigid=Stack_Correlation(IM_Source_prime,IM_Target,Source_StackPositions_nonrigid',Target_StackPositions)
    
    
    % corr_old
    % corr_Translation
    % corr_Rigid
    % corr_Affine
    % corr_NonRigid
    
    save([files(i).folder,'\',files(i).name,'\','corr.mat'],'corr_old','corr_Translation','corr_Rigid','corr_Affine','corr_NonRigid');
end

k = 1;
correllation = {};
for i = 3:size(files,1)-1
    correllation{k} = load([files(i).folder,'\',files(i).name,'\','corr.mat']);
    k=k+1;
end

for i=1:4
    O(i) = correllation{i}.corr_old;
    T(i) = correllation{i}.corr_Translation;
    R(i) = correllation{i}.corr_Rigid;
    A(i) = correllation{i}.corr_Affine;
    NR(i) = correllation{i}.corr_NonRigid;
    
end

boxplot([O(:),T(:),R(:),...
    A(:),NR(:)],'Whisker',inf)
axis square, box on
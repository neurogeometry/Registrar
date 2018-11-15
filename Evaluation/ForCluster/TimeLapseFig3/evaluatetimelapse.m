function evaluatetimelapse(ID)
ID = 1
ppm = 1;
affine = 1;
addpath('NeuronTracerV20');
addpath('Functions');


pixelSize = [0.26 0.26 0.8];

% load('../data/StackData_TimeLapse_Holtmaat.mat');
load('DL083-001-Matches.mat');

mu = [0,2.^(0:0.5:25)];

CutLength=100;

FeaturePositions_NR = load('MatchedPoints_Non-Rigid_mu1.mat');
T_Names = {'B','C','D','E','F','G','H','I','J','K','L','M','N'};

n= 4;
nxyz = ceil([1024;1024;312]./n);

sourceID = ID ;
targetID = sourceID + 1;
Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
fname_First = dir(['DL083',T_Names{sourceID},'001-A0*']);
fname_First={fname_First.name}';
fname_Second = dir(['DL083',T_Names{targetID},'001-A0*']);
fname_Second={fname_Second.name}';

%Use Boutons

%         B1 = Boutons{sourceID,1};
%         BSourcePoints = B1.r1;
%         BTargetPoints = B1.r2;
%         BoutonsDistanceOriginal = mean(mean((BTargetPoints-BSourcePoints).^2,2).^0.5)
%

for i=1:size(fname_First,1)
    sourcePath = [fname_First{i}];%'E:\Datasets\TimeLaps\Matches\Traces\DL083B001-A001.swc';
    targetPath = [fname_Second{i}];%'E:\Datasets\TimeLaps\Matches\Traces\DL083C001-A001.swc';
    
    [AM_Source_temp,r_Source_temp,~]=swc2AM_fast(sourcePath);
    [AM_Target_temp,r_Target_temp,~]=swc2AM_fast(targetPath);
    [AM_Source_temp,r_Source_temp,~] = AdjustPPM_fast(AM_Source_temp,r_Source_temp,zeros(size(r_Source_temp,1),1),ppm);
    [AM_Target_temp,r_Target_temp,~] = AdjustPPM_fast(AM_Target_temp,r_Target_temp,zeros(size(r_Target_temp,1),1),ppm);
    
    SourcePoints{i} = r_Source_temp;
    TargetPoints{i} = r_Target_temp;
    
    Ncuts=fix(size(AM_Source_temp,1)./CutLength);
    inds=ceil(size(AM_Source_temp,1)./(Ncuts+1)).*(1:Ncuts);
    
    AM_Source_temp(inds,:)=0;
    AM_Source_temp(:,inds)=0;
    AM_Source{i}=LabelBranchesAM(AM_Source_temp);
    
    AM_Target_temp(inds,:)=0;
    AM_Target_temp(:,inds)=0;
    AM_Target{i}=LabelBranchesAM(AM_Target_temp);
    
    [~,Dis_Original_voxel_temp] = TraceDistance(AM_Source{i}, r_Source_temp, AM_Target{i}, r_Target_temp,pixelSize,0);
    Dis_Original_voxel_temp=(Dis_Original_voxel_temp*diff([1,inds,size(AM_Source{i},1)])')./sum(diff([1,inds,size(AM_Source{i},1)]));
    TraceDistancesOriginal(i)=Dis_Original_voxel_temp;
    
    b=Optimal_Translation_Transform(Global_Matched_Source,Global_Matched_Target);
%     [SourcePoints_Translation_temp,~]=Perform_Linear_Transform(r_Source_temp,[],[],b);
    SourcePoints_Translation_temp=r_Source_temp+b';
    [~,Distances_Translation_voxels_temp] = TraceDistance(AM_Source{i}, SourcePoints_Translation_temp, AM_Target{i}, r_Target_temp,pixelSize,0);
    Distances_Translation_voxels_temp=(Distances_Translation_voxels_temp*diff([1,inds,size(AM_Source{i},1)])')./sum(diff([1,inds,size(AM_Source{i},1)]));
    TraceDistancesTranslation(i)=Distances_Translation_voxels_temp;
    
    [L,b]=Optimal_Rigid_Transform(Global_Matched_Source,Global_Matched_Target);
    [SourcePoints_Rigid_temp,~]=Perform_Linear_Transform(r_Source_temp,[],L,b);
    [~,Distances_Rigid_voxels_temp] = TraceDistance(AM_Source{i}, SourcePoints_Rigid_temp, AM_Target{i}, r_Target_temp,pixelSize,0);
    Distances_Rigid_voxels_temp=(Distances_Rigid_voxels_temp*diff([1,inds,size(AM_Source{i},1)])')./sum(diff([1,inds,size(AM_Source{i},1)]));
    TraceDistancesRigid(i)=Distances_Rigid_voxels_temp;

end

TraceDistancesNR = [];
TraceDistancesAffine = [];
for nummu = 1:size(mu,2)
    
    mu(nummu)
    [~,LAffine,bAffine]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target,mu(nummu));
    [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu(nummu));
    
    
    for i=1:size(fname_First,1)
        
        [SourcePoints_Affine_temp,~]=Perform_Linear_Transform(SourcePoints{i},[],LAffine,bAffine);
        
        [SourcePoints_NR_temp,~]=Perform_Bspline_Transform(SourcePoints{i}',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
        
        Ncuts=fix(size(AM_Source{i},1)./CutLength);
        inds=ceil(size(AM_Source{i},1)./(Ncuts+1)).*(1:Ncuts);
        
        AM_Source{i}(inds,:)=0;
        AM_Source{i}(:,inds)=0;
        AM_Source{i}=LabelBranchesAM(AM_Source{i});
        
        AM_Target{i}(inds,:)=0;
        AM_Target{i}(:,inds)=0;
        AM_Target{i}=LabelBranchesAM(AM_Target{i});
        
        [~,Distances_Affine_voxels_temp] = TraceDistance(AM_Source{i}, SourcePoints_Affine_temp, AM_Target{i}, TargetPoints{i},pixelSize,0);
        Distances_Affine_voxels_temp=(Distances_Affine_voxels_temp*diff([1,inds,size(AM_Source{i},1)])')./sum(diff([1,inds,size(AM_Source{i},1)]));
        TraceDistancesAffine(nummu,i)=Distances_Affine_voxels_temp;
        
        
        [~,Dis_NonRigid_voxel] = TraceDistance(AM_Source{i}, SourcePoints_NR_temp', AM_Target{i}, TargetPoints{i},pixelSize,0);
        Dis_NonRigid_voxel=(Dis_NonRigid_voxel*diff([1,inds,size(AM_Source{i},1)])')./sum(diff([1,inds,size(AM_Source{i},1)]));
        TraceDistancesNR(nummu,i)=Dis_NonRigid_voxel;
    
    end
    
end

save(['ID_',num2str(ID),'.mat'],'TraceDistancesOriginal','TraceDistancesTranslation','TraceDistancesRigid','TraceDistancesAffine','TraceDistancesNR')


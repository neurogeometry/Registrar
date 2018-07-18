clear all;
clc
close all;
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

pixelSize = [0.26 0.26 0.8];
GTpath = '..\..\RegistrationEvaluation\TimeLaps\';
Functionspath = '..\';
resultpath = '..\..\RegistrationEvaluation\Results-TimeLapse_Holtmaat_StackList\';

addpath([Functionspath,'Functions']);
% load('../data/StackData_TimeLapse_Holtmaat.mat');
load([GTpath,'Matches\DL083-001-Matches.mat']);
showTranceonImage = 0;



% mu = 0:10:1480; %0:20:2000;
mu = 0;

Dis_NonRigid_voxelall = [];
CutLength=100;
% FeaturePositions_NR = load([resultpath,'MatchedPoints_Non-Rigid_mu1.mat']);
FeaturePositions_NR = load('MatchedPoints_Non-Rigid_mu1.mat');
T_Names = {'B','C','D','E','F','G','H','I','J','K','L','M','N'};


% FeaturePositions_NR = load('MatchedPoints_Affine_mu1020.mat'); 1 and 12
% T_Names = {'B','N'};

n= 4;
nxyz = ceil([1024;1024;312]./n);
Boutons = struct2cell(Matches);

% figure(2)
% ylabel({'Distance in Pixels',''});
% xlabel('\mu');
% % xlim([0 max(mu)])
% % ylim([0 5])
% hold on
% drawnow

for ID = 1:12
    sourceID = ID ;
    targetID = sourceID + 1;
    Global_Matched_Source = FeaturePositions_NR.Matched{sourceID,targetID}(:,1:3)';
    Global_Matched_Target = FeaturePositions_NR.Matched{sourceID,targetID}(:,4:6)';
    fname_First = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
    fname_First={fname_First.name}';
    fname_Second = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
    fname_Second={fname_Second.name}';

    
    for i=1:size(fname_First,1)
        sourcePath = [GTpath,'Matches\Traces\',fname_First{i}];%'E:\Datasets\TimeLaps\Matches\Traces\DL083B001-A001.swc';
        targetPath = [GTpath,'Matches\Traces\',fname_Second{i}];%'E:\Datasets\TimeLaps\Matches\Traces\DL083C001-A001.swc';
        
        [AM_Source_temp,r_Source_temp,~]=swc2AM(sourcePath);
        [AM_Target_temp,r_Target_temp,~]=swc2AM(targetPath);
        [AM_Source_temp,r_Source_temp,~] = AdjustPPM(AM_Source_temp,r_Source_temp,zeros(size(r_Source_temp,1),1),ppm);
        [AM_Target_temp,r_Target_temp,~] = AdjustPPM(AM_Target_temp,r_Target_temp,zeros(size(r_Target_temp,1),1),ppm);
        
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
        TraceDistancesOriginal(ID,i)=Dis_Original_voxel_temp;
        
        b=Optimal_Translation_Transform(Global_Matched_Source,Global_Matched_Target);
        [SourcePoints_Translation_temp,~]=Perform_Linear_Transform(r_Source_temp,[],[],b);
        [~,Distances_Translation_voxels_temp] = TraceDistance(AM_Source{i}, SourcePoints_Translation_temp, AM_Target{i}, r_Target_temp,pixelSize,0);
        Distances_Translation_voxels_temp=(Distances_Translation_voxels_temp*diff([1,inds,size(AM_Source{i},1)])')./sum(diff([1,inds,size(AM_Source{i},1)]));
        TraceDistancesTranslation(ID,i)=Distances_Translation_voxels_temp;
        
        [L,b]=Optimal_Rigid_Transform(Global_Matched_Source,Global_Matched_Target);
        [SourcePoints_Rigid_temp,~]=Perform_Linear_Transform(r_Source_temp,[],L,b);
        [~,Distances_Rigid_voxels_temp] = TraceDistance(AM_Source{i}, SourcePoints_Rigid_temp, AM_Target{i}, r_Target_temp,pixelSize,0);
        Distances_Rigid_voxels_temp=(Distances_Rigid_voxels_temp*diff([1,inds,size(AM_Source{i},1)])')./sum(diff([1,inds,size(AM_Source{i},1)]));
        TraceDistancesRigid(ID,i)=Distances_Rigid_voxels_temp;
        
        
        [~,LAffine,bAffine]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target,mu);
        [SourcePoints_Affine_temp,~]=Perform_Linear_Transform(SourcePoints{i},[],LAffine,bAffine);
        [~,Distances_Affine_voxels_temp] = TraceDistance(AM_Source{i}, SourcePoints_Affine_temp, AM_Target{i}, TargetPoints{i},pixelSize,0);
        Distances_Affine_voxels_temp=(Distances_Affine_voxels_temp*diff([1,inds,size(AM_Source{i},1)])')./sum(diff([1,inds,size(AM_Source{i},1)]));
        TraceDistancesAffine(ID,i)=Distances_Affine_voxels_temp;
        Distances_Affine_voxels_temp

        [~,L,b,Cxyz,Nxyz,nxyz,Grid_start]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
        [SourcePoints_NR_temp,~]=Perform_Bspline_Transform(SourcePoints{i}',[],L,b,Cxyz,Nxyz,nxyz,Grid_start,affine);
        [~,Dis_NonRigid_voxel] = TraceDistance(AM_Source{i}, SourcePoints_NR_temp', AM_Target{i}, TargetPoints{i},pixelSize,0);
        Dis_NonRigid_voxel=(Dis_NonRigid_voxel*diff([1,inds,size(AM_Source{i},1)])')./sum(diff([1,inds,size(AM_Source{i},1)]));
        TraceDistancesNR(ID,i)=Dis_NonRigid_voxel;        
       
        
    end
%     plot([mu(1),mu(end)],mean(TraceDistancesOriginal).*[1,1],'r-')
%     plot([mu(1),mu(end)],mean(TraceDistancesTranslation).*[1,1],'m-')
%     plot([mu(1),mu(end)],mean(TraceDistancesRigid).*[1,1],'c-')
%     plot([mu(1),mu(end)],mean(TraceDistancesAffine).*[1,1],'k-')
%     drawnow
    
end
TraceDistancesOriginalAll = TraceDistancesOriginal(:);
TraceDistancesTranslationAll = TraceDistancesTranslation(:);
TraceDistancesRigidAll = TraceDistancesRigid(:);
TraceDistancesAffineAll = TraceDistancesAffine(:);
TraceDistancesNRAll = TraceDistancesNR(:);

figure,hold on
xx=[0.5:1:13.5];
[temp,~]=hist(TraceDistancesOriginalAll,xx,'r-');
plot(xx,temp)
[temp,~]=hist(TraceDistancesTranslationAll,xx,'m-');
plot(xx,temp)
[temp,~]=hist(TraceDistancesRigidAll,xx,'c-');
plot(xx,temp)
[temp,~]=hist(TraceDistancesAffineAll,xx,'b-');
plot(xx,temp)
[temp,~]=hist(TraceDistancesNRAll,xx,'g-');
plot(xx,temp)
axis square, box on
    

clear all
clc
close all
% folder = '10_6_NoAffine/';
folder = '101/New/';
% folder = '10_6_minus5Traces/';
mu = [0,2.^(0:0.5:25)];
% for ID = 1:12
% 
% load([folder,'ID_',num2str(ID),'.mat'])
% 
% figure
% ylabel({'Distance in Pixels',''});
% xlabel('\mu');
% % xlim([0 max(mu)])
% ylim([0 17])
% hold on
% drawnow
% plot([0,log2(mu(2:end))],mean(TraceDistancesAffine,2),'b-')
% plot([0,log2(mu(2:end))],mean(TraceDistancesNR,2),'g-')
% plot([0,log2(mu(end))],mean(TraceDistancesTranslation(:)).*[1,1],'m-')
% plot([0,log2(mu(end))],mean(TraceDistancesRigid(:)).*[1,1],'c-')
% axis square
% box on
% plot([0,log2(mu(end))],(mean(TraceDistancesOriginal(:))).*[1,1],'r-')
% title(['ID = ',num2str(ID)])
% end
% Affine1 = mean(TraceDistancesAffine,1);
% NR1 = mean(TraceDistancesNR,1);
% figure,hold on
% boxplot([TraceDistancesOriginal(:),TraceDistancesTranslation(:),TraceDistancesRigid(:),...
%     Affine1(:),NR1(:)],'Whisker',inf)
% axis square, box on
% title(['ID = ',num2str(ID)])


j = 1;
%  [1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17]
% [7,9,10,11,12,13,14,15,16,17]
% [10,11,12,13,14,15,16,17]
for i = 1:13
    R{i} = load([folder,'ID_',num2str(i),'.mat']);
    Translation(:,j) = R{i}.TraceDistancesTranslation;
    Rigid(:,j) = R{i}.TraceDistancesRigid;
    Before(:,j) = R{i}.TraceDistancesOriginal;
    
    Affine(:,:,j) = R{i}.TraceDistancesAffine;
    NR(:,:,j) = R{i}.TraceDistancesNR;
    
    j = j + 1;
end
TraceDistancesOriginal = mean(Before,2);
TraceDistancesTranslation = mean(Translation,2);
TraceDistancesRigid = mean(Rigid,2);
TraceDistancesAffine = mean(Affine,3);
TraceDistancesNR  = mean(NR,3);

figure(50)
ylabel({'Distance in Pixels',''});
xlabel('\mu');
xlim([0 max(log2(mu(2:end)))])
ylim([2 7])
hold on
drawnow
plot([0,log2(mu(2:end))],mean(TraceDistancesAffine,2),'b-')
plot([0,log2(mu(2:end))],mean(TraceDistancesNR,2),'g-')
plot([0,log2(mu(end))],mean(TraceDistancesTranslation(:)).*[1,1],'m-')
plot([0,log2(mu(end))],mean(TraceDistancesRigid(:)).*[1,1],'c-')
axis square
box on
plot([0,log2(mu(end))],mean(TraceDistancesOriginal(:))-11.*[1,1],'r-')

% for mouse 86
% BestAffine = Affine(20,:,:)
% BestNR = NR(20,:,:)

% for mouse 101
BestAffine = Affine(11,:,:);
BestNR = NR(10,:,:);

BestAffine1 = mean(BestAffine,3);
BestNR1 = mean(BestNR,3);
mean(BestAffine1)
mean(BestNR1)

figure(51),hold on
boxplot([TraceDistancesOriginal(:),TraceDistancesTranslation(:),TraceDistancesRigid(:),...
    BestAffine1(:),BestNR1(:)],'Whisker',inf)
axis square, box on

% affine = mean(TraceDistancesAffine,1);
% NR = mean(NR,1);
% 
% figure(51),hold on
% boxplot([TraceDistancesOriginal(:),TraceDistancesTranslation(:),TraceDistancesRigid(:),...
%     affine(:),NR(:)],'Whisker',inf)
% axis square, box on







% mean_original = mean(TraceDistancesOriginal(:))
% error_original = std(TraceDistancesOriginal(:))/sqrt(size(TraceDistancesOriginal,1))
% 
% mean_translation = mean(TraceDistancesTranslation(:))
% error_translation = std(TraceDistancesTranslation(:))/sqrt(size(TraceDistancesTranslation,1))
% 
% mean_rigid = mean(TraceDistancesRigid(:))
% error_rigid = std(TraceDistancesRigid(:))/sqrt(size(TraceDistancesRigid,1))
% 
% 
% 
% 
% mean_affine = mean(BestAffine1(:))
% error_affine = std(BestAffine1(:))/sqrt(size(BestAffine1,1))
% 
% mean_nr = mean(BestNR1(:))
% error_nr = std(BestNR1(:))/sqrt(size(BestNR1,1))




% affine = mean(TraceDistancesAffine,2);
% mean_affine = min(affine(:))
% error_affine = std(affine(:))/sqrt(size(affine,1))
% 
% NR = mean(TraceDistancesNR,2);
% mean_nr = min(NR(:))
% error_nr = std(NR(:))/sqrt(size(NR,1))






% affine = sort(mean(TraceDistancesAffine,2));
% NR = sort(mean(TraceDistancesNR,2));
% figure(51),hold on
% boxplot([TraceDistancesOriginal(:),TraceDistancesTranslation(:),TraceDistancesRigid(:),...
%     affine(15:36),NR(1:22)],'Whisker',inf)
% axis square, box on

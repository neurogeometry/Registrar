clear all;
close all;
clc;
load('../data/MatchedPoints3_NewFeatures.mat');
load('../data/StackData.mat');
%test
load('../data/StackPositions_Registered.mat');
resolution=StackSizes_mm(1,:)./1000./StackSizes_pixels(1,1:3); % in micromiters
P1 = 271;
P2 = 272;
P3 = 265;
P4 = 266;

%Matched_Points=[];DX_M=[];DY_M=[];DZ_M=[];

Matched_Points = cell2mat(Matched(P1,P2));
for i = 1: size(Matched_Points,1)
        DX_M(i) = Matched_Points(i,4) - Matched_Points(i,1);
        DY_M(i) = Matched_Points(i,5) - Matched_Points(i,2);
        DZ_M(i) = Matched_Points(i,6) - Matched_Points(i,3);
end
P1_P2_DX = mean(DX_M);%866
P1_P2_DY = mean(DY_M);%-10

Matched_Points=[];DX_M=[];DY_M=[];DZ_M=[];
Matched_Points = cell2mat(Matched(P3,P4));
for i = 1: size(Matched_Points,1)
        DX_M(i) = Matched_Points(i,4) - Matched_Points(i,1);
        DY_M(i) = Matched_Points(i,5) - Matched_Points(i,2);
        DZ_M(i) = Matched_Points(i,6) - Matched_Points(i,3);
end
P3_P4_DX = mean(DX_M);%866
P3_P4_DY = mean(DY_M);%-9

Matched_Points=[];DX_M=[];DY_M=[];DZ_M=[];
Matched_Points = cell2mat(Matched(P3,P1));
for i = 1: size(Matched_Points,1)
        DX_M(i) = Matched_Points(i,4) - Matched_Points(i,1);
        DY_M(i) = Matched_Points(i,5) - Matched_Points(i,2);
        DZ_M(i) = Matched_Points(i,6) - Matched_Points(i,3);
end
P3_P1_DX = mean(DX_M);%4
P3_P1_DY = mean(DY_M);%1453




StackPositions_Registered = StackPositions;



StackPositions_Registered(271,:) = StackPositions(271,:);
StackPositions_Registered(272,:) = StackPositions(271,:)+[P1_P2_DX*resolution(1)*1000,P1_P2_DY*resolution(2)*1000,0];
StackPositions_Registered(265,:) = StackPositions(271,:)-[P3_P1_DX*resolution(1)*1000,P3_P1_DY*resolution(2)*1000,0];
StackPositions_Registered(266,:) = StackPositions_Registered(265,:)+[P3_P4_DX*resolution(1)*1000,P3_P4_DY*resolution(2)*1000,0];
save('../data/StackPositions_Registered.mat','StackPositions_Registered');
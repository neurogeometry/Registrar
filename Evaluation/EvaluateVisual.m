close all;
clear all;
clc;
addpath ('NeuronTracerV20');
addpath ('../functions');

% IM_Original=ImportStack('C:\Users\Seyed\Desktop\FusedAll.tif');
% figure,imshow(max(IM_Original,[],3),[0 max(IM_Original(:))]);

% load('../data/StackData.mat','StackPositions','StackSizes_mm','StackSizes_pixels');

% load('../data/StackPositions_Registered_NEW2.mat'); %NEW_ = old global registration
% load('../data/StackPositions_Registered_NEW4.mat'); %NEW_ = old global registration

% load('../data/StackPositions_Registered.mat');
% load('../data/MatchedPoints3_NewFeatures.mat');

% load('../data/MatchedPoints3_NewFeatures');
% load('../data/MatchedPoints3_NewFeatures');

k = 1;
load('../data/StackData_Visual.mat');
load('../data/MatchedPoints_Visual.mat');
load('../data/StackPositions_Registered_Visual.mat');
showImage = 0;
showTraces = 0;
ppm=3;
F_Global = 1;
% useGlobal = 1;
% if useGlobal
%     StackPositions_Registered = StackPositions_Registered_NEW;
% end

% testNum = 6;

% OLD_StackPositions_pixels=StackPositions_pixels; % global stack positions in pixels
% NEW_StackPositions_pixels=StackPositions_Registered_NEW;%./1000./resolution; % global stack positions in pixels

% OLD_StackPositions_pixels(:,1) = max(OLD_StackPositions_pixels(:,1))-OLD_StackPositions_pixels(:,1);
% OLD_StackPositions_pixels(:,2) = max(OLD_StackPositions_pixels(:,2))-OLD_StackPositions_pixels(:,2);
% OLD_StackPositions_pixels(:,3) = max(OLD_StackPositions_pixels(:,2))-OLD_StackPositions_pixels(:,3);

% NEW_StackPositions_pixels(:,1) = max(NEW_StackPositions_pixels(:,1))-NEW_StackPositions_pixels(:,1);
% NEW_StackPositions_pixels(:,2) = max(NEW_StackPositions_pixels(:,2))-NEW_StackPositions_pixels(:,2);
% NEW_StackPositions_pixels(:,3) = max(NEW_StackPositions_pixels(:,2))-NEW_StackPositions_pixels(:,3);


NEW_StackPositions_pixels=StackPositions_Registered;%./1000./resolution; % global stack positions in pixels
NEW_StackPositions_pixels(:,1) = max(NEW_StackPositions_pixels(:,1))-NEW_StackPositions_pixels(:,1);
NEW_StackPositions_pixels(:,2) = max(NEW_StackPositions_pixels(:,2))-NEW_StackPositions_pixels(:,2);
NEW_StackPositions_pixels(:,3) = max(NEW_StackPositions_pixels(:,2))-NEW_StackPositions_pixels(:,3);


for testNum = 1:4
    if testNum == 1
        sourceID = 2;
        targetID = 4;
        sourcePath = '../data/evaluation/Visual/2-4/2_opt1.swc';
        targetPath = '../data/evaluation/Visual/2-4/4_opt1.swc';
        
        dx_GT = 910;
        dy_GT = 66;
        dz_GT = 0;
        
        dx_before = 910;
        dy_before = 66;
        dz_before = 0;
        
        %new
        P1_P2_DX = 910.168831168831;%1722.63636363636;%
        P1_P2_DY = 66.1038961038961;%12;%
        P1_P2_DZ = -2.32467532467532;% 9.09090909090909;
        
        
        % Fiji time = 270930ms - 68515ms --  Globally optimal stitching of tiled 3D microscopic image acquisitions
%         dx_Fiji = 1724;
%         dy_Fiji = 12;
%         dz_Fiji = 10;
%         
%         dx_XUV = 1725;
%         dy_XUV = 17.2222;
%         dz_XUV = 9.5; %  Z Ignored Z = 16.9375
%         
%         dx_Tara = 368;
%         dy_Tara = -6;
%         dz_Tara = 10.1667;
        
        if F_Global %3-4
            
            %1-2
            dx_Fiji = 908.8310546875;
            dy_Fiji = 72.0141220092773;
            dz_Fiji = 0 ; %  Z Ignored Z = 20.13704046153848;
            
            
            dx_XUV = 909.472;
            dy_XUV = 74.6021;
            dz_XUV = 1.93519;
            
            dx_Tara = 910;
            dy_Tara = 66;
            dz_Tara = 2;
            
        end
        
    elseif testNum == 2
        sourceID = 3;
        targetID = 5;
        sourcePath = '../data/evaluation/Visual/3-5/3_opt.swc';
        targetPath = '../data/evaluation/Visual/3-5/5_opt.swc';
        
        
        dx_GT = 909;
        dy_GT = 76;
        dz_GT = 0; % Z Ignored
        
        dx_before = 909;
        dy_before = 76;
        dz_before = 0; % Z Ignored
        
        
        P1_P2_DX = 909.172839506173;%1711.61111111111;
        P1_P2_DY = 75.962962962963;%3.77777777777778;
        P1_P2_DZ =  -1.90123456790123;%-2;
        
%         dx_Fiji = 909;
%         dy_Fiji = 76;
%         dz_Fiji = 0;
%         
%         dx_XUV = 909;
%         dy_XUV = 76;
%         dz_XUV = 0; % Z Ignored Z=-10.625;
%         
%         dx_Tara = 909;
%         dy_Tara = 76;
%         dz_Tara = 0;
        
        if F_Global %3-4
            dx_Fiji = 908.8310546875;
            dy_Fiji = 72.0141220092773;
            dz_Fiji = 0; % Z Ignored Z=-10.2407377608233;
            
            dx_XUV = 909.9721;
            dy_XUV = 67.621;
            dz_XUV = 2.43519 ; %  Z Ignored Z = 19.0417
            
            dx_Tara = 910;
            dy_Tara = 76;
            dz_Tara = 0;
            
        end
        
    elseif testNum == 3
        
        sourceID = 2;
        targetID = 3;
        sourcePath = '../data/evaluation/Visual/3-2/2_opt.swc';
        targetPath = '../data/evaluation/Visual/3-2/3_opt.swc';
        
        
        dx_GT = -76;
        dy_GT = 943;
        dz_GT = 0;
        
        dx_before = -76;
        dy_before = 942.533333333333;
        dz_before = 0;
        
        
        P1_P2_DX = -75.8666666666667;%-16.8157894736842;
        P1_P2_DY = 864.076923076923;%864.105263157895;
        P1_P2_DZ =  -1.73333333333333;%-12.8421052631579;
        
%         dx_Fiji = -76;%
%         dy_Fiji = 943;%
%         dz_Fiji = 0;%
%         
%         dx_XUV = -76;
%         dy_XUV = 943;
%         dz_XUV = 0;
%         
%         dx_Tara = -76;
%         dy_Tara = 943;
%         dz_Tara = 0;
        
        if F_Global  % 2-4
            
            dx_Fiji = -75.9191970825199;%
            dy_Fiji = 949.160644531249;%
            dz_Fiji = 0;%
            
            dx_XUV = -76.102;
            dy_XUV = 948.559;
            dz_XUV = 0.25;
            
            dx_Tara = -76;
            dy_Tara = 943;
            dz_Tara = 0;
        end
        
    elseif testNum == 4
        
        sourceID = 4;
        targetID = 5;
        sourcePath = '../data/evaluation/Visual/5-4/4_opt.swc';
        targetPath = '../data/evaluation/Visual/5-4/5_opt.swc';
        
        
        dx_GT = -77;
        dy_GT = 953;
        dz_GT = 0;
        
        dx_before = -77;
        dy_before = 953;
        dz_before = 0;
        
        %old
        P1_P2_DX = -76.9090909090909;%12;%469.8;
        P1_P2_DY = 952.409090909091;%877;%-14.3;
        P1_P2_DZ = -0.5;%0.666666666666667;%-2.1;
        
        
%         dx_Fiji = -77; %   Failed
%         dy_Fiji = 953;
%         dz_Fiji = 0;
%         
%         dx_XUV = -77; %   Failed
%         dy_XUV = 953;
%         dz_XUV = 0;
%         
%         dx_Tara = -77;
%         dy_Tara = 953;
%         dz_Tara = 0;
        
        if F_Global %2-4
            dx_Fiji = -75.9191970825199;
            dy_Fiji = 949.160644531249;
            dz_Fiji = 0;
            
            dx_XUV = -75.6019;
            dy_XUV = 941.5779;
            dz_XUV = 0.75;
            
            
            
            dx_Tara = -77;
            dy_Tara = 953;
            dz_Tara = 0;
        end
        
        
        
    end
    
    DX_M = [];
    DY_M = [];
    DZ_M = [];
    
    
    %     for i = 1: size(Matched_Points,1)
    %         DX_M(i) = Matched_Points(i,4) - Matched_Points(i,1);
    %         DY_M(i) = Matched_Points(i,5) - Matched_Points(i,2);
    %         DZ_M(i) = Matched_Points(i,6) - Matched_Points(i,3);
    %     end
    %     P1_P2_DX = median(DX_M);%866
    %     P1_P2_DY = median(DY_M);%-10
    %     P1_P2_DZ = median(DZ_M);%-10
    
    
    [AM_G,r_G,R_G]=swc2AM(sourcePath);
    [AM_A,r_A,R_A]=swc2AM(targetPath);
    
    
    [AM_G,r_G,R_G] = AdjustPPM(AM_G,r_G,R_G,ppm);
    [AM_A,r_A,R_A] = AdjustPPM(AM_A,r_A,R_A,ppm);
    
    Trace_AM_G{testNum} = AM_G;
    Trace_r_G{testNum} = r_G;
    Trace_AM_A{testNum} = AM_A;
    Trace_r_A{testNum} = r_A;
    
    r_G_old = r_G;% + [OLD_StackPositions_pixels(sourceID,2),OLD_StackPositions_pixels(sourceID,1),OLD_StackPositions_pixels(sourceID,3)]-1;
    r_A_old = r_A +[dy_before,dx_before,dz_before];% + [OLD_StackPositions_pixels(targetID,2),OLD_StackPositions_pixels(targetID,1),OLD_StackPositions_pixels(targetID,3)]-1;
    
    
    r_G_Global = r_G + [NEW_StackPositions_pixels(sourceID,2),NEW_StackPositions_pixels(sourceID,1),-NEW_StackPositions_pixels(sourceID,3)]-1;
    r_A_Global = r_A + [NEW_StackPositions_pixels(targetID,2),NEW_StackPositions_pixels(targetID,1),-NEW_StackPositions_pixels(targetID,3)]-1;
    
    r_G_Matching = r_G;% + [P1_P2_DX,P1_P2_DY,0];
    r_A_Matching = r_A +[P1_P2_DY,P1_P2_DX,P1_P2_DZ];
    %+ [-10,864,0];
    
    % Fiji
    r_G_Fiji = r_G;% + [P1_P2_DX,P1_P2_DY,0];
    r_A_Fiji = r_A +[dy_Fiji,dx_Fiji,dz_Fiji];
    %+ [-10,864,0];
    
    
    r_G_Tara = r_G;% + [P1_P2_DX,P1_P2_DY,0];
    r_A_Tara = r_A +[dy_Tara,dx_Tara,dz_Tara];
    %+ [-10,864,0];
    
    
    r_G_XUV = r_G;% + [P1_P2_DX,P1_P2_DY,0];
    r_A_XUV = r_A +[dy_XUV,dx_XUV,dz_XUV];
    %+ [-10,864,0];
    
    % Fiji
    r_G_GT = r_G;% + [P1_P2_DX,P1_P2_DY,0];
    r_A_GT = r_A +[dy_GT,dx_GT,dz_GT];
    %+ [-10,864,0];
    
    if showTraces
        figure(2);
        subplot(3,3,1);
        PlotAM(AM_G,r_G_old,'r');
        PlotAM(AM_A,r_A_old,'g');
        title('Before Registration');
        axis equal
        
        hold on;
        subplot(3,3,2);
        PlotAM(AM_G,r_G_Global,'r');
        PlotAM(AM_A,r_A_Global,'g');hold on;
        title('After Registration');
        axis equal
        
        hold on;
        subplot(3,3,3);
        PlotAM(AM_G,r_G_Matching,'r');
        PlotAM(AM_A,r_A_Matching,'g');hold on;
        title('After Matching');
        axis equal
        
        hold on;
        subplot(3,3,4);
        PlotAM(AM_G,r_G_Fiji,'r');
        PlotAM(AM_A,r_A_Fiji,'g');hold on;
        title('Fiji');
        axis equal
        
        hold on;
        subplot(3,3,5);
        PlotAM(AM_G,r_G_Tara,'r');
        PlotAM(AM_A,r_A_Tara,'g');hold on;
        title('Tara Stitcher');
        axis equal
        
        hold on;
        subplot(3,3,6);
        PlotAM(AM_G,r_G_XUV,'r');
        PlotAM(AM_A,r_A_XUV,'g');hold on;
        title('XUV Tools');
        axis equal
        
        hold on;
        subplot(3,3,7);
        PlotAM(AM_G,r_G_GT,'r');
        PlotAM(AM_A,r_A_GT,'g');hold on;
        title('GT');
        axis equal
    end
    
    if F_Global
        [Distances_Proposed, AllDistances_Proposed] = TraceDistance(AM_G, r_G_Global, AM_A, r_A_Global,showTraces);
    else
        [Distances_Proposed, AllDistances_Proposed] = TraceDistance(AM_G, r_G_Matching, AM_A, r_Matching,showTraces);
    end
    
    [Distances_Fiji,AllDistances_Fiji] = TraceDistance(AM_G, r_G_Fiji, AM_A, r_A_Fiji,showTraces);
    [Distances_Tara,AllDistances_Tara] = TraceDistance(AM_G, r_G_Tara, AM_A, r_A_Tara,showTraces);
    [Distances_XUV,AllDistances_XUV] = TraceDistance(AM_G, r_G_XUV, AM_A, r_A_XUV,showTraces);
    [Distances_GT,AllDistances_GT] = TraceDistance(AM_G, r_G_GT, AM_A, r_A_GT,showTraces);
    [Distances_old,AllDistances_old] = TraceDistance(AM_G, r_G_old, AM_A, r_A_old,showTraces);
    
    
    All_Mean_Proposed(testNum,1:2) = [mean(Distances_Proposed),std(AllDistances_Proposed)/sqrt(length(AllDistances_Proposed))]
    All_Mean_Fiji(testNum,1:2) = [mean(Distances_Fiji),std(AllDistances_Fiji)/sqrt(length(AllDistances_Fiji))]
    All_Mean_Tera(testNum,1:2) = [mean(Distances_Tara),std(AllDistances_Tara)/sqrt(length(AllDistances_Tara))]
    All_Mean_XUV(testNum,1:2) = [mean(Distances_XUV),std(AllDistances_XUV)/sqrt(length(AllDistances_XUV))]
    All_Mean_GT(testNum,1:2) = [mean(Distances_GT),std(AllDistances_GT)/sqrt(length(AllDistances_GT))]
    All_Mean_Old(testNum,1:2) = [mean(Distances_old),std(AllDistances_old)/sqrt(length(AllDistances_old))]
    
    
   
end

Proposed = mean(All_Mean_Proposed)
Fiji = mean(All_Mean_Fiji)
Tera = mean(All_Mean_Tera)
XUV = mean(All_Mean_XUV)
GT = mean(All_Mean_GT)
BeforeRegister = mean(All_Mean_Old)

All = [All_Mean_Proposed,All_Mean_Fiji,All_Mean_Tera,All_Mean_XUV,All_Mean_GT,All_Mean_Old]

if showImage
    Dataset = 'Visual';
    
    Image1 = 2;
    Image2 = 4;
    Image3 = 3;
    Image4 = 5;
    
    
    P2_Proposed = [910,66,-2];P3_Proposed = [909,76,-2]; P4_Proposed = [-76,864,-2];P5_Proposed = [-77,943,0];
    stitchedImage_Proposed = StitichAll(Image1,Image2,Image3,Image4,Dataset,P2_Proposed,P3_Proposed,P5_Proposed);
    im3dscroll(stitchedImage_Proposed,[]);caxis([0 30000])
    figure(),imshow(max(stitchedImage_Proposed,[],3),[0 max(stitchedImage_Proposed(:))]);title('Proposed');
    hold on;
    PlotAM(Trace_AM_G{1},Trace_r_G{1}+1,'c');
    PlotAM(Trace_AM_A{1},Trace_r_A{1}+P2_Proposed([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{4},Trace_r_G{4}+1,'r');
    PlotAM(Trace_AM_A{4},Trace_r_A{4}+P4_Proposed([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{4},Trace_r_G{4}+P2_Proposed([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{4},Trace_r_A{4}+P2_Proposed([2,1,3])+P5_Proposed([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{2},Trace_r_G{2}+P4_Proposed([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{2},Trace_r_A{2}+P4_Proposed([2,1,3])+P3_Proposed([2,1,3])+1,'g');
    
    
    P2_Fiji = [819,24,0];P3_Fiji = [908,72,0]; P4_Fiji = [-165,901,0];P5_Fiji = [-76,949,0];
    stitchedImage_Fiji = StitichAll(Image1,Image2,Image3,Image4,Dataset,P2_Proposed,P3_Proposed,P4_Proposed);
    figure(),imshow(max(stitchedImage_Fiji,[],3),[0 max(stitchedImage_Fiji(:))]);title('Fiji');
    hold on;
    PlotAM(Trace_AM_G{1},Trace_r_G{1}+1,'r');
    PlotAM(Trace_AM_A{1},Trace_r_A{1}+P2_Fiji([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{3},Trace_r_G{3}+1,'r');
    PlotAM(Trace_AM_A{3},Trace_r_A{3}+P4_Fiji([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{4},Trace_r_G{4}+P2_Fiji([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{4},Trace_r_A{4}+P2_Fiji([2,1,3])+P5_Fiji([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{2},Trace_r_G{2}+P4_Fiji([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{2},Trace_r_A{2}+P4_Fiji([2,1,3])+P3_Fiji([2,1,3])+1,'g');
    
    P2_XUV = [909,75,2];P3_XUV = [910,68,2]; P4_XUV = [-76,949,0];P5_XUV = [-76,949,1];
    stitchedImage_XUV = StitichAll(Image1,Image2,Image3,Image4,Dataset,P2_XUV,P3_XUV,P4_XUV);
    figure(),imshow(max(stitchedImage_XUV,[],3),[0 max(stitchedImage_XUV(:))]);title('XUV');
    hold on;
    PlotAM(Trace_AM_G{1},Trace_r_G{1}+1,'r');
    PlotAM(Trace_AM_A{1},Trace_r_A{1}+P2_XUV([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{3},Trace_r_G{3}+1,'r');
    PlotAM(Trace_AM_A{3},Trace_r_A{3}+P4_XUV([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{4},Trace_r_G{4}+P2_XUV([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{4},Trace_r_A{4}+P2_XUV([2,1,3])+P5_XUV([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{2},Trace_r_G{2}+P4_XUV([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{2},Trace_r_A{2}+P4_XUV([2,1,3])+P3_XUV([2,1,3])+1,'g');
    
    P2_Tera = [910,66,2];P3_Tera = [910,76,0]; P4_Tera = [-76,943,0];P5_Tera = [-77,953,0];
    stitchedImage_Tera = StitichAll(Image1,Image2,Image3,Image4,Dataset,P2_Tera,P3_Tera,P4_Tera);
    figure(),imshow(max(stitchedImage_Tera,[],3),[0 max(stitchedImage_Tera(:))]);title('Tera');
    hold on;
    PlotAM(Trace_AM_G{1},Trace_r_G{1}+1,'r');
    PlotAM(Trace_AM_A{1},Trace_r_A{1}+P2_Tera([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{3},Trace_r_G{3}+1,'r');
    PlotAM(Trace_AM_A{3},Trace_r_A{3}+P4_Tera([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{4},Trace_r_G{4}+P2_Tera([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{4},Trace_r_A{4}+P2_Tera([2,1,3])+P5_Tera([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{2},Trace_r_G{2}+P4_Tera([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{2},Trace_r_A{2}+P4_Tera([2,1,3])+P3_Tera([2,1,3])+1,'g');
end


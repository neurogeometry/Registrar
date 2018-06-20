close all;
clear all;
clc;
addpath ('NeuronTracerV20');
addpath ('../functions');
% load('../data/StackData.mat','StackPositions','StackSizes_mm','StackSizes_pixels');

% load('../data/StackPositions_Registered_NEW2.mat'); %NEW_ = old global registration
% load('../data/StackPositions_Registered_NEW4.mat'); %NEW_ = old global registration

% load('../data/StackPositions_Registered.mat');
% load('../data/MatchedPoints3_NewFeatures.mat');

% load('../data/MatchedPoints3_NewFeatures');
% load('../data/MatchedPoints3_NewFeatures');

k = 1;
load('../data/StackData_Neuromuscular1.mat');
load('../data/MatchedPoints_Neuromuscular1.mat');
load('../data/StackPositions_Registered_Neuromuscular1.mat');
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


for testNum = 1:1
    if testNum == 1
        sourceID = 106;
        targetID = 107;
        sourcePath = '../data/evaluation/DIADEM_Neuro/106-107/106_opt.swc';
        targetPath = '../data/evaluation/DIADEM_Neuro/106-107/107_opt.swc';
        
        dx_GT = 902.8;
        dy_GT = 2.28;
        dz_GT = 4;
        
        dx_before = 902;
        dy_before = 4;
        dz_before = 0;
        
        %new
        P1_P2_DX = 902.318181818182;%
        P1_P2_DY = 2.29545454545455;%
        P1_P2_DZ =  4;
        
        
        % Fiji time = 270930ms - 68515ms --  Globally optimal stitching of tiled 3D microscopic image acquisitions
        dx_Fiji = 902;
        dy_Fiji = 2;
        dz_Fiji = 4;
        
        dx_XUV = 903;
        dy_XUV = -2;
        dz_XUV = 4; %  Z Ignored Z = 16.9375
        
        dx_Tara = 368;
        dy_Tara = -6;
        dz_Tara = 10.1667;
        
        if F_Global
            
            %1-2
            dx_Fiji = 903.0025325402191;
            dy_Fiji = 2.4733750156627448;
            dz_Fiji = 0 ; %  Z Ignored Z = 20.13704046153848;
            
            dx_XUV = 904.97055;
            dy_XUV = -4.46000000000004;
            dz_XUV = 0 ; %  Z Ignored Z = 19.0417
            
            dx_Tara = 902;
            dy_Tara = 2;
            dz_Tara = 17;
            
        end
        
        P2_Proposed = [902,2];
         
        
    elseif testNum == 2
        sourceID = 114;
        targetID = 115;
        sourcePath = '../data/evaluation/DIADEM_Neuro/114-115/00114_OPT.swc';
        targetPath = '../data/evaluation/DIADEM_Neuro/114-115/00115_OPT.swc';
        
        
        dx_GT = 900.5;
        dy_GT = 3.9;
        dz_GT = -27; % Z Ignored
        
        dx_before = 902;
        dy_before = 4;
        dz_before = -25; % Z Ignored
        
        
        P1_P2_DX = 901.714285714286;
        P1_P2_DY = 2.92857142857143;
        P1_P2_DZ =  -26.6428571428571;
        
        dx_Fiji = 901;
        dy_Fiji = 3;
        dz_Fiji = -11;
        
        dx_XUV = 901;
        dy_XUV = 5;
        dz_XUV = -26; % Z Ignored Z=-10.625;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global
            
            % 3-4
            dx_Fiji = 902.13028612211;
            dy_Fiji = 2.97961839599088;
            dz_Fiji = -20; % Z Ignored Z=-10.2407377608233;
            
            dx_XUV = 900.637;
            dy_XUV = 1.79259;
            dz_XUV = -20; % Z Ignored Z=-9.625;
            
            dx_Tara = 902;
            dy_Tara = 2;
            dz_Tara = -17;
            
        end
        
        P3_Proposed = [901,3]; 
        
    elseif testNum == 3
        
        sourceID = 106;
        targetID = 114;
        sourcePath = '../data/evaluation/DIADEM_Neuro/106_114/00106_opt.swc';
        targetPath = '../data/evaluation/DIADEM_Neuro/106_114/00114_opt.swc';
        
        
        dx_GT = -8.2;
        dy_GT = 920.3;
        dz_GT = 24.2;
        
        dx_before = -8;
        dy_before = 920;
        dz_before = 19;
        
        
        P1_P2_DX = -6.05555555555556;
        P1_P2_DY = 919.111111111111;
        P1_P2_DZ =  24.5;
        
        dx_Fiji = -7;%
        dy_Fiji = 920;%
        dz_Fiji = 36;%
        
        dx_XUV = -8;
        dy_XUV = 919;
        dz_XUV = 33.4167;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global
            dx_Fiji = -6.589166922781885;
            dy_Fiji = 920.2399368957483;
            dz_Fiji = 20; % Z Ignored Z=-36.869147996071675;
            
            dx_XUV = -5.91845;
            dy_XUV = 920.93341;
            dz_XUV = 20; % Z Ignored Z=-36.0833;
            
            dx_Tara = -7;
            dy_Tara = 920;
            dz_Tara = 4;
        end
        
        P4_Proposed = [-6,919];
        
    elseif testNum == 4
        
        sourceID = 107;
        targetID = 115;
        sourcePath = '../data/evaluation/DIADEM_Neuro/107-115/00107_opt.swc';
        targetPath = '../data/evaluation/DIADEM_Neuro/107-115/00115_opt.swc';
        
        
        dx_GT = -7.7;
        dy_GT = 919.1;
        dz_GT = -6.7;
        
        dx_before = -8;
        dy_before = 920;
        dz_before = 0;
        
        %old
        P1_P2_DX = -7.72727272727273;%469.8;
        P1_P2_DY = 919.363636363636;%-14.3;
        P1_P2_DZ = -6.68181818181818;%-2.1;
        
        
        dx_Fiji = -7;
        dy_Fiji = 920;
        dz_Fiji = -4;
        
        dx_XUV = -8;
        dy_XUV = 922.667;
        dz_XUV = -2;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global
            dx_Fiji = -7.46141334089134;
            dy_Fiji = 920.746180276076;
            dz_Fiji = -6.49136977370992;
            
            dx_XUV = -10.2520000000001;
            dy_XUV = 918.266;
            dz_XUV = -7.4166;
            
            dx_Tara = -7;
            dy_Tara = 920;
            dz_Tara = -4;
        end
        
        
        
    end
    
    DX_M = [];
    DY_M = [];
    DZ_M = [];
    
    
    
    
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
        
        hold on;
        subplot(3,3,2);
        PlotAM(AM_G,r_G_Global,'r');
        PlotAM(AM_A,r_A_Global,'g');hold on;
        title('After Registration');
        
        hold on;
        subplot(3,3,3);
        PlotAM(AM_G,r_G_Matching,'r');
        PlotAM(AM_A,r_A_Matching,'g');hold on;
        title('After Matching');
        
        hold on;
        subplot(3,3,4);
        PlotAM(AM_G,r_G_Fiji,'r');
        PlotAM(AM_A,r_A_Fiji,'g');hold on;
        title('Fiji');
        
        hold on;
        subplot(3,3,5);
        PlotAM(AM_G,r_G_Tara,'r');
        PlotAM(AM_A,r_A_Tara,'g');hold on;
        title('Tara Stitcher');
        
        hold on;
        subplot(3,3,6);
        PlotAM(AM_G,r_G_XUV,'r');
        PlotAM(AM_A,r_A_XUV,'g');hold on;
        title('XUV Tools');
        
        hold on;
        subplot(3,3,7);
        PlotAM(AM_G,r_G_GT,'r');
        PlotAM(AM_A,r_A_GT,'g');hold on;
        title('GT');
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
    Dataset = 'Neuromuscular1';
     
    Image1 = 106;
    Image2 = 107;
    Image3 = 114;
    Image4 = 115;
    
    
    P2_Proposed = [902,2,4];P3_Proposed = [901,3,-27]; P4_Proposed = [-6,919,24]; P5_Proposed = [-8,919,-7];
    stitchedImage_Proposed = StitichAll(Image1,Image2,Image3,Image4,Dataset,P2_Proposed,P3_Proposed,P4_Proposed);
    figure(),imshow(max(stitchedImage_Proposed,[],3),[0 max(stitchedImage_Proposed(:))]);title('Proposed');
    hold on;
    PlotAM(Trace_AM_G{1},Trace_r_G{1}+1,'r');
    PlotAM(Trace_AM_A{1},Trace_r_A{1}+P2_Proposed([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{3},Trace_r_G{3}+1,'r');
    PlotAM(Trace_AM_A{3},Trace_r_A{3}+P4_Proposed([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{4},Trace_r_G{4}+P2_Proposed([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{4},Trace_r_A{4}+P2_Proposed([2,1,3])+P5_Proposed([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{2},Trace_r_G{2}+P4_Proposed([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{2},Trace_r_A{2}+P4_Proposed([2,1,3])+P3_Proposed([2,1,3])+1,'g');

    
    P2_Fiji = [903,2,0];P3_Fiji = [902,3,-20]; P4_Fiji = [-6,920,20];P5_Fiji = [-7,921,-6];
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
    
    P2_XUV = [904,-4,0];P3_XUV = [900,2,-20]; P4_XUV = [-6,921,20];P5_XUV = [-10,918,-7];
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
    
    P2_Tera = [902,2,17];P3_Tera = [902,2,-17]; P4_Tera = [-7,920,4];P5_Tera = [-7,920,-4];
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


close all;
clear all;
clc;
addpath ('NeuronTracerV20');
addpath ('../functions');

k = 1;
load('../data/StackData_Holtmaat_ch0.mat');
load('../data/MatchedPoints_Holtmaat_ch0.mat');
load('../data/StackPositions_Registered_Holtmaat_ch0.mat');
showImage = 0;
showTraces = 0;
ppm=1;
F_Global = 1;


NEW_StackPositions_pixels=StackPositions_Registered;%./1000./resolution; % global stack positions in pixels
NEW_StackPositions_pixels(:,1) = max(NEW_StackPositions_pixels(:,1))-NEW_StackPositions_pixels(:,1);
NEW_StackPositions_pixels(:,2) = max(NEW_StackPositions_pixels(:,2))-NEW_StackPositions_pixels(:,2);
NEW_StackPositions_pixels(:,3) = max(NEW_StackPositions_pixels(:,2))-NEW_StackPositions_pixels(:,3);


for testNum = 4:4
    if testNum == 1
        sourceID = 8;
        targetID = 1;
        sourcePath = '../data/evaluation/Holtmaat/8-1/8_opt.swc';
        targetPath = '../data/evaluation/Holtmaat/8-1/1_opt.swc';
        
        dx_GT = 1722.76923;
        dy_GT = 10;
        dz_GT = 18;
        
        dx_before = 1730.76923;
        dy_before = 0;
        dz_before = 0.25;
        
        %new
        P1_P2_DX = 1728.75;%1722.63636363636;%
        P1_P2_DY = 13.25;%12;%
        P1_P2_DZ = 8.875;% 9.09090909090909;
        
        
        % Fiji time = 270930ms - 68515ms --  Globally optimal stitching of tiled 3D microscopic image acquisitions
        dx_Fiji = 1724;
        dy_Fiji = 12;
        dz_Fiji = 10;
        
        dx_XUV = 1725;
        dy_XUV = 17.2222;
        dz_XUV = 9.5; %  Z Ignored Z = 16.9375
        
        dx_Tara = 368;
        dy_Tara = -6;
        dz_Tara = 10.1667;
        
        if F_Global
            
            %1-2
            dx_Fiji = 1722;
            dy_Fiji = 12;
            dz_Fiji = 10 ; %  Z Ignored Z = 20.13704046153848;
            
            dx_XUV = 1719;
            dy_XUV = 6;
            dz_XUV = 3.25 ; %  Z Ignored Z = 19.0417
            
            dx_Tara = 1809;
            dy_Tara = 4;
            dz_Tara = -2;
            
        end
        
    elseif testNum == 2
        sourceID = 7;
        targetID = 6;
        sourcePath = '../data/evaluation/Holtmaat/7-6/7_opt.swc';
        targetPath = '../data/evaluation/Holtmaat/7-6/6_opt.swc';
        
        
        dx_GT = 1730.76923;
        dy_GT = 0;
        dz_GT = 3.5; % Z Ignored
        
        dx_before = 1730.76923;
        dy_before = 0;
        dz_before = 3.5; % Z Ignored
        
        
        P1_P2_DX = 1711.625;%1711.61111111111;
        P1_P2_DY = 3.75;%3.77777777777778;
        P1_P2_DZ =  -1.9375;%-2;
        
        dx_Fiji = 1714;
        dy_Fiji = 4;
        dz_Fiji = -2;
        
        dx_XUV = 1710.67;
        dy_XUV = 4.11111;
        dz_XUV = -9.5; % Z Ignored Z=-10.625;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global
        dx_Fiji = 1707.354259790076;
        dy_Fiji = -5.34170391603038;
        dz_Fiji = -6.683407832060783; % Z Ignored Z=-10.2407377608233;
        
        dx_XUV = 1715.5;
        dy_XUV = -3;
        dz_XUV = -3; 
        
        dx_Tara = 1808;
        dy_Tara = 12;
        dz_Tara = 11;
        
        end
        
    elseif testNum == 3
        
        sourceID = 7;
        targetID = 8;
        sourcePath = '../data/evaluation/Holtmaat/7-8/7_opt.swc';
        targetPath = '../data/evaluation/Holtmaat/7-8/8_opt.swc';
        
        
        dx_GT = 0;
        dy_GT = 865.385000000009;
        dz_GT = 3;
        
        dx_before = 0;
        dy_before = 865.385000000009;
        dz_before = 3;
        
        
        P1_P2_DX = -16.8717948717949;%-16.8157894736842;
        P1_P2_DY = 864.076923076923;%864.105263157895;
        P1_P2_DZ =  -12.8974358974359;%-12.8421052631579;
        
        dx_Fiji = -17;%  
        dy_Fiji = 864;%
        dz_Fiji = -13;%
        
        dx_XUV = -17;
        dy_XUV = 859.556;
        dz_XUV = -16;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global 
            
        dx_Fiji = -15.979837737814364;
        dy_Fiji = 863.5919350951258;
        dz_Fiji = -13.816129809748475; % Z Ignored Z=-36.869147996071675;
        
        dx_XUV = -17.5;
        dy_XUV = 861;
        dz_XUV = -1; % Z Ignored Z=-36.0833;
        
        dx_Tara = -12;
        dy_Tara = 904;
        dz_Tara = -14;
        end
        
    elseif testNum == 4
        
        sourceID = 6;
        targetID = 1;
        sourcePath = '../data/evaluation/Holtmaat/6-1/00006_opt.swc';
        targetPath = '../data/evaluation/Holtmaat/6-1/00001_opt.swc';
        
        
        dx_GT = 0;
        dy_GT = 865;
        dz_GT = 0.75;
        
        dx_before = 0;
        dy_before = 865;
        dz_before = 0.75;
        
        %old
        P1_P2_DX = 10;%12;%469.8;
        P1_P2_DY = 876;%877;%-14.3;
        P1_P2_DZ = 0.3333333333333;%0.666666666666667;%-2.1;
        
        
        dx_Fiji = -7; %   Failed
        dy_Fiji = 920;
        dz_Fiji = -4;
        
        dx_XUV = -8; %   Failed
        dy_XUV = 922.667;
        dz_XUV = -2;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global
            dx_Fiji = 1.33409752789021;
            dy_Fiji = 880.933639011156;
            dz_Fiji = -2.86727802231231;
            
            dx_XUV = 14;
            dy_XUV = 870;
            dz_XUV = -1.25;
            
            dx_Tara = 12;
            dy_Tara = 904;
            dz_Tara = 1;
        end
        
    
        
    end
    
    DX_M = [];
    DY_M = [];
    DZ_M = [];
    
    
    
    [AM_G,r_G,R_G]=swc2AM(sourcePath);
    [AM_A,r_A,R_A]=swc2AM(targetPath);
    
    
    [AM_G,r_G,R_G] = AdjustPPM(AM_G,r_G,R_G,ppm);
    [AM_A,r_A,R_A] = AdjustPPM(AM_A,r_A,R_A,ppm);
    
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
%     temp = sourceID;
%     sourceID = targetID;
%     targetID = temp;
Matched_Points = cell2mat(Matched(sourceID,targetID));
    if size(Matched_Points,1)>0
        tifFile_1 = dir([char(StackList(sourceID,2)) '/*.tif' ]);
        Source_Stack_Folder = [tifFile_1(1).folder,'/'];
        Source_Stack_File = tifFile_1(1).name;
        disp(['Reading ',Source_Stack_File]);
        IM_Source=ImportStack([Source_Stack_Folder,Source_Stack_File]);
        
        tifFile_2 = dir([char(StackList(targetID,2)) '/*.tif' ]);
        Target_Stack_Folder = [tifFile_2(1).folder,'/'];
        Target_Stack_File = tifFile_2(1).name;
        disp(['Reading ',Target_Stack_File]);
        IM_Target=ImportStack([Target_Stack_Folder,Target_Stack_File]);
        
        ShowPairedImages(IM_Source,IM_Target,sourceID,targetID,Matched_Points);
        
    end
end


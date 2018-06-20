close all;
clear all;
clc;
addpath ('NeuronTracerV20');
addpath ('../functions');


k = 1;

load('../data/StackData_DIADEM.mat');
load('../data/MatchedPoints_DIADEM.mat');
load('../data/StackPositions_Registered_DIADEM.mat');

showImage = 0;
showTraces = 0;
ppm=3;
F_Global = 1;

% useGlobal = 1;
% if useGlobal
%     StackPositions_Registered = StackPositions_Registered_NEW;
% end

% testNum = 6;

NEW_StackPositions_pixels=StackPositions_Registered;%./1000./resolution; % global stack positions in pixels

% NEW_StackPositions_pixels(:,1) = max(NEW_StackPositions_pixels(:,1))-NEW_StackPositions_pixels(:,1);
% NEW_StackPositions_pixels(:,2) = max(NEW_StackPositions_pixels(:,2))-NEW_StackPositions_pixels(:,2);
% NEW_StackPositions_pixels(:,3) = max(NEW_StackPositions_pixels(:,2))-NEW_StackPositions_pixels(:,3);

for testNum = 2:2
    if testNum == 1
        sourceID = 1;
        targetID = 2;
        sourcePath = '../data/evaluation/DIADEM/1-2/00001_opt.swc';
        targetPath = '../data/evaluation/DIADEM/1-2/00002_opt.swc';
        
        dx_GT = 453.5;
        dy_GT = -23;
        dz_GT = 19.5;
        
        dx_before = 453;
        dy_before = -23;
        dz_before = 16;
        
        %new
        P1_P2_DX = 454;%
        P1_P2_DY = -22;%
        P1_P2_DZ =  20;
        
        
        % Fiji time = 270930ms - 68515ms --  Globally optimal stitching of tiled 3D microscopic image acquisitions
        dx_Fiji = 453;
        dy_Fiji = -23;
        dz_Fiji = 19;
        
        dx_XUV = 452;
        dy_XUV = -24;
        dz_XUV = 19; % Z is ignored, otherwise is failed. dz_XUV = 5;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global % # 5 ,6 in Fiji
            dx_Fiji = 453; %Failed , before registration replaced
            dy_Fiji = -23; %Failed , before registration replaced
            dz_Fiji =  16; %Failed , before registration replaced
            
            dx_Tara = 453;
            dy_Tara = -23;
            dz_Tara = 2;
            
            %             dx_XUV = 861;
            %             dy_XUV = -10.5;
            %             dz_XUV = 1.125;
        end
        
    elseif testNum == 2
        sourceID = 2;
        targetID = 3;
        sourcePath = '../data/evaluation/DIADEM_Neo1/2-3/00002_opt.swc';
        targetPath = '../data/evaluation/DIADEM_Neo1/2-3/00003_opt.swc';
        pixelSize = [0.3 0.3 1.02];
        
        dx_GT = 427.5;
        dy_GT = -22.6;
        dz_GT = -31.5;
        
        dx_before = 426;
        dy_before = -22;
        dz_before = -32;
        
        
        P1_P2_DX = 427;
        P1_P2_DY = -23;
        P1_P2_DZ =  -30.5;
        
        dx_Fiji = 426;
        dy_Fiji = -24;
        dz_Fiji = -32;
        
        dx_XUV = 425;
        dy_XUV = -24;
        dz_XUV = -26.6667;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global % # 5,4 in Fiji
            dx_Fiji = 427;
            dy_Fiji = -21;
            dz_Fiji =  -23;
            
            dx_Tara = 453;
            dy_Tara = -23;
            dz_Tara = -2;
            
            %             dx_XUV = 861;
            %             dy_XUV = -10.5;
            %             dz_XUV = 1.125;
        end
        
    elseif testNum == 3
        
        sourceID = 5;
        targetID = 4;
        sourcePath = '../data/evaluation/DIADEM/5-4/00005_opt.swc';
        targetPath = '../data/evaluation/DIADEM/5-4/00004_opt.swc';
        
        
        dx_GT = 465.2;
        dy_GT = 14.4;
        dz_GT = -16.7;
        
        dx_before = 456;
        dy_before = 17;
        dz_before = -18;
        
        
        P1_P2_DX = 467.5;
        P1_P2_DY = 16.5;
        P1_P2_DZ =  -17;
        
        dx_Fiji = 456;% 399; %Failed
        dy_Fiji = 17;% -346;
        dz_Fiji = -16;% 1;
        
        dx_XUV = 456;
        dy_XUV = 17;
        dz_XUV = -18;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global % 2,3 in Fiji
            dx_Fiji = 455.999999999999 ;
            dy_Fiji = 5.99999999999946;
            dz_Fiji =  -13;
            
            
            dx_Tara = 450;
            dy_Tara = 18;
            dz_Tara = -2;
            
            
            %             dx_XUV = 861;
            %             dy_XUV = -10.5;
            %             dz_XUV = 1.125;
        end
        
        
        
        
    elseif testNum == 4
        
        sourceID = 6;
        targetID = 5;
        sourcePath = '../data/evaluation/DIADEM/6-5/6.swc';
        targetPath = '../data/evaluation/DIADEM/6-5/00005.swc';
        
        
        dx_GT = 469.8;
        dy_GT = -14.56;
        dz_GT = -1.6;
        
        dx_before = 468;
        dy_before = -14;
        dz_before = -1;
        
        %old
        P1_P2_DX = 470.333333333333;%469.8;
        P1_P2_DY = -14.5555555555556;%-14.3;
        P1_P2_DZ = -1.66666666666667;%-2.1;
        
        
        dx_Fiji = 471;
        dy_Fiji = -14;
        dz_Fiji = -2;
        
        dx_XUV = 470;
        dy_XUV = -14;
        dz_XUV = -0.333333;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global % # 1,2 in Fiji
            dx_Fiji = 470.999999999999;
            dy_Fiji = -14.00000000000071;
            dz_Fiji = -7.064278826347438E-15;
            
            dx_Tara = 474;
            dy_Tara = 0;
            dz_Tara = -2;
            
            %             dx_XUV = 861;
            %             dy_XUV = -10.5;
            %             dz_XUV = 1.125;
        end
        
    elseif testNum == 5
        
        sourceID = 4;
        targetID = 3;
        sourcePath = '../data/evaluation/DIADEM/4-3/4_opt.swc';
        targetPath = '../data/evaluation/DIADEM/4-3/3_opt.swc';
        
        
        dx_GT = 27;
        dy_GT = 478;
        dz_GT = 1;
        
        dx_before = 28;
        dy_before = 459;
        dz_before = -2;
        
        P1_P2_DX = 25.3333333333333;%469.8;
        P1_P2_DY = 466.25;%-14.3;
        P1_P2_DZ = -0.416666666666667;%-2.1;
        
        % not yet
        dx_Fiji = 28;% Failed = 211;
        dy_Fiji = 459;%266;
        dz_Fiji = -5;%-7;
        
        dx_XUV = 25;
        dy_XUV = 466.667;
        dz_XUV = -4;
        
        % not yet
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        
        if F_Global
            dx_Fiji = 455.689189910888;
            dy_Fiji = -5.64257812499961;
            dz_Fiji =  12.7318849563599;
            
            %             dx_XUV = 861;
            %             dy_XUV = -10.5;
            %             dz_XUV = 1.125;
        end
        
    end
    
    DX_M = [];
    DY_M = [];
    DZ_M = [];
    
    %     Matched_Points = cell2mat(Matched(targetID,sourceID));
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
    
    r_G_old = r_G;% + [OLD_StackPositions_pixels(sourceID,2),OLD_StackPositions_pixels(sourceID,1),OLD_StackPositions_pixels(sourceID,3)]-1;
    r_A_old = r_A +[dy_before,dx_before,dz_before];% + [OLD_StackPositions_pixels(targetID,2),OLD_StackPositions_pixels(targetID,1),OLD_StackPositions_pixels(targetID,3)]-1;
    
    
    r_G_Global = r_G + [NEW_StackPositions_pixels(targetID,2),NEW_StackPositions_pixels(targetID,1),-NEW_StackPositions_pixels(targetID,3)]-1;
    r_A_Global = r_A + [NEW_StackPositions_pixels(sourceID,2),NEW_StackPositions_pixels(sourceID,1),-NEW_StackPositions_pixels(sourceID,3)]-1;
    
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
        [Distances_Proposed, AllDistances_Proposed] = TraceDistance(AM_G, r_G_Global, AM_A, r_A_Global,pixelSize,showTraces);
    else
        [Distances_Proposed, AllDistances_Proposed] = TraceDistance(AM_G, r_G_Matching, AM_A, r_Matching,pixelSize,showTraces);
    end
    
    [Distances_Fiji,AllDistances_Fiji] = TraceDistance(AM_G, r_G_Fiji, AM_A, r_A_Fiji,pixelSize,showTraces);
    [Distances_Tara,AllDistances_Tara] = TraceDistance(AM_G, r_G_Tara, AM_A, r_A_Tara,pixelSize,showTraces);
    [Distances_XUV,AllDistances_XUV] = TraceDistance(AM_G, r_G_XUV, AM_A, r_A_XUV,pixelSize,showTraces);
    [Distances_GT,AllDistances_GT] = TraceDistance(AM_G, r_G_GT, AM_A, r_A_GT,pixelSize,showTraces);
    [Distances_old,AllDistances_old] = TraceDistance(AM_G, r_G_old, AM_A, r_A_old,pixelSize,showTraces);
    
    
    All_Mean_Proposed(testNum,1:2) = [mean(Distances_Proposed),std(AllDistances_Proposed)/sqrt(length(AllDistances_Proposed))]
    All_Mean_Fiji(testNum,1:2) = [mean(Distances_Fiji),std(AllDistances_Fiji)/sqrt(length(AllDistances_Fiji))]
    All_Mean_Tera(testNum,1:2) = [mean(Distances_Tara),std(AllDistances_Tara)/sqrt(length(AllDistances_Tara))]
    All_Mean_XUV(testNum,1:2) = [mean(Distances_XUV),std(AllDistances_XUV)/sqrt(length(AllDistances_XUV))]
    All_Mean_GT(testNum,1:2) = [mean(Distances_GT),std(AllDistances_GT)/sqrt(length(AllDistances_GT))]
    All_Mean_Old(testNum,1:2) = [mean(Distances_old),std(AllDistances_old)/sqrt(length(AllDistances_old))]
    
    
    
    
    %     mean(Distances_old)
    %     mean(Distances_Proposed)
    %     mean(Distances_Fiji)
    %     mean(Distances_Tara)
    %     mean(Distances_XUV)
    %     mean(Distances_GT)
    %
    %        k
    %     RMean(1,k:k+2)= Distances_old;
    %     RMean(2,k:k+2)= Distances_Proposed;
    %     RMean(3,k:k+2)= Distances_Fiji;
    %     RMean(4,k:k+2)= Distances_Tara;
    %     RMean(5,k:k+2)= Distances_XUV;
    %     RMean(6,k:k+2)= Distances_GT;
    
    
    
    
    %     h_length=7;
    %     IndexStr_old = CompareTraces_byLength(AM_G, r_G_old, AM_A, r_A_old, h_length);
    %     IndexStr_Global = CompareTraces_byLength(AM_G, r_G_Global, AM_A, r_A_Global, h_length);
    %     IndexStr_Matching = CompareTraces_byLength(AM_G, r_G_Matching, AM_A, r_A_Matching, h_length);
    %     IndexStr_Fiji = CompareTraces_byLength(AM_G, r_G_Fiji, AM_A, r_A_Fiji, h_length);
    %     IndexStr_Tara = CompareTraces_byLength(AM_G, r_G_Tara, AM_A, r_A_Tara, h_length);
    %     IndexStr_XUV = CompareTraces_byLength(AM_G, r_G_XUV, AM_A, r_A_XUV, h_length);
    %     IndexStr_GT = CompareTraces_byLength(AM_G, r_G_GT, AM_A, r_A_GT, h_length);
    %
    %
    %
    %     Dis = 10;
    %     IndexStr_old = IndexStr_old.Dag_trace_full(IndexStr_old.Dag_trace_full<Dis);
    %     IndexStr_Global = IndexStr_Global.Dag_trace_full(IndexStr_Global.Dag_trace_full<Dis);
    %     IndexStr_Matching = IndexStr_Matching.Dag_trace_full(IndexStr_Matching.Dag_trace_full<Dis);
    %     IndexStr_Fiji = IndexStr_Fiji.Dag_trace_full(IndexStr_Fiji.Dag_trace_full<Dis);
    %     IndexStr_Tara = IndexStr_Tara.Dag_trace_full(IndexStr_Tara.Dag_trace_full<Dis);
    %     IndexStr_XUV = IndexStr_XUV.Dag_trace_full(IndexStr_XUV.Dag_trace_full<Dis);
    %     IndexStr_GT = IndexStr_GT.Dag_trace_full(IndexStr_GT.Dag_trace_full<Dis);
    %
    %     PlotLimit = 500;
    %     figure;
    %     plot(sort(IndexStr_old(1:min(size(IndexStr_old,2),PlotLimit))));
    %     hold on
    %     plot(sort(IndexStr_Matching(1:min(size(IndexStr_Matching,2),PlotLimit))));
    %     hold on
    %     plot(sort(IndexStr_Global(1:min(size(IndexStr_Global,2),PlotLimit))));
    %     hold on
    %     plot(sort(IndexStr_Fiji(1:min(size(IndexStr_Fiji,2),PlotLimit))));
    %     hold on
    %     plot(sort(IndexStr_Tara(1:min(size(IndexStr_Tara,2),PlotLimit))));
    %     hold on
    %     plot(sort(IndexStr_XUV(1:min(size(IndexStr_XUV,2),PlotLimit))));
    %     %     hold on
    %     %     plot(sort(IndexStr_GT(1:PlotLimit)));
    %     legend('Distances before registraion','Distances After Matching','Distances After Global Registration','Fiji','TeraStitcher','XUVTools','Ground Truth');
    %
    %     figure;
    %     bin_centers=0.2:0.4:10;
    %     IndexStr_oldD_hist=hist(IndexStr_old,bin_centers);
    %     IndexStr_MatchingD_hist=hist(IndexStr_Matching,bin_centers);
    %     IndexStr_GlobalD_hist=hist(IndexStr_Global,bin_centers);
    %     IndexStr_FijiD_hist=hist(IndexStr_Fiji,bin_centers);
    %     IndexStr_Tara_hist=hist(IndexStr_Tara,bin_centers);
    %     IndexStr_XUV_hist=hist(IndexStr_XUV,bin_centers);
    %     IndexStr_GT_hist=hist(IndexStr_GT,bin_centers);
    %     plot(bin_centers,IndexStr_oldD_hist);
    %     hold on
    %     plot(bin_centers,IndexStr_MatchingD_hist);
    %     hold on
    %     plot(bin_centers,IndexStr_GlobalD_hist);
    %     hold on
    %     plot(bin_centers,IndexStr_FijiD_hist);
    %     hold on
    %     plot(bin_centers,IndexStr_Tara_hist);
    %     hold on
    %     plot(bin_centers,IndexStr_XUV_hist);
    %     hold on
    %     plot(bin_centers,IndexStr_GT_hist);
    %     legend('Distances before registraion','Distances After Matching','Distances After Global Registration','fiji','Tera','GT');
    %
    %     GroundTruth=[mean(IndexStr_GT),std(IndexStr_GT)/sqrt(length(IndexStr_GT))];
    %     before_registraion = [mean(IndexStr_old)-mean(IndexStr_GT),std(IndexStr_old)/sqrt(length(IndexStr_old))];
    %     after_matching=[mean(IndexStr_Matching)-mean(IndexStr_GT),std(IndexStr_Matching)/sqrt(length(IndexStr_Matching))];
    %     after_global=[mean(IndexStr_Global)-mean(IndexStr_GT),std(IndexStr_Global)/sqrt(length(IndexStr_Global))];
    %     Fiji = [mean(IndexStr_Fiji)-mean(IndexStr_GT),std(IndexStr_Fiji)/sqrt(length(IndexStr_Fiji))];
    %     Tara=[mean(IndexStr_Tara)-mean(IndexStr_GT),std(IndexStr_Tara)/sqrt(length(IndexStr_Tara))];
    %     XUV=[mean(IndexStr_XUV)-mean(IndexStr_GT),std(IndexStr_XUV)/sqrt(length(IndexStr_XUV))];
    %
    %     before_registraion = [mean(IndexStr_old),std(IndexStr_old)/sqrt(length(IndexStr_old))]
    %
    %     if F_Global
    %         after_global=[mean(IndexStr_Global),std(IndexStr_Global)/sqrt(length(IndexStr_Global))]
    %     else
    %         after_matching=[mean(IndexStr_Matching),std(IndexStr_Matching)/sqrt(length(IndexStr_Matching))]
    %     end
    %     Fiji = [mean(IndexStr_Fiji),std(IndexStr_Fiji)/sqrt(length(IndexStr_Fiji))]
    %     Tara=[mean(IndexStr_Tara),std(IndexStr_Tara)/sqrt(length(IndexStr_Tara))]
    %     XUV=[mean(IndexStr_XUV),std(IndexStr_XUV)/sqrt(length(IndexStr_XUV))]
    %     GroundTruth=[mean(IndexStr_GT),std(IndexStr_GT)/sqrt(length(IndexStr_GT))]
    %
    %     k
    %     RMean(1,k:k+1)= before_registraion;
    %     RMean(2,k:k+1)= after_matching;
    %     RMean(3,k:k+1)= after_global;
    %     RMean(4,k:k+1)= Fiji;
    %     RMean(5,k:k+1)= Tara;
    %     RMean(6,k:k+1)= XUV;
    %     RMean(7,k:k+1)= GroundTruth;
    %     k = k + 2;
    %     RMean;
end
Proposed = mean(All_Mean_Proposed)
Fiji = mean(All_Mean_Fiji)
Tera = mean(All_Mean_Tera)
XUV = mean(All_Mean_XUV)
GT = mean(All_Mean_GT)
BeforeRegister = mean(All_Mean_Old)

All = [All_Mean_Proposed,All_Mean_Fiji,All_Mean_Tera,All_Mean_XUV,All_Mean_GT,All_Mean_Old]


if showImage
    temp = sourceID;
    sourceID = targetID;
    targetID = temp;
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


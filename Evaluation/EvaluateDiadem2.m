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
load('../data/StackData_DIADEM2.mat');
load('../data/MatchedPoints_DIADEM2_.mat');
load('../data/StackPositions_Registered_DIADEM2_.mat');
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


for testNum = 3:6
    if testNum == 1
        sourceID = 6;
        targetID = 5;
        sourcePath = '../data/evaluation/DIADEM2/6-5/6_opt.swc';
        targetPath = '../data/evaluation/DIADEM2/6-5/5_opt.swc';
        
        dx_GT = 368;
        dy_GT = -5.35;
        dz_GT = 14.9;
        
        dx_before = 368;
        dy_before = -7;
        dz_before = 14;
        
        %new
        P1_P2_DX = 367.943661971831;%
        P1_P2_DY = -5.80281690140845;%
        P1_P2_DZ =  14.830985915493;
        
        
        % Fiji time = 270930ms - 68515ms --  Globally optimal stitching of tiled 3D microscopic image acquisitions
        dx_Fiji = 368;
        dy_Fiji = -6;
        dz_Fiji = 15;
        
        dx_XUV = 368;
        dy_XUV = -6;
        dz_XUV = 10.1667; % Z is ignored, otherwise is failed. dz_XUV = 5;
        
        %         if F_Global
        %             dx_Fiji = 867;
        %             dy_Fiji = -9;
        %             dz_Fiji = 0;
        %
        %             dx_XUV = 862;
        %             dy_XUV = -9.5;
        %             dz_XUV = 1.625;
        %         end
        
        dx_Tara = 368;
        dy_Tara = -6;
        dz_Tara = 10.1667;
        
    elseif testNum == 2
        sourceID = 5;
        targetID = 4;
        sourcePath = '../data/evaluation/DIADEM2/5-4/5_OPT.swc';
        targetPath = '../data/evaluation/DIADEM2/5-4/4_OPT.swc';
        
        
        dx_GT = 445;
        dy_GT = -10;
        dz_GT = 14;
        
        dx_before = 445;
        dy_before = -10;
        dz_before = 14;
        
        
        P1_P2_DX = 445.024390243902;
        P1_P2_DY = -7.95121951219512;
        P1_P2_DZ =  13.219512195122;
        
        dx_Fiji = 445;
        dy_Fiji = -8;
        dz_Fiji = 13;
        
        dx_XUV = 446;
        dy_XUV = -9;
        dz_XUV = 12.5333;
        
        %         if F_Global
        %             dx_Fiji = 867;
        %             dy_Fiji = -9;
        %             dz_Fiji = 0;
        %
        %             dx_XUV = 862;
        %             dy_XUV = -9.5;
        %             dz_XUV = 1.625;
        %         end
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
    elseif testNum == 3
        
        sourceID = 4;
        targetID = 3;
        sourcePath = '../data/evaluation/DIADEM2/4-3/4_opt.swc';
        targetPath = '../data/evaluation/DIADEM2/4-3/3_opt.swc';
        
        
        dx_GT = 470;
        dy_GT = -8;
        dz_GT = -7;
        
        dx_before = 470;
        dy_before = -8;
        dz_before = -7;
        
        
        P1_P2_DX = 471.590909090909;
        P1_P2_DY = -7.95454545454545;
        P1_P2_DZ =  -6.90909090909091;
        
        dx_Fiji = 472;% 399; %Failed
        dy_Fiji = -8;% -346;
        dz_Fiji = -7;% 1;
        
        dx_XUV = 471;
        dy_XUV = -7;
        dz_XUV = -2.66667;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global %6,7 in Fiji
            dx_Fiji = 472.000000001735;
            dy_Fiji = -7.83841816861599;
            dz_Fiji = -0.726282805207044;
            
            %             dx_XUV = 862;
            %             dy_XUV = -9.5;
            %             dz_XUV = 1.625;
            
            
            
            dx_Tara = 470;
            dy_Tara = -8;
            dz_Tara = 4;
            
            
        end
        
        
        
    elseif testNum == 4
        
        sourceID = 1;
        targetID = 2;
        sourcePath = '../data/evaluation/DIADEM2/1-2/1_opt.swc';
        targetPath = '../data/evaluation/DIADEM2/1-2/2_opt.swc';
        
        
        dx_GT = 5.5;
        dy_GT = 443.3;
        dz_GT = 15;
        
        dx_before = 6;
        dy_before = 443;
        dz_before = 15;
        
        %old
        P1_P2_DX = 5.11111111111111;%469.8;
        P1_P2_DY = 443.259259259259;%-14.3;
        P1_P2_DZ = 15.047619047619;%-2.1;
        
        
        dx_Fiji = 5;
        dy_Fiji = 443;
        dz_Fiji = 15;
        
        dx_XUV = 6;
        dy_XUV = 443;
        dz_XUV = 11.5;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global %4,8 in Fiji
            dx_Fiji = 4.99999999871716;
            dy_Fiji = 443.999999999774;
            dz_Fiji = 7.00000000007902;
            
            %             dx_XUV = 862;
            %             dy_XUV = -9.5;
            %             dz_XUV = 1.625;
            dx_Tara = 6;
            dy_Tara = 475;
            dz_Tara = 7;
            
        end
        
        
    elseif testNum == 5
        
        sourceID = 9;
        targetID = 1;
        sourcePath = '../data/evaluation/DIADEM2/9-1/9_opt.swc';
        targetPath = '../data/evaluation/DIADEM2/9-1/1_opt.swc';
        
        
        dx_GT = 237;
        dy_GT = -67;
        dz_GT = -24.2;
        
        dx_before = 233;
        dy_before = -62;
        dz_before = -29;
        
        
        P1_P2_DX = 237;%469.8;
        P1_P2_DY = -66.85;%-14.3;
        P1_P2_DZ = -24.8625;%-2.1;
        
        % not yet
        dx_Fiji = 237;% Failed = 211;
        dy_Fiji = -67;%266;
        dz_Fiji = -24;%-7;
        
        dx_XUV = 238;
        dy_XUV = -67;
        dz_XUV = -24.4;
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        if F_Global %3,4 in Fiji
            dx_Fiji = 236.999999894512;
            dy_Fiji = -66.6317058234223;
            dz_Fiji = -9.10488256668778;
            
            %             dx_XUV = 862;
            %             dy_XUV = -9.5;
            %             dz_XUV = 1.625;
            dx_Tara = 450;
            dy_Tara = -21;
            dz_Tara = -5;
        end
        
        % not yet
        
        
        
    elseif testNum == 6
        
        sourceID = 8;
        targetID = 9;
        sourcePath = '../data/evaluation/DIADEM2/8-9/8_opt.swc';
        targetPath = '../data/evaluation/DIADEM2/8-9/9_opt.swc';
        
        
        dx_GT = 471;
        dy_GT = -8.6;
        dz_GT = 27;
        
        dx_before = 473;
        dy_before = -8;
        dz_before = 27;
        
        
        P1_P2_DX = 470.25;%469.8;
        P1_P2_DY = -8.5;%-14.3;
        P1_P2_DZ = 26.8333333333333;%-2.1;
        
        % not yet
        dx_Fiji = 471;% Failed = 211;
        dy_Fiji = -8;%266;
        dz_Fiji = 27;%-7;
        
        dx_XUV = 470;
        dy_XUV = -9;
        dz_XUV = 27;% Z ignored 13.1667
        
        
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        
        if F_Global %2,3 in Fiji
            dx_Fiji = 469.568842082536 ;
            dy_Fiji = -7.392560132593;
            dz_Fiji = 3.00000000009108;
            
            %             dx_XUV = 862;
            %             dy_XUV = -9.5;
            %             dz_XUV = 1.625;
            dx_Tara = 475;
            dy_Tara = -6;
            dz_Tara = 7;
            
        end
        
        % not yet
        
        
    elseif testNum == 7
        
        sourceID = 7;
        targetID = 8;
        sourcePath = '../data/evaluation/DIADEM2/7-8/7_opt.swc';
        targetPath = '../data/evaluation/DIADEM2/7-8/8_opt.swc';
        
        
        dx_GT = 438;
        dy_GT = -11.7;
        dz_GT = -39;
        
        dx_before = 435;
        dy_before = -12;
        dz_before = -37;
        
        
        P1_P2_DX = 438.117647058824;%469.8;
        P1_P2_DY = -11.6470588235294;%-14.3;
        P1_P2_DZ = -38.9411764705882;%-2.1;
        
        % not yet
        dx_Fiji = 438;% Failed = 211;
        dy_Fiji = -12;%266;
        dz_Fiji = -38;%-7;
        
        dx_XUV = 438.613;
        dy_XUV = -12.9467;
        dz_XUV = 19.7;% Z ignored 13.1667
        
        % not yet
        dx_Tara = 454;
        dy_Tara = -8;
        dz_Tara = 19;
        
        
    elseif testNum == 8
        
        sourceID = 10;
        targetID = 9;
        sourcePath = '../data/evaluation/DIADEM2/9-10/10_opt.swc';
        targetPath = '../data/evaluation/DIADEM2/9-10/9_opt.swc';
        
        
        dx_GT = 217;
        dy_GT = 276;
        dz_GT = 2;
        
        dx_before = 217;
        dy_before = 276;
        dz_before = 2;
        
        
        P1_P2_DX = 213.727272727273;%469.8;
        P1_P2_DY = 279.045454545455;%-14.3;
        P1_P2_DZ = 2.13636363636364;%-2.1;
        
        % not yet
        dx_Fiji = 214;% Failed = 211;
        dy_Fiji = 279;%266;
        dz_Fiji = 2;%-7;
        
        dx_XUV = 214;
        dy_XUV = 279;
        dz_XUV = 2.93333;% Z ignored 13.1667
        
        % not yet
        dx_Tara = 214;
        dy_Tara = 279;
        dz_Tara = 12.93333;
        
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
    % %         hold on
    % %         plot(sort(IndexStr_GT(1:PlotLimit)));
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
    %     if F_Global
    %         after_global=[mean(IndexStr_Global),std(IndexStr_Global)/sqrt(length(IndexStr_Global))]
    %     else
    %         after_matching=[mean(IndexStr_Matching),std(IndexStr_Matching)/sqrt(length(IndexStr_Matching))]
    %     end
    %
    %     Fiji = [mean(IndexStr_Fiji),std(IndexStr_Fiji)/sqrt(length(IndexStr_Fiji))]
    %     Tara=[mean(IndexStr_Tara),std(IndexStr_Tara)/sqrt(length(IndexStr_Tara))]
    %     XUV=[mean(IndexStr_XUV),std(IndexStr_XUV)/sqrt(length(IndexStr_XUV))]
    %     GroundTruth=[mean(IndexStr_GT),std(IndexStr_GT)/sqrt(length(IndexStr_GT))]
    
    
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


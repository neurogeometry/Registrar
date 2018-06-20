close all;
clear all;
clc;
addpath ('NeuronTracerV20');
addpath ('../functions');

do_Rigid = 0;
do_Affine = 0;
do_nonrigid = 0;

% load('../data/StackData.mat','StackPositions','StackSizes_mm','StackSizes_pixels');

load('../../GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList/MatchedPoints_Translation.mat');
Matched_Translation = Matched;

if do_Rigid
load('../../GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList/MatchedPoints_Rigid.mat');
Matched_Rigid = Matched;
end

if do_Affine
load('../../GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList/MatchedPoints_Affine.mat');
Matched_Affine = Matched;
end
if do_nonrigid
load('../../GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList/MatchedPoints_Non-Rigid.mat');
Matched_Non_Rigid = Matched;
end

% load('../data/StackPositions_Registered_NEW2.mat'); %NEW_ = old global registration
% load('../data/StackPositions_Registered.mat'); %NEW_ = old global registration
load('../../GUI\NCT_Registration\MicroscopeFiles\Results-MouseLight_StackList/StackPositions_Registered.mat');


 
% load('../data/StackPositions_Registered.mat');
% load('../data/MatchedPoints3_NewFeatures.mat');

% load('../data/MatchedPoints3_NewFeatures');
load('../data/MatchedPoints');

k = 1;
load('../data/StackData.mat');
showImage = 0;
showTraces = 0;
ppm=3;
F_Global = 1;
% useGlobal = 1;
% if useGlobal
%     StackPositions_Registered = StackPositions_Registered_NEW;
% end

% testNum = 6;

pixelSize = [0.377607421875 0.277486979166667 0.99601593625498];
for testNum = 3:3
    
    if testNum == 1
        %         sourceID = 272;%
        %         targetID = 271;%
        sourceID = 4;%
        targetID = 3;%
        
        matches = Matched_Translation{targetID,sourceID}; 
        Proposed_Translationdt = mean(matches(:,1:3)-matches(:,4:6));
        
        if do_Rigid
        matches = Matched_Rigid{targetID,sourceID}; 
        Proposed_Rigiddt = mean(matches(:,1:3)-matches(:,4:6));
        end
        if do_Affine
        matches = Matched_Affine{targetID,sourceID}; 
        Proposed_Affinedt = mean(matches(:,1:3)-matches(:,4:6));
        end
        if do_nonrigid
        matches = Matched_Non_Rigid{targetID,sourceID}; 
        Proposed_NonRigiddt = mean(matches(:,1:3)-matches(:,4:6));
        end
        sourcePath = '../data/evaluation/3-759-760/00760_opt.swc';%00760_Opt.swc';%
        targetPath = '../data/evaluation/3-759-760/00759_opt.swc';%00759_Opt.swc';%
        
        % Affine 
%         b = [870.01269389456,-9.98965658776069,-1.2764822418988];
%         L = [1.02973779822887,-0.00682260741690419,-0.000531048778935273;...
%       -0.00583676478044099,1.00071801306223,0.00234666815299084;...
%        0.00601755826814938,-3.05168471925705e-05,0.998545168416539];
   
        b = [869.38844665471344796970, -8.72560301850990072126, -1.56615136884207117873];
        L = [1.00106649436125416663, -0.00387553564739811195, 0.00108013562941820969;-0.00096061363126749176, 0.99938179206001598320, 0.00128218295923291986;0.00215744862187259266, 0.00034991872081132616, 0.99905636163179734055];

        H = [0.13292809148806428166, -0.03012522873112585819, -0.07447015678899361613;0.00831927735635549365, 1.00106221498589920138, 0.00070624540981573006;0.00217421362852968959, -0.00072958730143453160, 0.99963073483302700151];
%    H = [0.0588963669298282,0.023574001503752,-0.061565823213923;...
%         0.0123712072404139,0.998224738612834,0.000225281841702108;...
%        0.00368982628084379,-0.00188375372777412,1.00028918357206];
%    
%    H = [0.0352664215990992,0.0362871149137769,0.00984783148367658;...
%        0.00895077542327651,1.00136554782449,-0.00447880831514867;...
%        0.00109158122139807,-2.42861043571437e-05,0.999992822087483];
        
        dx_GT = 865;
        dy_GT = -9.8;
        dz_GT = 1;
        
        % 3  00760 00759
        P1_P2_DX = 864.97;% 864.756756756757;
        P1_P2_DY = -9.34;%-9.37837837837838;%-10
        P1_P2_DZ = 1.12;%1.05405405405405;%-10
        
        dx_Fiji = 870;%pairwise
        dy_Fiji = -9;
        dz_Fiji = 1;
        
        
        dx_Tara = 859;
        dy_Tara = -9;
        dz_Tara = 1;
        
        
        dx_XUV = 432*2; %Pairwise
        dy_XUV = 5*2;
        dz_XUV = 1*2;
        
        if F_Global %1,2 Fiji
            dx_Fiji = 867;
            dy_Fiji = -9;
            dz_Fiji = 0;
            
            dx_XUV = 862;
            dy_XUV = -9.5;
            dz_XUV = 1.625;
            
            dx_Tara = 859;
            dy_Tara = -9;
            dz_Tara = 1;
        end
        
    elseif testNum == 2
        %         sourceID = 272;%272;%
        %         targetID = 266;%271;%
        sourceID = 4;%272;%
        targetID = 2;%271;%
        
        matches = Matched_Translation{targetID,sourceID}; 
        Proposed_Translationdt = mean(matches(:,1:3)-matches(:,4:6));
        if do_Rigid
        matches = Matched_Rigid{targetID,sourceID}; 
        Proposed_Rigiddt = mean(matches(:,1:3)-matches(:,4:6));
        end
        if do_Affine
        matches = Matched_Affine{targetID,sourceID}; 
        Proposed_Affinedt = mean(matches(:,1:3)-matches(:,4:6));
        end
        if do_nonrigid
        matches = Matched_Non_Rigid{targetID,sourceID}; 
        Proposed_NonRigiddt = mean(matches(:,1:3)-matches(:,4:6));
        end
        %     sourcePath = '../data/evaluation/271_272_265_266/00760_Opt.swc';%00760_Opt.swc';%
        %     targetPath = '../data/evaluation/271_272_265_266/00736_Opt.swc';%00759_Opt.swc';%
        sourcePath = '../data/evaluation/4-266-272/00760_opt.swc';%00760_Opt.swc';%
        targetPath = '../data/evaluation/4-266-272/00736_opt.swc';%00759_Opt.swc';%
        
        % Affine
        b = [2.42849382542079,1452.44259512299,0.999999999999972];
        L = [1.00287905066583,0.00238110175766704,-0.0014342897726266;...
            0.000764426572516798,1.01492888571262,-0.00233351121320761;...
            2.81519041975262e-17,2.97770238213684e-16,1];
        H = [0.99660685343146726556, -0.00147272604409197936, 0.00184968625008551883;-0.00475086195550167116, 0.02757086960094735570, 0.03630372339343584942;0.00108371068609024288, -0.00127383167607467260, 1.00015588391790544165];
        dy_GT = 1453;
        dx_GT = 4;
        dz_GT = -1;
        
        % 4 00736 00760
        P1_P2_DX = 4.17307692307692;%4.14583333333333;
        P1_P2_DY = 1453.03846153846;%1452.71875;%-10
        P1_P2_DZ = -0.961538461538462;%-0.96875;%-10
        
        dy_Fiji = 1452; %pairwise 31.9536332029551         0.118097446159476
        dx_Fiji = 4;
        dz_Fiji = -1;
        
        
        dy_Tara = 1435;%214310
        dx_Tara = 4;
        dz_Tara = -1;
        
        dy_XUV = 1452; %Pairwise
        dx_XUV = 3;
        dz_XUV = 1;
        
        if F_Global
            dy_Fiji = 1382; %Faile
            dx_Fiji = 3;
            dz_Fiji = 0;
            
            dy_XUV = 1453.5;
            dx_XUV = 3;
            dz_XUV = 0.875;
            
            dy_Tara = 1453;
            dx_Tara = 3;
            dz_Tara = 1;
        end
        
        %     ppm=4;
    elseif testNum == 3
        sourceID = 3;%271;%272;%
        targetID = 1;%265;%271;%
        
       matches = Matched_Translation{targetID,sourceID}; 
        Proposed_Translationdt = mean(matches(:,1:3)-matches(:,4:6));
        if do_Rigid
        matches = Matched_Rigid{targetID,sourceID}; 
        Proposed_Rigiddt = mean(matches(:,1:3)-matches(:,4:6));
        end
        if do_Affine
        matches = Matched_Affine{targetID,sourceID}; 
        Proposed_Affinedt = mean(matches(:,1:3)-matches(:,4:6));
        end
        if do_nonrigid
        matches = Matched_Non_Rigid{targetID,sourceID}; 
        Proposed_NonRigiddt = mean(matches(:,1:3)-matches(:,4:6));
        end
        %     sourcePath = '../data/evaluation/271_272_265_266/00760_Opt.swc';%00760_Opt.swc';%
        %     targetPath = '../data/evaluation/271_272_265_266/00736_Opt.swc';%00759_Opt.swc';%
        sourcePath = '../data/evaluation/MouseLight/5-735-759/00759_opt.swc';%00760_Opt.swc';%
        targetPath = '../data/evaluation/MouseLight/5-735-759/00735_opt.swc';%00759_Opt.swc';%
        
        % Affine
        b = [1.78807501569065,1453.17889548856,1.2760208311022];
        L = [1.00328702640887,-0.000782328906235278,0.00177522496503123;...
             0.000936677285862418,0.977610964922664,0.00140197229619826;...
            -0.000175504110150646,-0.00267832312914809,0.998708138615279];
        H = [0.99724860423516503705, -0.00074737948434199484, -0.00874199477500223893;0.00167964348910853035, 0.02794429684578999054, -0.02146993241819282125;0.00000156439581075919, -0.00066966705548163869, 0.99998889903804977219];
        
        dy_GT = 1452;
        dx_GT = 3;
        dz_GT = -1;
        
        % 5  00759  00735
        P1_P2_DX = 3.20547945205479;%3.13;
        P1_P2_DY = 1452.80821917808;%1453.06;%-10
        P1_P2_DZ = -1.0958904109589;%-0.95;%-10
        
        dy_Fiji = 1453;%pairwise
        dx_Fiji = 3;
        dz_Fiji = -1;
        
        
        dy_Tara = 1446;
        dx_Tara = 5;
        dz_Tara = 0;
        
        dy_XUV = 1455;
        dx_XUV = 3;
        dz_XUV = -2;
        if F_Global
            dy_Fiji = 1382;
            dx_Fiji = 9;
            dz_Fiji = 0;
            
            dy_XUV = 1452.5;
            dx_XUV = 2;
            dz_XUV = -1.375;
            
            dy_Tara = 1453;
            dx_Tara = 3;
            dz_Tara = 1;
            
        end
        
        %     ppm=4;
    elseif testNum == 4
        sourceID = 2;%266;%
        targetID = 1;%265;%
        
        matches = Matched_Translation{targetID,sourceID}; 
        Proposed_Translationdt = mean(matches(:,1:3)-matches(:,4:6));
        if do_Rigid
        matches = Matched_Rigid{targetID,sourceID}; 
        Proposed_Rigiddt = mean(matches(:,1:3)-matches(:,4:6));
        end
        if do_Affine
        matches = Matched_Affine{targetID,sourceID}; 
        Proposed_Affinedt = mean(matches(:,1:3)-matches(:,4:6));
        end
        if do_nonrigid
        matches = Matched_Non_Rigid{targetID,sourceID}; 
        Proposed_NonRigiddt = mean(matches(:,1:3)-matches(:,4:6));
        end
        sourcePath = '../data/evaluation/6-735-736/00736.swc';%00760_Opt.swc';%
        targetPath = '../data/evaluation/6-735-736/00735.swc';%00759_Opt.swc';%
        
        % Affine
        b = [855.598910569938,-11.963701008333,-0.985026014654011];
        L = [1.02094651702822,0.0146025503207548,0.0040186136236321;...
       0.0126374152632418,1.00244501258628,0.00498693480638159;...
        0.00239144897828268,0.000581299592459287,0.998633446807071];
        H = [0.07920950585785151155, 0.01221697600785342716, 0.00817004332509204713;0.01119508736276470309, 0.99900413453487757476, -0.00244983425217287169;0.00106852201250063326, -0.00002245616030709077, 1.00000250869501217110];
        dx_GT = 867;
        dy_GT = -12;
        dz_GT = 1;
        
        % 6 00735 00736
        P1_P2_DX = 861.717948717949;%861;
        P1_P2_DY = -9.48717948717949;%-9.7;%-10
        P1_P2_DZ = 1.28205128205128;%1.1;%-10
        
        dx_Fiji = 870;%Pairwise %Failed
        dy_Fiji = -9;
        dz_Fiji = 1;
        
        
        dx_Tara = 859;
        dy_Tara = -9;
        dz_Tara = 1;
        
        dx_XUV = 865;
        dy_XUV = 9;
        dz_XUV = 1;
        if F_Global
            %         dx_Fiji = 801;
            %         dy_Fiji = -9;
            %         dz_Fiji = 1;
            dx_Fiji = 850;
            dy_Fiji = -9;
            dz_Fiji = 1;
            
            dx_XUV = 864.919440349651; %Failed
            dy_XUV = 0;
            dz_XUV = 0;
            
            dx_Tara = 859;
            dy_Tara = -9;
            dz_Tara = 1;
        end
    end
    
    
    
    
    OLD_StackPositions_pixels=StackPositions_pixels; % global stack positions in pixels
    NEW_StackPositions_pixels=StackPositions_Registered;%./1000./resolution; % global stack positions in pixels
    
    OLD_StackPositions_pixels(:,1) = max(OLD_StackPositions_pixels(:,1))-OLD_StackPositions_pixels(:,1);
    OLD_StackPositions_pixels(:,2) = max(OLD_StackPositions_pixels(:,2))-OLD_StackPositions_pixels(:,2);
    OLD_StackPositions_pixels(:,3) = max(OLD_StackPositions_pixels(:,2))-OLD_StackPositions_pixels(:,3);
    
    NEW_StackPositions_pixels(:,1) = max(NEW_StackPositions_pixels(:,1))-NEW_StackPositions_pixels(:,1);
    NEW_StackPositions_pixels(:,2) = max(NEW_StackPositions_pixels(:,2))-NEW_StackPositions_pixels(:,2);
    NEW_StackPositions_pixels(:,3) = max(NEW_StackPositions_pixels(:,3))-NEW_StackPositions_pixels(:,3);
    
    DX_M = [];
    DY_M = [];
    DZ_M = [];
    
    Matched_Points = cell2mat(Matched(targetID,sourceID));
    for i = 1: size(Matched_Points,1)
        DX_M(i) = Matched_Points(i,4) - Matched_Points(i,1);
        DY_M(i) = Matched_Points(i,5) - Matched_Points(i,2);
        DZ_M(i) = Matched_Points(i,6) - Matched_Points(i,3);
    end
    if testNum == 1 ||  testNum == 2
        P1_P2_DX = median(DX_M);%866
        P1_P2_DY = median(DY_M);%-10
        P1_P2_DZ = median(DZ_M);%-10
    end
    
    
    [AM_G,r_G,R_G]=swc2AM(sourcePath);
    [AM_A,r_A,R_A]=swc2AM(targetPath);
    
    
    [AM_G,r_G,R_G] = AdjustPPM(AM_G,r_G,R_G,ppm);
    [AM_A,r_A,R_A] = AdjustPPM(AM_A,r_A,R_A,ppm);
    
    Trace_AM_G{testNum} = AM_G;
    Trace_r_G{testNum} = r_G;
    Trace_AM_A{testNum} = AM_A;
    Trace_r_A{testNum} = r_A;
    
    r_G_Affine = r_G(:,[2,1,3])';%L*r_G'+b'*ones(1,size(r_G',2));
    b = b';
    r_A_Affine = L'*r_A(:,[2,1,3])'+b*ones(1,size(r_A(:,[2,1,3])',2));

     r_A_Projective = r_A(:,[2,1,3])';
     r_G_Projective = H*r_G(:,[2,1,3])';
    
    r_G_old = r_G + [OLD_StackPositions_pixels(sourceID,2),OLD_StackPositions_pixels(sourceID,1),OLD_StackPositions_pixels(sourceID,3)]-1;
    r_A_old = r_A + [OLD_StackPositions_pixels(targetID,2),OLD_StackPositions_pixels(targetID,1),OLD_StackPositions_pixels(targetID,3)]-1;
    
    
    r_G_Global = r_G + [NEW_StackPositions_pixels(sourceID,2),NEW_StackPositions_pixels(sourceID,1),NEW_StackPositions_pixels(sourceID,3)]-1;
    r_A_Global = r_A  +[NEW_StackPositions_pixels(targetID,2),NEW_StackPositions_pixels(targetID,1),NEW_StackPositions_pixels(targetID,3)]-1;
    
    r_G_Matching = r_G;% + [P1_P2_DX,P1_P2_DY,0];
    r_A_Matching = r_A +[P1_P2_DY,P1_P2_DX,P1_P2_DZ];
    %+ [-10,864,0];
    
    
    
    % Proposed Translation
    r_G_Proposed_Translation = r_G;% + [P1_P2_DX,P1_P2_DY,0];
    r_A_Proposed_Translation = r_A +[Proposed_Translationdt(2),Proposed_Translationdt(1),Proposed_Translationdt(3)];
    %+ [-10,864,0];
    
    if do_Rigid
        % Proposed Rigid
    r_G_Proposed_Rigid = r_G;% + [P1_P2_DX,P1_P2_DY,0];
    r_A_Proposed_Rigid = r_A +[Proposed_Rigiddt(2),Proposed_Rigiddt(1),Proposed_Rigiddt(3)];
    %+ [-10,864,0];
    end
    if do_Affine
        % Proposed Affine
    r_G_Proposed_Affine = r_G;% + [P1_P2_DX,P1_P2_DY,0];
    r_A_Proposed_Affine = r_A +[Proposed_Affinedt(2),Proposed_Affinedt(1),Proposed_Affinedt(3)];
    %+ [-10,864,0];
    end
    if do_nonrigid
        % Proposed NonRigid
    r_G_Proposed_NonRigid = r_G;% + [P1_P2_DX,P1_P2_DY,0];
    r_A_Proposed_NonRigid = r_A +[Proposed_NonRigiddt(2),Proposed_NonRigiddt(1),Proposed_NonRigiddt(3)];
    %+ [-10,864,0];
    end
    
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
    
    
    [Distances_Affine,AllDistances_Affine] = TraceDistance(AM_G, r_G_Affine', AM_A, r_A_Affine',pixelSize,showTraces);
     mean(Distances_Affine);
    
     [Distances_Projective,AllDistances_Projective] = TraceDistance(AM_G, r_G_Projective', AM_A, r_A_Projective',pixelSize,showTraces);
     mean(Distances_Projective);
    if F_Global
        [Distances_Proposed, AllDistances_Proposed] = TraceDistance(AM_G, r_G_Global, AM_A, r_A_Global,pixelSize,showTraces);
    else
        [Distances_Proposed, AllDistances_Proposed] = TraceDistance(AM_G, r_G_Matching, AM_A, r_Matching,pixelSize,showTraces);
    end
%     mean(Distances_Proposed)

    [Distances_Proposed_Translation,AllDistances_Proposed_Translation] = TraceDistance(AM_G, r_G_Proposed_Translation, AM_A, r_A_Proposed_Translation,pixelSize,showTraces);
    if do_Rigid
    [Distances_Proposed_Rigid,AllDistances_Proposed_Rigid] = TraceDistance(AM_G, r_G_Proposed_Rigid, AM_A, r_A_Proposed_Rigid,pixelSize,showTraces);
    end
    if do_Affine
    [Distances_Proposed_Affine,AllDistances_Proposed_Affine] = TraceDistance(AM_G, r_G_Proposed_Affine, AM_A, r_A_Proposed_Affine,pixelSize,showTraces);
    end
    if do_nonrigid
    [Distances_Proposed_NonRigid,AllDistances_Proposed_NonRigid] = TraceDistance(AM_G, r_G_Proposed_NonRigid, AM_A, r_A_Proposed_NonRigid,pixelSize,showTraces);
    end
    
    [Distances_Fiji,AllDistances_Fiji] = TraceDistance(AM_G, r_G_Fiji, AM_A, r_A_Fiji,pixelSize,showTraces);
    [Distances_Tara,AllDistances_Tara] = TraceDistance(AM_G, r_G_Tara, AM_A, r_A_Tara,pixelSize,showTraces);
    [Distances_XUV,AllDistances_XUV] = TraceDistance(AM_G, r_G_XUV, AM_A, r_A_XUV,pixelSize,showTraces);
    [Distances_GT,AllDistances_GT] = TraceDistance(AM_G, r_G_GT, AM_A, r_A_GT,pixelSize,showTraces);
    [Distances_old,AllDistances_old] = TraceDistance(AM_G, r_G_old, AM_A, r_A_old,pixelSize,showTraces);
    
    All_Mean_Proposed_Rigid = 0;
    All_Mean_Proposed_Affine = 0;
    All_Mean_Proposed_NonRigid = 0;
    
    All_Mean_Proposed_Translation(testNum,1:2) = [mean(Distances_Proposed_Translation),std(AllDistances_Proposed_Translation)/sqrt(length(AllDistances_Proposed_Translation))]
    
    if do_Rigid
    All_Mean_Proposed_Rigid(testNum,1:2) = [mean(Distances_Proposed_Rigid),std(AllDistances_Proposed_Rigid)/sqrt(length(AllDistances_Proposed_Rigid))]
    end
    if do_Affine
    All_Mean_Proposed_Affine(testNum,1:2) = [mean(Distances_Proposed_Affine),std(AllDistances_Proposed_Affine)/sqrt(length(AllDistances_Proposed_Affine))]
    end
    if do_nonrigid 
        All_Mean_Proposed_NonRigid(testNum,1:2) = [mean(Distances_Proposed_NonRigid),std(AllDistances_Proposed_NonRigid)/sqrt(length(AllDistances_Proposed_NonRigid))]
    end
    
    
    All_Mean_Proposed(testNum,1:2) = [mean(Distances_Proposed),std(AllDistances_Proposed)/sqrt(length(AllDistances_Proposed))]
    All_Mean_Fiji(testNum,1:2) = [mean(Distances_Fiji),std(AllDistances_Fiji)/sqrt(length(AllDistances_Fiji))]
    All_Mean_Tera(testNum,1:2) = [mean(Distances_Tara),std(AllDistances_Tara)/sqrt(length(AllDistances_Tara))]
    All_Mean_XUV(testNum,1:2) = [mean(Distances_XUV),std(AllDistances_XUV)/sqrt(length(AllDistances_XUV))]
    All_Mean_GT(testNum,1:2) = [mean(Distances_GT),std(AllDistances_GT)/sqrt(length(AllDistances_GT))]
    All_Mean_Old(testNum,1:2) = [mean(Distances_old),std(AllDistances_old)/sqrt(length(AllDistances_old))]
    All_Mean_Affine(testNum,1:2) = [mean(Distances_Affine),std(AllDistances_Affine)/sqrt(length(AllDistances_Affine))]
    All_Mean_Projective(testNum,1:2) = [mean(Distances_Projective),std(Distances_Projective)/sqrt(length(Distances_Projective))]
   
    AllDistances_Proposed = sort(AllDistances_Proposed);
    AllDistances_Fiji = sort(AllDistances_Fiji);
    AllDistances_Tara = sort(AllDistances_Tara);
    AllDistances_Proposed = sort(AllDistances_XUV);
    
   figure();plot(AllDistances_Proposed(1:500),'g');hold on;
   plot(AllDistances_Fiji(1:500),'r');hold on;
   plot(AllDistances_Tara(1:500),'c');hold on;
   plot(AllDistances_XUV(1:500),'b');hold on;

   
end

Proposed = mean(All_Mean_Proposed)
Fiji = mean(All_Mean_Fiji)
Tera = mean(All_Mean_Tera)
XUV = mean(All_Mean_XUV)
GT = mean(All_Mean_GT)
BeforeRegister = mean(All_Mean_Old)
Affine = mean(All_Mean_Affine)

All = [All_Mean_Proposed_Translation,All_Mean_Proposed,All_Mean_Fiji,All_Mean_Tera,All_Mean_XUV,All_Mean_GT,All_Mean_Old,All_Mean_Affine,All_Mean_Projective]
% All = [All_Mean_Proposed_Translation,All_Mean_Proposed_Rigid,All_Mean_Proposed_Affine,All_Mean_Proposed_NonRigid,All_Mean_Proposed,All_Mean_Fiji,All_Mean_Tera,All_Mean_XUV,All_Mean_GT,All_Mean_Old,All_Mean_Affine,All_Mean_Projective]

if showImage
    load('Evaluate.mat');
    Dataset = 'Mine';
    Image1 = 4;
    Image2 = 3;
    Image3 = 2;
    Image4 = 1;
    
        P3_Proposed = round(NEW_StackPositions_pixels(Image4,:)-NEW_StackPositions_pixels(Image3,:));
        P2_Proposed = round(NEW_StackPositions_pixels(Image2,:)-NEW_StackPositions_pixels(Image1,:));
        P5_Proposed = round(NEW_StackPositions_pixels(Image3,:)-NEW_StackPositions_pixels(Image1,:));
        P4_Proposed = round(NEW_StackPositions_pixels(Image4,:)-NEW_StackPositions_pixels(Image2,:));
    
    %     PlotAM(Trace_AM_G{1},Trace_r_G{1}+[NEW_StackPositions_pixels(Image1,2),NEW_StackPositions_pixels(Image1,1),NEW_StackPositions_pixels(Image1,3)]-1,'r');
    %     PlotAM(Trace_AM_A{1},Trace_r_A{1}+[NEW_StackPositions_pixels(Image2,2),NEW_StackPositions_pixels(Image2,1),NEW_StackPositions_pixels(Image2,3)]-1,'g');
    %     hold on;
    
%     P2_Proposed = [864,-9,-1];P3_Proposed = [862,-10,-1]; P4_Proposed = [5,1453,1]; P5_Proposed = [3,1453,1];
    stitchedImage_Proposed = StitichAll(Image1,Image2,Image3,Image4,Dataset,P2_Proposed,P3_Proposed,P4_Proposed);
    figure(),imshow(max(stitchedImage_Proposed,[],3),[0 max(stitchedImage_Proposed(:))]);title('Proposed');
    hold on;
    PlotAM(Trace_AM_G{1},Trace_r_G{1}+1,'r');
    PlotAM(Trace_AM_A{1},Trace_r_A{1}+P2_Proposed([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{2},Trace_r_G{2}+1,'r');
    PlotAM(Trace_AM_A{2},Trace_r_A{2}+P4_Proposed([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{3},Trace_r_G{3}+abs(P2_Proposed([2,1,3]))+1,'r');
    PlotAM(Trace_AM_A{3},Trace_r_A{3}+abs(P2_Proposed([2,1,3]))+P5_Proposed([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{4},Trace_r_G{4}+P4_Proposed([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{4},Trace_r_A{4}+P4_Proposed([2,1,3])+P3_Proposed([2,1,3])+1,'g');
    
    AM1 = Trace_AM_G{1};
    AM2 = Trace_AM_G{2};
    AM3 = Trace_AM_G{3};
    AM4 = Trace_AM_G{4};
    AMG=AM1;
    AMG(end+1:end+size(AM2,1),end+1:end+size(AM2,2)) = AM2;
    AMG(end+1:end+size(AM3,1),end+1:end+size(AM3,2)) = AM3;
    AMG(end+1:end+size(AM4,1),end+1:end+size(AM4,2)) = AM4;
    AM1 = Trace_AM_A{1};
    AM2 = Trace_AM_A{2};
    AM3 = Trace_AM_A{3};
    AM4 = Trace_AM_A{4};
    AMA=AM1;
    AMA(end+1:end+size(AM2,1),end+1:end+size(AM2,2)) = AM2;
    AMA(end+1:end+size(AM3,1),end+1:end+size(AM3,2)) = AM3;
    AMA(end+1:end+size(AM4,1),end+1:end+size(AM4,2)) = AM4;
    
    r_G = [Trace_r_G{1};Trace_r_G{2};Trace_r_G{3}+P2_Proposed([2,1,3]);Trace_r_G{4}+P4_Proposed([2,1,3])];
    r_A = [Trace_r_A{1}+P2_Proposed([2,1,3])+1;Trace_r_A{2}+P4_Proposed([2,1,3]);Trace_r_A{3}+P2_Proposed([2,1,3])+P5_Proposed([2,1,3]);Trace_r_A{4}+P4_Proposed([2,1,3])+P3_Proposed([2,1,3])];
    
    
    AM = AMG;
    AM(end+1:end+size(AMA,1),end+1:end+size(AMA,2)) = AMA; 
    IM = [];
    Original = uint16(stitchedImage_Proposed);
    r = [r_G;r_A];
    R = zeros(size(r,1),1);
    reduction_x = 1;
    reduction_y = 1;
    reduction_z = 1;
    save('C:\Users\Seyed\Desktop\Test\visual.mat','AM','IM','Original','r','R','reduction_x','reduction_y','reduction_z','-v7.3');
    
    
    P2_Fiji = [867,-9,0];P3_Fiji = [850,-9,1]; P4_Fiji = [3,1382,0];P5_Fiji=[9,1382,0];
    stitchedImage_Fiji = StitichAll(Image1,Image2,Image3,Image4,Dataset,P2_Fiji,P3_Fiji,P4_Fiji);
    figure(),imshow(max(stitchedImage_Fiji,[],3),[0 max(stitchedImage_Fiji(:))]);title('Fiji');
    hold on;
    PlotAM(Trace_AM_G{1},Trace_r_G{1}+1,'r');
    PlotAM(Trace_AM_A{1},Trace_r_A{1}+P2_Fiji([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{2},Trace_r_G{2}+1,'r');
    PlotAM(Trace_AM_A{2},Trace_r_A{2}+P4_Fiji([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{3},Trace_r_G{3}+P2_Fiji([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{3},Trace_r_A{3}+P2_Fiji([2,1,3])+P5_Fiji([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{4},Trace_r_G{4}+P4_Fiji([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{4},Trace_r_A{4}+P4_Fiji([2,1,3])+P3_Fiji([2,1,3])+1,'g');
    
    
    P2_XUV = [862,-9,0];P3_XUV = [864,0,0]; P4_XUV = [3,1453,1];P5_XUV=[2,1452,-1];
    stitchedImage_XUV = StitichAll(Image1,Image2,Image3,Image4,Dataset,P2_XUV,P3_XUV,P4_XUV);
    figure(),imshow(max(stitchedImage_XUV,[],3),[0 max(stitchedImage_XUV(:))]);title('XUV');
    hold on;
    PlotAM(Trace_AM_G{1},Trace_r_G{1}+1,'r');
    PlotAM(Trace_AM_A{1},Trace_r_A{1}+P2_XUV([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{2},Trace_r_G{2}+1,'r');
    PlotAM(Trace_AM_A{2},Trace_r_A{2}+P4_XUV([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{3},Trace_r_G{3}+P2_XUV([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{3},Trace_r_A{3}+P2_XUV([2,1,3])+P5_XUV([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{4},Trace_r_G{4}+P4_XUV([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{4},Trace_r_A{4}+P4_XUV([2,1,3])+P3_XUV([2,1,3])+1,'g');
    
    P2_Tera = [859,-9,1];P3_Tera = [859,-9,1]; P4_Tera = [3,1453,1];P5_Tera=[3,1453,1];
    stitchedImage_Tera = StitichAll(Image1,Image2,Image3,Image4,Dataset,P2_Tera,P3_Tera,P4_Tera);
    figure(),imshow(max(stitchedImage_Tera,[],3),[0 max(stitchedImage_Tera(:))]);title('Tera');
    hold on;
    PlotAM(Trace_AM_G{1},Trace_r_G{1}+1,'r');
    PlotAM(Trace_AM_A{1},Trace_r_A{1}+P2_Tera([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{2},Trace_r_G{2}+1,'r');
    PlotAM(Trace_AM_A{2},Trace_r_A{2}+P4_Tera([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{3},Trace_r_G{3}+P2_Tera([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{3},Trace_r_A{3}+P2_Tera([2,1,3])+P5_Tera([2,1,3])+1,'g');
    hold on;
    PlotAM(Trace_AM_G{4},Trace_r_G{4}+P4_Tera([2,1,3])+1,'r');
    PlotAM(Trace_AM_A{4},Trace_r_A{4}+P4_Tera([2,1,3])+P3_Tera([2,1,3])+1,'g');
    
end
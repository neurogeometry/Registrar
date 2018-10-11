
close all;
clear all;
% clc;
addpath ('NeuronTracerV20');
addpath ('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\Functions');
CalculateOtherMEthods = 1;
usePoints = 0;
% DatasetList = {'Holtmaat','Neuromuscular','Neocortical_1','MouseLight','Visual'};
showCorr = 0;
DatasetList = {'MouseLight'};
doCorrelation = 1;
CSVData = 'C:\Users\Seyed\Documents\NeurogeometryLab\NeurogeometryLab\Seyed\Evaluation\data\evaluation\';
db = 1;
for Dataset = DatasetList
    fprintf('%s\n',Dataset{1});
    
    if strcmp(Dataset,'Holtmaat')
        pixelSize = [0.26 0.26 0.8];
        load('data/StackData_Holtmaat_ch0.mat');
        Evaluation_csv_pth = [CSVData,'Holtmaat\Holtmaat_Evaluation.csv'];
        pad = [0 0 0];
        
        ResultFolder = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\Results-Holtmaat_StackList - Full\';
        
        StackPositions_pixels = StackPositions_pixels(:,[2,1,3]);
        
    elseif strcmp(Dataset,'Neuromuscular')
        pixelSize = [0.1 0.1 0.21];
        load('data/StackData_Neuromuscular1.mat');
        Evaluation_csv_pth = [CSVData,'DIADEM_Neuro\Neuro_Evaluation.csv'];
        
        StackPositions_pixels_4Stacks(1,:) = StackPositions_pixels(106,:);
        StackPositions_pixels_4Stacks(2,:) = StackPositions_pixels(107,:);
        StackPositions_pixels_4Stacks(3,:) = StackPositions_pixels(114,:);
        StackPositions_pixels_4Stacks(4,:) = StackPositions_pixels(115,:);
        StackList_4Stacks(1,:) = StackList(106,:);
        StackList_4Stacks(2,:) = StackList(107,:);
        StackList_4Stacks(3,:) = StackList(114,:);
        StackList_4Stacks(4,:) = StackList(115,:);
        StackList = StackList_4Stacks;
        pad = [0 0 0];
        
        ResultFolder = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\Results-Neuromuscular_StackList_4Stacks\';
        
        StackPositions_pixels = StackPositions_pixels_4Stacks(:,[2,1,3]);
    elseif strcmp(Dataset,'Neocortical_1')
        pixelSize = [0.3 0.3 1.02];
        load('data/StackData_DIADEM.mat');
        Evaluation_csv_pth = [CSVData,'DIADEM_Neo1\Neo1_Evaluation.csv'];
        pad = [0 0 0];
        
        ResultFolder = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\Results-Neocortical1_StackList\';
        StackPositions_pixels = StackPositions_pixels(:,[2,1,3]);
        
    elseif strcmp(Dataset,'Neocortical_2')
        pixelSize = [0.3 0.3 1.02];
        load('data/StackData_DIADEM2.mat');
        Evaluation_csv_pth = [CSVData,'DIADEM_Neo2\Neo2_Evaluation.csv'];
        
        ResultFolder = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\Results-Neocortical2_StackList\';
        
        pad = [0 0 0];
        
        StackPositions_pixels = StackPositions_pixels(:,[2,1,3]);
    elseif strcmp(Dataset,'MouseLight')
        pixelSize = [0.377607421875 0.277486979166667 0.99601593625498];
        load('data/StackData.mat');
        Dataset = 'MouseLight';
        pad = [0 30 0];
        
        Evaluation_csv_pth = [CSVData,'MouseLight\MouseLight_Evaluation.csv'];
        ResultFolder = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\Results-MouseLight_StackList\';
        
        
        temp=StackPositions_pixels;
        StackPositions_pixels(:,1) = max(temp(:,2))-temp(:,2);
        StackPositions_pixels(:,2) = max(temp(:,1))-temp(:,1);
        StackPositions_pixels(:,3) = max(temp(:,3))-temp(:,3); % Not Sure - Need to check
    elseif strcmp(Dataset,'Visual')
        pixelSize = [0.3 0.3 0.879];
        load('data/StackData_Visual.mat');
        Evaluation_csv_pth = [CSVData,'Visual\Visual_Evaluation.csv'];
        pad = [0 0 0];
        
        ResultFolder = 'C:\Users\Seyed\Documents\DatasetTests\MicroscopeFiles\Results-Visual_StackList\';
        
        StackPositions_pixels = StackPositions_pixels(:,[2,1,3]);
    end
    
    
    FeaturePositions_T = load([ResultFolder,'MatchedPoints_Translation.mat']);
    FeaturePositions_R = load([ResultFolder,'MatchedPoints_Rigid.mat']);
    FeaturePositions_A = load([ResultFolder,'MatchedPoints_Affine.mat']);
    
    Transformation_T =  load([ResultFolder,'T_Translation.mat']);
    Transformation_R = load([ResultFolder,'T_Rigid.mat']);
    Transformation_A = load([ResultFolder,'T_Affine.mat']);
    EvaluationInfo = table2cell(readtable(Evaluation_csv_pth,'Delimiter',','));
    
    for k = 1: size(EvaluationInfo,1)
        sourceID = EvaluationInfo{k,1};
        targetID = EvaluationInfo{k,2};
        
        if strcmp(Dataset,'MouseLight')
            Source_Stack_File = [char(StackList(sourceID,2)),'\',char(StackList(sourceID,1)),'-ngc.0.tif'];
            Target_Stack_File = [char(StackList(targetID,2)),'\',char(StackList(targetID,1)),'-ngc.0.tif'];
        else
            Source_Stack_File = [char(StackList(sourceID,2)),'\',char(StackList(sourceID,1)),'.tif'];
            Target_Stack_File = [char(StackList(targetID,2)),'\',char(StackList(targetID,1)),'.tif'];
            
        end
        IM_Source=ImportStack(char(Source_Stack_File));
        IM_Target=ImportStack(char(Target_Stack_File));
        
        
        temp=StackList{sourceID,1};
        sourcePath = ['data/evaluation/',char(EvaluationInfo{k,3}),'/',temp,'_opt.swc'];
        temp=StackList{targetID,1};
        targetPath = ['data/evaluation/',char(EvaluationInfo{k,3}),'/',temp,'_opt.swc'];
        
        if CalculateOtherMEthods
            dx_Fiji = EvaluationInfo{k,4};
            dy_Fiji = EvaluationInfo{k,5};
            dz_Fiji = EvaluationInfo{k,6};
            
            dx_XUV = EvaluationInfo{k,7};
            dy_XUV = EvaluationInfo{k,8};
            dz_XUV = EvaluationInfo{k,9};
            
            dx_Tara = EvaluationInfo{k,10};
            dy_Tara = EvaluationInfo{k,11};
            dz_Tara = EvaluationInfo{k,12};
        end
        ppm=EvaluationInfo{k,13};
        
        
        Source_StackPositions = round(StackPositions_pixels(sourceID,:));
        Target_StackPositions = round(StackPositions_pixels(targetID,:));
        
        [AM_Source,r_Source,R_Source]=swc2AM(sourcePath);
        [AM_Target,r_Target,R_Target]=swc2AM(targetPath);
        
        [AM_Source,r_Source,~] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
        [AM_Target,r_Target,~] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);
        
        
        % % Translation
        if usePoints
            Global_Matched_Source = FeaturePositions_T.Matched{sourceID,targetID}(:,1:3)'+Source_StackPositions'-1;
            Global_Matched_Target = FeaturePositions_T.Matched{sourceID,targetID}(:,4:6)'+Target_StackPositions'-1;
            b=Optimal_Translation_Transform(Global_Matched_Source,Global_Matched_Target);
            r_Source_TranslationPairs = r_Source+Source_StackPositions+b'-1;
            r_Target_TranslationPairs = r_Target+Target_StackPositions-1;
        else
            b = Transformation_T.T.b;
            r_Source_TranslationPairs = r_Source+Source_StackPositions+b(:,sourceID)'-1;
            r_Target_TranslationPairs = r_Target+Target_StackPositions+b(:,targetID)'-1;
        end
        %                     figure
        %                         PlotAM(AM_Source,r_Source_TranslationPairs,'r')
        %                         axis equal
        %                         PlotAM(AM_Target,r_Target_TranslationPairs,'g')
        %                         axis equal
        [Distances_Translation_um,Distances_Translation_voxel] = TraceDistance(AM_Source, r_Source_TranslationPairs, AM_Target, r_Target_TranslationPairs,pixelSize,0);
        Translation_voxel = mean(Distances_Translation_voxel);
        Translation_um = mean(Distances_Translation_um);
        if doCorrelation
            
           corr_Translation = caclulateCorr(Source_StackPositions+b(:,sourceID)',Target_StackPositions,IM_Source,IM_Target,pad,showCorr);
  %          CC=Stack_Correlation(IM_Source,IM_Target,Source_StackPositions+b(:,sourceID)',Target_StackPositions)
        end
        
        % Rigid
        if usePoints
            Global_Matched_Source = FeaturePositions_R.Matched{sourceID,targetID}(:,1:3)'+Source_StackPositions'-1;
            Global_Matched_Target = FeaturePositions_R.Matched{sourceID,targetID}(:,4:6)'+Target_StackPositions'-1;
            [R,b]=Optimal_Rigid_Transform(Global_Matched_Source,Global_Matched_Target);
            r_Source_Rigid = (R*(r_Source+Source_StackPositions-1)'+b)';
            r_Target_Rigid = r_Target+Target_StackPositions-1;
        else
            b = Transformation_R.T.b;
            R_source = Transformation_R.T.L(:,(3*sourceID)-2:3*sourceID);
            R_target = Transformation_R.T.L(:,(3*targetID)-2:3*targetID);
            r_Source_Rigid = (R_source*(r_Source+Source_StackPositions-1)'+b(:,sourceID))';
            r_Target_Rigid = (R_target*(r_Target+Target_StackPositions-1)'+b(:,targetID))';
        end
        %                 figure
        %                 PlotAM(AM_Source,r_Source_Rigid,'r')
        %                 axis equal
        %                 PlotAM(AM_Target,r_Target_Rigid,'g')
        %                 axis equal
        [Distances_Rigid_um,Distances_Rigid_voxel] = TraceDistance(AM_Source, r_Source_Rigid, AM_Target, r_Target_Rigid,pixelSize,0);
        
        Rigid_voxel = mean(Distances_Rigid_voxel);
        Rigid_um = mean(Distances_Rigid_um);
        
        if doCorrelation
            [IM_Source_prime,StackPosition_prime]=Perform_Linear_Transform(IM_Source,(R_source*Source_StackPositions'+b(:,sourceID))',R_source,b(:,sourceID));
%             figure;imshow(max(IM_Source_prime,[],3),[0 max(IM_Source_prime(:))]);
            [IM_Target_prime,StackPosition_prime]=Perform_Linear_Transform(IM_Target,(R_target*Target_StackPositions'+b(:,targetID))',R_target,b(:,targetID));
%             figure;imshow(max(IM_Target_prime,[],3),[0 max(IM_Target_prime(:))]);
            corr_Rigid = caclulateCorr((R_source*Source_StackPositions'+b(:,sourceID))',(R_target*Target_StackPositions'+b(:,targetID))',IM_Source_prime,IM_Target_prime,pad,showCorr);
%             CC=Stack_Correlation(IM_Source_prime,IM_Target_prime,(R_source*Source_StackPositions'+b(:,sourceID))',(R_target*Target_StackPositions'+b(:,targetID))')
        end
        
        % Affine
        if usePoints
            Global_Matched_Source = FeaturePositions_A.Matched{sourceID,targetID}(:,1:3)'+Source_StackPositions'-1;
            Global_Matched_Target = FeaturePositions_A.Matched{sourceID,targetID}(:,4:6)'+Target_StackPositions'-1;
            [L,b]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target);
            r_Source_Affine = (L*(r_Source+Source_StackPositions-1)'+b)';
            r_Target_Affine = r_Target+Target_StackPositions-1;
        else
            b = Transformation_A.T.b;
            L_source = Transformation_A.T.L(:,(3*sourceID)-2:3*sourceID);
            L_target = Transformation_A.T.L(:,(3*targetID)-2:3*targetID);
            r_Source_Affine = (L_source*(r_Source+Source_StackPositions-1)'+b(:,sourceID))';
            r_Target_Affine = (L_target*(r_Target+Target_StackPositions-1)'+b(:,targetID))';
        end
        
        %         figure
        %         PlotAM(AM_Source,r_Source_Affine,'r')
        %         axis equal
        %         PlotAM(AM_Target,r_Target_Affine,'g')
        %         axis equal
        [Distances_Affine_um,Distances_Affine_voxel] = TraceDistance(AM_Source, r_Source_Affine, AM_Target, r_Target_Affine,pixelSize,0);
        
        Affine_voxel = mean(Distances_Affine_voxel);
        Affine_um = mean(Distances_Affine_um);
        if doCorrelation
            
            [IM_Source_prime,StackPosition_prime]=Perform_Linear_Transform(IM_Source,(L_source*Source_StackPositions'+b(:,sourceID))',L_source,b(:,sourceID));
%             figure;imshow(max(IM_Source_prime,[],3),[0 max(IM_Source_prime(:))]);
            [IM_Target_prime,StackPosition_prime]=Perform_Linear_Transform(IM_Target,(L_target*Target_StackPositions'+b(:,targetID))',L_target,b(:,targetID));
%             figure;imshow(max(IM_Target_prime,[],3),[0 max(IM_Target_prime(:))]);
            corr_Affine = caclulateCorr((L_source*Source_StackPositions'+b(:,sourceID))',(L_target*Target_StackPositions'+b(:,targetID))',IM_Source_prime,IM_Target_prime,pad,showCorr);
%         CC=Stack_Correlation(IM_Source_prime,IM_Target_prime,(L_source*Source_StackPositions'+b(:,sourceID))',(L_target*Target_StackPositions'+b(:,targetID))');
%             corr_Affine = caclulateCorr((L*Source_StackPositions'+b)',Target_StackPositions,IM_Source,IM_Target,pad,showCorr);
        end
        % AffineNonRigid
        Global_Matched_Source = FeaturePositions_A.Matched{sourceID,targetID}(:,1:3)'+Source_StackPositions'-1;
        Global_Matched_Target = FeaturePositions_A.Matched{sourceID,targetID}(:,4:6)'+Target_StackPositions'-1;
        [~,L,b]=Optimal_Affine_Transform(Global_Matched_Source,Global_Matched_Target,0.1);
        
        r_Source_Affine = (L*(r_Source+Source_StackPositions-1)'+b)';
        r_Target_Affine = r_Target+Target_StackPositions-1;
        
        N_L=1;
        temp = (L*Global_Matched_Source+b);
        Min=min([min(temp,[],2),min(Global_Matched_Target,[],2)],[],2);
        Max=max([max(temp,[],2),max(Global_Matched_Target,[],2)],[],2);
        [~,XYZlmn,N_L,Min,Max]=Optimal_Nonrigid_Transform(temp,Global_Matched_Target,N_L,Min,Max);
        r_Source_Affine=Perform_Nonrigid_Transform(r_Source_Affine',XYZlmn,N_L,Min,Max)';
        
        %         figure
        %         PlotAM(AM_Source,r_Source_Affine,'r')
        %         axis equal
        %         PlotAM(AM_Target,r_Target_Affine,'g')
        %         axis equal
        [Distances_AffineNonRigid_um,Distances_AffineNonRigid_voxel] = TraceDistance(AM_Source, r_Source_Affine, AM_Target, r_Target_Affine,pixelSize,0);
        
        AffineNonRigid_voxel = mean(Distances_AffineNonRigid_voxel);
        AffineNonRigid_um = mean(Distances_AffineNonRigid_um);
        
        
        if doCorrelation
            Pos_Source_Affine = (L*(Source_StackPositions-1)'+b)';
            Source_StackPositions_nonrigid = Perform_Nonrigid_Transform(Pos_Source_Affine',XYZlmn,N_L,Min,Max);
            corr_NonRigid = caclulateCorr(Source_StackPositions_nonrigid',Target_StackPositions,IM_Source,IM_Target,pad,showCorr);
        end
        
        showTraces = 0;
        r_Source_old = r_Source + [StackPositions_pixels(sourceID,1),StackPositions_pixels(sourceID,2),StackPositions_pixels(sourceID,3)]-1;
        r_Target_old = r_Target + [StackPositions_pixels(targetID,1),StackPositions_pixels(targetID,2),StackPositions_pixels(targetID,3)]-1;
        
        if CalculateOtherMEthods
            r_Source_Fiji = r_Source;
            r_Target_Fiji = r_Target +[dy_Fiji,dx_Fiji,dz_Fiji];
            
            r_Source_Tara = r_Source;%
            r_Target_Tara = r_Target +[dy_Tara,dx_Tara,dz_Tara];
            
            r_Source_XUV = r_Source;
            r_Target_XUV = r_Target +[dy_XUV,dx_XUV,dz_XUV];
            
            if doCorrelation
                corr_Fiji = caclulateCorr([0,0,0],[dy_Fiji,dx_Fiji,dz_Fiji],IM_Source,IM_Target,pad,0);
                corr_Tera = caclulateCorr([0,0,0],[dy_Tara,dx_Tara,dz_Tara],IM_Source,IM_Target,pad,showCorr);
                corr_XUV = caclulateCorr([0,0,0],[dy_XUV,dx_XUV,dz_XUV],IM_Source,IM_Target,pad,showCorr);
            end
            [Distances_Fiji_um,Distances_Fiji_voxel] = TraceDistance(AM_Source, r_Source_Fiji, AM_Target, r_Target_Fiji,pixelSize,showTraces);
            [Distances_Tera_um,Distances_Tera_voxel] = TraceDistance(AM_Source, r_Source_Tara, AM_Target, r_Target_Tara,pixelSize,showTraces);
            [Distances_XUV_um,Distances_XUV_voxel] = TraceDistance(AM_Source, r_Source_XUV, AM_Target, r_Target_XUV,pixelSize,showTraces);
        end
        
        [Distances_Old_um,Distances_Old_voxel] = TraceDistance(AM_Source, r_Source_old, AM_Target, r_Target_old,pixelSize,showTraces);
        corr_old = caclulateCorr(StackPositions_pixels(sourceID,:),StackPositions_pixels(targetID,:),IM_Source,IM_Target,pad,showCorr);
        
        M_Translation(db,k,1) = mean(Distances_Translation_voxel);
        M_Rigid(db,k,1) = mean(Distances_Rigid_voxel);
        M_Affine(db,k,1) = mean(Distances_Affine_voxel);
        M_NonRigid(db,k,1) = mean(Distances_AffineNonRigid_voxel);
        if CalculateOtherMEthods
            M_Fiji(db,k,1) = mean(Distances_Fiji_voxel);
            M_Tera(db,k,1) = mean(Distances_Tera_voxel);
            M_XUV(db,k,1) = mean(Distances_XUV_voxel);
        else
            M_Fiji(db,k,1) = 0;
            M_Tera(db,k,1) = 0;
            M_XUV(db,k,1) = 0;
        end
        M_Old(db,k,1) = mean(Distances_Old_voxel);
        
        M_Translation(db,k,2) = mean(Distances_Translation_um);
        M_Affine(db,k,2) = mean(Distances_Affine_um);
        M_Rigid(db,k,2) = mean(Distances_Rigid_um);
        M_NonRigid(db,k,2) = mean(Distances_AffineNonRigid_um);
        if CalculateOtherMEthods
            M_Fiji(db,k,2) = mean(Distances_Fiji_um);
            M_Tera(db,k,2) = mean(Distances_Tera_um);
            M_XUV(db,k,2) = mean(Distances_XUV_um);
        else
            M_Fiji(db,k,2) =0;
            M_Tera(db,k,2) = 0;
            M_XUV(db,k,2) = 0;
        end
        M_Old(db,k,2) = mean(Distances_Old_um);
        
        
        
        if doCorrelation
            M_Translation(db,k,3) = corr_Translation;
            M_Rigid(db,k,3) = corr_Rigid;
            M_Affine(db,k,3) = corr_Affine;
            M_NonRigid(db,k,3) = corr_NonRigid;
            if CalculateOtherMEthods
                M_Fiji(db,k,3) = corr_Fiji;
                M_Tera(db,k,3) = corr_Tera;
                M_XUV(db,k,3) = corr_XUV;
            else
                M_Fiji(db,k,3) = 0;
                M_Tera(db,k,3) = 0;
                M_XUV(db,k,3) = 0;
            end
            M_Old(db,k,3) = corr_old;
        end
        k
        [M_Old(db,k,1),M_Translation(db,k,1),M_Rigid(db,k,1),M_Affine(db,k,1),M_NonRigid(db,k,1),M_Fiji(db,k,1),M_XUV(db,k,1),M_Tera(db,k,1)]
        [M_Old(db,k,2),M_Translation(db,k,2),M_Rigid(db,k,2),M_Affine(db,k,2),M_NonRigid(db,k,2),M_Fiji(db,k,2),M_XUV(db,k,2),M_Tera(db,k,2)]
        
        
    end

    All_voxels(db,:) = [mean(M_Old(db,:,1)),mean(M_Translation(db,:,1)),mean(M_Rigid(db,:,1)),mean(M_Affine(db,:,1)),mean(M_NonRigid(db,:,1)),mean(M_Fiji(db,:,1)),mean(M_XUV(db,:,1)),mean(M_Tera(db,:,1))];
    All_um(db,:) = [mean(M_Old(db,:,2)),mean(M_Translation(db,:,2)),mean(M_Rigid(db,:,2)),mean(M_Affine(db,:,2)),mean(M_NonRigid(db,:,2)),mean(M_Fiji(db,:,2)),mean(M_XUV(db,:,2)),mean(M_Tera(db,:,2))]
    if doCorrelation
        All_corr(db,:) = [mean(M_Old(db,:,3)),mean(M_Translation(db,1:1,3)),mean(M_Rigid(db,:,3)),mean(M_Affine(db,:,3)),mean(M_NonRigid(db,:,3)),mean(M_Fiji(db,:,3)),mean(M_XUV(db,:,3)),mean(M_Tera(db,:,3))];
    end

end


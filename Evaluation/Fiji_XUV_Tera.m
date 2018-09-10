close all;
clear all;
clc;
addpath ('NeuronTracerV20');
addpath ('C:\Users\Seyed\Documents\DatasetTests\GUI\NCT_Registration\Functions');

usePoints = 0;
DatasetList = {'Holtmaat','Neuromuscular','Neocortical_1','MouseLight','Visual'};
showCorr = 0;
showTraces = 0;
% DatasetList = {'MouseLight'};
doCorrelation = 0;
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
        
        
        dx_Fiji = EvaluationInfo{k,4};
        dy_Fiji = EvaluationInfo{k,5};
        dz_Fiji = EvaluationInfo{k,6};
        
        dx_XUV = EvaluationInfo{k,7};
        dy_XUV = EvaluationInfo{k,8};
        dz_XUV = EvaluationInfo{k,9};
        
        dx_Tara = EvaluationInfo{k,10};
        dy_Tara = EvaluationInfo{k,11};
        dz_Tara = EvaluationInfo{k,12};

        ppm=EvaluationInfo{k,13};
        
        Source_StackPositions = round(StackPositions_pixels(sourceID,:));
        Target_StackPositions = round(StackPositions_pixels(targetID,:));
        
        [AM_Source,r_Source,R_Source]=swc2AM(sourcePath);
        [AM_Target,r_Target,R_Target]=swc2AM(targetPath);
        
        [AM_Source,r_Source,~] = AdjustPPM(AM_Source,r_Source,R_Source,ppm);
        [AM_Target,r_Target,~] = AdjustPPM(AM_Target,r_Target,R_Target,ppm);
        
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
        
        M_Fiji(db,k,1) = mean(Distances_Fiji_voxel);
        M_Tera(db,k,1) = mean(Distances_Tera_voxel);
        M_XUV(db,k,1) = mean(Distances_XUV_voxel);
        
        M_Fiji(db,k,2) = mean(Distances_Fiji_um);
        M_Tera(db,k,2) = mean(Distances_Tera_um);
        M_XUV(db,k,2) = mean(Distances_XUV_um);
        
        if doCorrelation
            M_Fiji(db,k,3) = corr_Fiji;
            M_Tera(db,k,3) = corr_Tera;
            M_XUV(db,k,3) = corr_XUV;
        end
        
        k
    end
    
    All_voxels(db,:) = [mean(M_Fiji(db,:,1)),mean(M_XUV(db,:,1)),mean(M_Tera(db,:,1))];
    All_um(db,:) = [mean(M_Fiji(db,:,2)),mean(M_XUV(db,:,2)),mean(M_Tera(db,:,2))]
    if doCorrelation
        All_corr(db,:) = [mean(M_Fiji(db,:,3)),mean(M_XUV(db,:,3)),mean(M_Tera(db,:,3))];
    end
    db = db +1;
    
end
function registeration (varargin)
% ============================== About ====================================
% -------------------------------------------------------------------------
%
% Purpose: Main Registration function to call other functions
% Input:
%   k                  : 1*1, The folder address of stack
%   Mesh                    : 1*1, The file name of stack
%   StackList_csv_pth          : 1*1, The address to save the features
%   filterValue          : 1*1, The address to save the features
%   TransformationValue          : 1*1, The address to save the features
%   Seq_Par          : 1*1, The address to save the features
%   Par_workers          : 1*1, The address to save the features
%   blendingSID          : 1*1, The address to save the features
%
% -------------------------------------------------------------------------
% Author: Seyed Mostafa Mousavi Kahaki, Armen Stepanyants
% Northeastern University, USA
% kahaki@neu.edu, a.stepanyants.neu.edu
% =========================================================================
% -------------------------------------------------------------------------

% Load Parameters
parameters;

% Define Initial Variables
StackList_csv_pth = varargin{1};
TransformationValue = varargin{2};
Seq_Par = varargin{3};
Par_workers = varargin{4};
blendingSID = varargin{5};
handles = varargin{6};
runFeatureExtraction = handles.checkbox3.Value;
runFeatureMatching = handles.checkbox4.Value;
runGlobalOpt = handles.checkbox5.Value;
runBlending = handles.checkbox6.Value;
TypeOfRegistration = handles.v.Value;
runRetilling = handles.chkretilling.Value;
outputType = handles.popupretilling.Value;
debug = handles.chkdebug.Value;
AllGUI = NCT_Registration;
LogLine = 1;
stop = 0;


% Check if CSV file exist
if exist(StackList_csv_pth,'file') > 0
    StackList = table2cell(readtable(StackList_csv_pth,'Delimiter',','));
    
    [PathStr,FolderName]=fileparts(StackList_csv_pth);
    DataFolder=[PathStr,'/Results-',FolderName];
    
    handles.axes1.YLabel.String = '';
    handles.axes1.XLabel.String='';
    handles.axes1.Title.String = '';
    handles.axes5.YLabel.String = '';
    handles.axes5.XLabel.String='';
    handles.axes5.Title.String = '';
    handles.axes6.YLabel.String = '';
    handles.axes6.XLabel.String='';
    handles.axes6.Title.String = '';
    handles.slider2.set('Visible','off');
    cla(handles.axes5);
    cla(handles.axes6);
    axes_main = findall(AllGUI,'Tag','axes1');
    axes_main = axes_main(1);

    % index to remove invalid files
    errIndxs = [];
    StackSizes_pixels = zeros(size(StackList,1),3);
    for i = 1:size(StackList,1)
        try
            [filepath,~,ext] = fileparts(char(StackList(i,1)));
            if size(ext,2)>1
                InfoImage=imfinfo(char(StackList(i,1)));
                StackSizes_pixels(i,1) = InfoImage.Height;
                StackSizes_pixels(i,2) = InfoImage.Width;
                StackSizes_pixels(i,3) = size(InfoImage,1);
            else
                allfiles = dir(filepath);
                NumFiles=length(allfiles);
                imgIdx = 0;
                for j = 1:NumFiles
                    [~,~,ext] = fileparts(allfiles(j).name);
                    if strcmpi(ext,'.tif') || strcmpi(ext,'.jp2') || strcmpi(ext,'.png') || strcmpi(ext,'.jpeg')
                        InfoImage=imfinfo(char([allfiles(j).folder,'\',allfiles(j).name]));
                        StackSizes_pixels(i,1) = InfoImage.Height;
                        StackSizes_pixels(i,2) = InfoImage.Width;
                        imgIdx = imgIdx +1;
                    end
                end
                StackSizes_pixels(i,3) = imgIdx;
            end
        catch
            disp(['Can not read the file: ',char(StackList(i,1))]);
            errIndxs = [errIndxs, i];
        end
    end
    StackList (errIndxs,:) = [];
    StackSizes_pixels (errIndxs,:) = [];
    if size(StackList,2) == 7
        StackSizes_pixels = cell2mat(StackList(:,5:7));
    end
    
    % Get positions from the list
    StackPositions_pixels = cell2mat(StackList(:,2:4));%xlsread(StackPositions_pixels_csv_pth);
    
    % This part is not fixed yet, should be based on correct positions from csv
    if TypeOfRegistration == 1
        temp=StackPositions_pixels;
        StackPositions_pixels(:,1) = max(temp(:,2))-temp(:,2);
        StackPositions_pixels(:,2) = max(temp(:,1))-temp(:,1);
        StackPositions_pixels(:,3) = max(temp(:,3))-temp(:,3);
    else
        StackPositions_pixels = StackPositions_pixels(:,[2,1,3]);
    end
    
    % Find stacks overlap
    % overlap = round((params.FE.overlap/100) * StackSizes_pixels(1,:));
    handles.axes2.Visible = 'on';
    All_overlaps = FindOverlaps(StackPositions_pixels,StackSizes_pixels,StackList,debug);
    
    if TypeOfRegistration == 8 || TypeOfRegistration == 7
        All_overlaps = zeros(size(All_overlaps));
        All_overlaps(2:size(All_overlaps,1)+1:end) = 1;
        All_overlaps = All_overlaps';
    end
    
    % variable 'LogLine' is handling the log
    listbox_log{LogLine}  = 'Registration Process Started';
    handles.listbox1.String = listbox_log;drawnow
    LogLine = LogLine + 1;
    handles.listbox1.Value = LogLine-1;
    handles.axes1.Children.delete
    
    if runFeatureExtraction
        if(exist(DataFolder, 'dir')>0 )
            if(exist([DataFolder,'/features/'], 'dir')>0 )
                if ispc
                    dos_cmd = sprintf( 'rmdir /S /Q "%s"', [DataFolder,'/features/'] );
                    [ ~, ~ ] = system( dos_cmd );
                else
                    dos_cmd = sprintf( 'rmdir /S /Q "%s"', [DataFolder,'/features/'] );
                    [ ~, ~ ] = system( dos_cmd );
                end
            end
        else
            mkdir(DataFolder);
        end
        
        [~,listbox_log,LogLine,stop]=FeatureExtraction(LogLine,StackList,listbox_log,DataFolder,Seq_Par,Par_workers,All_overlaps,StackPositions_pixels,StackSizes_pixels,debug);
        
    end
    
    if runFeatureMatching && ~stop
        if exist([DataFolder,'/tmp'],'dir') > 0
            
            listbox_log{LogLine}  = 'Feature Matchig Started';
            handles.listbox1.String = listbox_log;drawnow
            LogLine = LogLine + 1;
            handles.listbox1.Value = LogLine-1;drawnow
            
            [~,Matched,listbox_log,LogLine,stop]=Generate_Reg_MatchedPoints(listbox_log,LogLine,All_overlaps,StackList,StackPositions_pixels,StackSizes_pixels,TransformationValue,DataFolder,Seq_Par,Par_workers,debug,params.FM.mu);
        else
            warndlg('The Features Folder is not Exists! Please run Feature Extraction.','!! Warning !!');
        end
        
    end
    runGlobal = 0;
    if runGlobalOpt && ~stop
        if exist([DataFolder,'/tmp'],'dir') > 0
            listbox_log{LogLine}  = 'Global Optimization Started';
            handles.listbox1.String = listbox_log;drawnow
            LogLine = LogLine + 1;
            handles.listbox1.Value = LogLine-1;drawnow
            switch TransformationValue
                case 1
                    if exist([DataFolder,params.FM.MatchedLocationsFile_Translation],'file')
                        load([DataFolder,params.FM.MatchedLocationsFile_Translation],'Matched');
                        runGlobal = 1;
                    else
                        warndlg('The matching file using Translation transform is not exists, please run feature matching using Translation.','!! Warning !!');
                    end
                case 2
                    
                    if exist([DataFolder,params.FM.MatchedLocationsFile_Rigid],'file')
                        load([DataFolder,params.FM.MatchedLocationsFile_Rigid],'Matched');
                        runGlobal = 1;
                    else
                        warndlg('The matching file using Rigid transform is not exists, please run feature matching using Rigid.','!! Warning !!');
                    end
                case 3
                    if exist([DataFolder,params.FM.MatchedLocationsFile_Affine],'file')
                        load([DataFolder,params.FM.MatchedLocationsFile_Affine],'Matched');
                        runGlobal = 1;
                    else
                        warndlg('The matching file using Affine transform is not exists, please run feature matching using Affine.','!! Warning !!');
                    end
                case 4
                    if exist([DataFolder,params.FM.MatchedLocationsFile_Non_Rigid],'file')
                        load([DataFolder,params.FM.MatchedLocationsFile_Non_Rigid],'Matched');
                        runGlobal = 1;
                    else
                        warndlg('The matching file using Non-Rigid transform is not exists, please run feature matching using Non-Rigid.','!! Warning !!');
                    end
                otherwise
                    if exist([DataFolder,params.FM.MatchedLocationsFile_Translation],'file')
                        load([DataFolder,params.FM.MatchedLocationsFile_Translation],'Matched');
                        runGlobal = 1;
                    else
                        warndlg('The matching file using Translation transform is not exists, please run feature matching using Translation.','!! Warning !!');
                    end
            end
            
            if runGlobal && ~stop
                if TypeOfRegistration == 8
                    StackPositions_pixels(:,3) = 0;
                end
                switch TransformationValue
                    case 1
                        Transform_Type = 'Translation';
                        %                       [SaveLocation,T]=Global_Optimal_Translation(StackPositions_pixels,Matched,DataFolder);
                        [SaveLocation,T]=Global_Linear_Transform(StackPositions_pixels,Matched,Transform_Type,DataFolder);
                        %                         T.L = 0;
                    case 2
                        Transform_Type = 'Rigid';
                        %                         [SaveLocation,StackPositions_Registered]=Global_Optimal_Translation(StackPositions_pixels,Matched,DataFolder);
                        [SaveLocation,T]=Global_Linear_Transform(StackPositions_pixels,Matched,Transform_Type,DataFolder);
                        
                    case 3
                        %                [SaveLocation,StackPositions_Registered]=Global_Optimal_Translation(StackPositions_pixels,Matched,DataFolder);
                        %                                                 [SaveLocation,StackPositions_Registered,L,b]=Global_Optimal_Affine(StackPositions_pixels,Matched,DataFolder);
                        %                                                 T.L = L;
                        %                                                 T.b = b;
                        %                                                 T.transform = 'Affine';
                        Transform_Type = 'Affine';
                        [SaveLocation,T]=Global_Linear_Transform(StackPositions_pixels,Matched,Transform_Type,DataFolder);
                    case 4
                        Transform_Type = 'Translation';
                        [SaveLocation,T]=Global_Linear_Transform(StackPositions_pixels,Matched,Transform_Type,DataFolder);
                        T.L = 0;     
                end
                
                listbox_log{LogLine}  = ['Registered Positions Data saved as: ',SaveLocation];
                handles.listbox1.String = listbox_log;drawnow
                LogLine = LogLine + 1;
                handles.listbox1.Value = LogLine-1;drawnow
            else
                warndlg('Global Optimization can not be run.','!! Warning !!');
            end
        else
            warndlg('The Matched Points file is not Exist! Please run Feature Matching.','!! Warning !!');
        end
    end
    
    if runBlending && ~stop
        %         if exist([DataFolder,params.GT.StackPositions_Registered]) > 0 || exist([DataFolder,params.GA.StackPositions_Registered]) > 0 || exist([DataFolder,params.GR.StackPositions_Registered]) > 0
        listbox_log{LogLine}  = 'Blending Started';
        handles.listbox1.String = listbox_log;drawnow
        LogLine = LogLine + 1;
        handles.listbox1.Value = LogLine-1;drawnow

        % 8 for Allignement of Stacks
        if TypeOfRegistration == 8
            [Tile3D_org,Tile3D,stop] = blending_stackreg(StackPositions_Registered,StackSizes_pixels,StackList,T.L,T.b,DataFolder);
        else      
            [Tile3D_org,Tile3D,stop] = blending(StackPositions_pixels,StackSizes_pixels,StackList,blendingSID,T.L,T.b,DataFolder); 
            %                 [Tile3D_org,Tile3D,stop] = blending_mine(StackPositions_pixels,StackSizes_pixels,StackList,blendingSID,T.L,T.b,DataFolder);
        end
        %         Tile3D(all(all(Tile3D == 0,3),2),:,:) = [];
        %         Tile3D(:,all(all(Tile3D == 0,3),1),:) = [];
        %         Tile3D_org(all(all(Tile3D_org == 0,3),2),:,:) = [];
        %         Tile3D_org(:,all(all(Tile3D_org == 0,3),1),:) = [];
        
        handles.radio_before.Enable = 'on';
        handles.radio_after.Enable = 'on';
        handles.z_projection.Enable = 'on';
        handles.radio_layerview.Enable = 'on';
        handles.checkbox10.Enable = 'on';
        handles.axes1.Children.delete
        
        t = findobj(NCT_Registration,'Tag', 'uibuttongroup2');
        u = findobj(NCT_Registration,'Tag', 'radio_after');
        set(t,'SelectedObject',u);
        
        t = findobj(NCT_Registration,'Tag', 'uibuttongroup3');
        u = findobj(NCT_Registration,'Tag', 'z_projection');
        set(t,'SelectedObject',u);
        
        h_im=imshow(max(Tile3D,[],3),[0 max(Tile3D(:))],'Parent',axes_main);
        axes_main=h_im.Parent;
        axes_main.Tag='axes1';
        axes_main.XLim = [0.5 size(h_im.CData,2)+0.5];
        axes_main.YLim = [0.5 size(h_im.CData,1)+0.5];
        listbox_log{LogLine}  = ['Registered Positions Data saved as: ',DataFolder];
        handles.listbox1.String = listbox_log;drawnow
        LogLine = LogLine + 1;
        handles.listbox1.Value = LogLine-1;drawnow
        if handles.checkbox10.Value
            volumeViewer(Tile3D);
        end
        assignin('base','Tile3D',Tile3D);
        assignin('base','Tile3D_org',Tile3D_org);
        %         else
        %             warndlg('The Registered Stacks file is not exist! Please run Global optimization.','!! Warning !!');
        %         end
        
    end
    
    % if Retilling checkbox is checked
    if runRetilling && ~stop
        
        %         TSize = strsplit(handles.editTilesSize.String,',');
        %         TilesSize = [str2double(TSize{1}),str2double(TSize{2}),str2double(TSize{3})];
        switch TransformationValue
            case 1
                load([DataFolder,params.GT.Transformation],'T');
            case 2
                load([DataFolder,params.GR.Transformation],'T');
            case 3
                load([DataFolder,params.GA.Transformation],'T');
            case 4
                load([DataFolder,params.GT.Transformation],'T');
        end
 
        Retiling({StackList{:,1}},StackPositions_pixels,StackSizes_pixels,T,DataFolder,outputType,Seq_Par,Par_workers,DBAddress);
        
    end
    
    if runRetilling && ~stop && (outputType == 1 || outputType == 2)
        ZL = 2;
        NumTiles = inf;
        
        while NumTiles > 1
            NumTiles = CreateZoomLevels(ZL,DataFolder,outputType,DBAddress);
            ZL = ZL * 2;
        end
    end
    
    if ~stop
        listbox_log{LogLine}  = 'Done!';
        handles.listbox1.String = listbox_log;drawnow
        LogLine = LogLine + 1;
        handles.listbox1.Value = LogLine-1;drawnow
    end
else
    warndlg('Please select a valid stack list file','!! Warning !!');
end

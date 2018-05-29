function [FeatureMatchingLog,Matched,listboxItems,v,stop]=Generate_Reg_MatchedPoints(listboxItems,v,All_overlaps,StackList,StackPositions_pixels,StackSizes_pixels,TransformationValue,DataFolder,Seq_Par,Par_workers,debug)
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
addpath('../Functions');
parameters;

tb11 = findobj(NCT_Registration,'Tag', 'pushbutton10');
set(tb11,'userdata',0);
stop = 0;

matchedLocationFile_Translation = [DataFolder,params.FM.MatchedLocationsFile_Translation];
matchedLocationFile_Rigid = [DataFolder,params.FM.MatchedLocationsFile_Rigid];
matchedLocationFile_Affine = [DataFolder,params.FM.MatchedLocationsFile_Affine];
matchedLocationFile_Non_Rigid = [DataFolder,params.FM.MatchedLocationsFile_Non_Rigid];

paramsFMminReqFeatures = 15; % define in parameters file gives error during parallel
paramsFMPad = 0; % define in parameters file gives error during parallel
FeatureMatchingLog = [];
tic
Matched=cell(size(All_overlaps));
Registrationtime = 0;
k = 1;

x = 0:0.1:100;
tb9 = findobj(NCT_Registration,'Tag', 'axes3');
if ~isempty(tb9)
    cla(tb9,'reset')
    tb9.XLim = [0 100];
    patch('XData',[0,0,x(40),x(40)],'YData',[0,20,20,0],'FaceColor','green','Parent',tb9);
    drawnow;
    hold on %hold('units','on')
end

tb_dataset = findobj(NCT_Registration,'Tag', 'v');
% if timelapse or slice registration
if tb_dataset.Value == 8 || tb_dataset.Value == 7
    %     All_overlaps = diag(diag(All_overlaps,+1),+1);
    All_overlaps = zeros(size(All_overlaps));
    All_overlaps(2:size(All_overlaps,1)+1:end) = 1;
    All_overlaps = All_overlaps';
end

% q = round(1000/size(All_overlaps,1));
% all_comp = size(All_overlaps(All_overlaps>0),1);
match = zeros(1,8);
if Seq_Par > 1
    parpool(Par_workers)
    parfor SourceID=1:size(All_overlaps,1)
        %         [~,SseedsFileName]=fileparts(char(StackList(SourceID,1)));
        %         if isempty(SseedsFileName)
        %             SseedsFileName = num2str(SourceID);
        %         end
        %         Source_seed = load([featuresFolder,SseedsFileName, '_seeds.mat' ]);
        Source_seed_r_seed = hdf5read([DataFolder,'\tmp\Feature_seeds',num2str(SourceID),'.h5'], '/dataset1');
        Source_seed_FeatureVector = hdf5read([DataFolder,'\tmp\Feature_vector',num2str(SourceID),'.h5'], '/dataset1');
        if ~isempty(Source_seed_r_seed)
            j_ind=find(All_overlaps(SourceID,:));
            for jj=1:length(j_ind)
                TargetID=j_ind(jj);
                
                %                 [~,TseedsFileName]=fileparts(char(StackList(TargetID,1)));
                %                 if isempty(TseedsFileName)
                %                     TseedsFileName = num2str(TargetID);
                %                 end
                %                 Target_seed = load([featuresFolder,TseedsFileName, '_seeds.mat' ]);
                Target_seed_r_seed = hdf5read([DataFolder,'\tmp\Feature_seeds',num2str(TargetID),'.h5'], '/dataset1');
                Target_seed_FeatureVector = hdf5read([DataFolder,'\tmp\Feature_vector',num2str(TargetID),'.h5'], '/dataset1');
                
                if ~isempty(Target_seed_r_seed)
                    
                    
                    [SourceSubFeatures,SourceSubFeaturesVectors] = OverlapRegion1(SourceID,TargetID,Source_seed_r_seed,Source_seed_FeatureVector,StackPositions_pixels,StackSizes_pixels,paramsFMPad);
                    [TargetSubFeatures,TargetSubFeaturesVectors] = OverlapRegion1(TargetID,SourceID,Target_seed_r_seed,Target_seed_FeatureVector,StackPositions_pixels,StackSizes_pixels,paramsFMPad);
                    
                    
                    if size(SourceSubFeatures,1) >= paramsFMminReqFeatures && size(TargetSubFeatures,1)>= paramsFMminReqFeatures%params.FM.minReqFeatures
                        
                        Source_StackPositions = StackPositions_pixels(SourceID,:);
                        Target_StackPositions = StackPositions_pixels(TargetID,:);
                        [~,MatchLocations,~,~,~,~,~] = Stitching_3D_Func(StackSizes_pixels(SourceID,:),StackSizes_pixels(TargetID,:),listboxItems,v,StackList,SourceID,TargetID,SourceSubFeatures,SourceSubFeaturesVectors,TargetSubFeatures,TargetSubFeaturesVectors,Source_StackPositions,Target_StackPositions,TransformationValue,Seq_Par,tb11,stop,debug,DataFolder);
                        
                        sourceCol = zeros(size(MatchLocations,1),1);
                        sourceCol(:,1) = SourceID;
                        TargetCol = zeros(size(MatchLocations,1),1);
                        TargetCol(:,1) = TargetID;
                        
                        match = [match;[MatchLocations,sourceCol,TargetCol]];
                    end
                end
            end
        end
        
    end
    parfor_progress(0);
    delete(gcp)
    
    for SourceID=1:size(All_overlaps,1)
        j_ind=find(All_overlaps(SourceID,:));
        for jj=1:length(j_ind)
            TargetID=j_ind(jj);
            Matched(SourceID,TargetID) = {match(find(match(:,7) == SourceID & match(:,8) == TargetID),1:6)};
        end
    end
else
    for SourceID=1:size(All_overlaps,1) %  [106,107]%
        if get(tb11,'userdata') || stop% stop condition
            disp(num2str(tb11.UserData));
            stop = 1;
            listboxItems{v}  = 'Process Stopped ';
            v = v + 1;
            tb = findobj(NCT_Registration,'Tag', 'listbox1');
            set(tb, 'String', listboxItems);drawnow
            tb.Value = v-1;drawnow
            break;
        end
        %                    SourceID = 3% For Testing
        %         [~,SseedsFileName]=fileparts(char(StackList(SourceID,1)));
        %         if isempty(SseedsFileName)
        %             SseedsFileName = num2str(SourceID);
        %         end
        %         Source_seed = load([featuresFolder,SseedsFileName, '_seeds.mat' ]);
        Source_seed_r_seed = hdf5read([DataFolder,'\tmp\Feature_seeds',num2str(SourceID),'.h5'], '/dataset1');
        Source_seed_FeatureVector = hdf5read([DataFolder,'\tmp\Feature_vector',num2str(SourceID),'.h5'], '/dataset1');
        if ~isempty(Source_seed_r_seed)
            j_ind=find(All_overlaps(SourceID,:));
            for jj=1:length(j_ind)
                
                
                if get(tb11,'userdata') || stop% stop condition
                    disp(num2str(tb11.UserData));
                    stop = 1;
                    listboxItems{v}  = 'Process Stopped ';
                    v = v + 1;
                    tb = findobj(NCT_Registration,'Tag', 'listbox1');
                    set(tb, 'String', listboxItems);drawnow
                    tb.Value = v-1;drawnow
                    break;
                end
                
                
                TargetID=j_ind(jj);
                %              TargetID = 4% For Testing
                %                 [~,TseedsFileName]=fileparts(char(StackList(TargetID,1)));
                %                 if isempty(TseedsFileName)
                %                     TseedsFileName = num2str(TargetID);
                %                 end
                %                 Target_seed = load([featuresFolder,TseedsFileName, '_seeds.mat' ]);
                Target_seed_r_seed = hdf5read([DataFolder,'\tmp\Feature_seeds',num2str(TargetID),'.h5'], '/dataset1');
                Target_seed_FeatureVector = hdf5read([DataFolder,'\tmp\Feature_vector',num2str(TargetID),'.h5'], '/dataset1');
                if ~isempty(Target_seed_r_seed)
                    
                    [SourceSubFeatures,SourceSubFeaturesVectors] = OverlapRegion1(SourceID,TargetID,Source_seed_r_seed,Source_seed_FeatureVector,StackPositions_pixels,StackSizes_pixels,paramsFMPad);
                    [TargetSubFeatures,TargetSubFeaturesVectors] = OverlapRegion1(TargetID,SourceID,Target_seed_r_seed,Target_seed_FeatureVector,StackPositions_pixels,StackSizes_pixels,paramsFMPad);
                    
                    FeatureMatchingLog{k,1} = SourceID;
                    FeatureMatchingLog{k,2} = TargetID;
                    
                    if size(SourceSubFeatures,1) >= paramsFMminReqFeatures && size(TargetSubFeatures,1)>= paramsFMminReqFeatures%params.FM.minReqFeatures
                        tb = findobj(NCT_Registration,'Tag', 'listbox1');
                        listboxItems{v}  = ['Finding Correspondence between: ',num2str(SourceID),' and ',num2str(TargetID)];
                        set(tb, 'String', listboxItems);drawnow
                        v = v + 1;
                        tb.Value = v-1;drawnow
                        Source_StackPositions = StackPositions_pixels(SourceID,:);
                        Target_StackPositions = StackPositions_pixels(TargetID,:);
                        [Registrationtime,MatchLocations,listboxItems,v,~,~,stop] = Stitching_3D_Func(StackSizes_pixels(SourceID,:),StackSizes_pixels(TargetID,:),listboxItems,v,StackList,SourceID,TargetID,SourceSubFeatures,SourceSubFeaturesVectors,TargetSubFeatures,TargetSubFeaturesVectors,Source_StackPositions,Target_StackPositions,TransformationValue,Seq_Par,tb11,stop,debug,DataFolder);
                        
                        Matched(SourceID,TargetID) = {MatchLocations};
                        %                     end
                        
                    end
                    FeatureMatchingLog{k,3} = Registrationtime;
                    k = k + 1;
                    
                    drawnow;
                end
            end
        end
    end
end
TotalTime = toc

if TransformationValue ==1
    save(matchedLocationFile_Translation,'Matched','-v7.3');
elseif TransformationValue ==2
    save(matchedLocationFile_Rigid,'Matched','-v7.3');
elseif TransformationValue ==3
    save(matchedLocationFile_Affine,'Matched','-v7.3');
elseif TransformationValue ==4
    save(matchedLocationFile_Non_Rigid,'Matched','-v7.3');
end
tb = findobj(NCT_Registration,'Tag', 'listbox1');
listboxItems{v}  = ['Matched Locations are saved as: ',DataFolder];
set(tb, 'String', listboxItems);drawnow
v = v + 1;
tb.Value = v-1;drawnow
disp(['Matched Locations saved as: ',DataFolder]);
if ~isempty(tb9)
    tb9.Tag='axes3';
end

end

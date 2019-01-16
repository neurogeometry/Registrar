function [stop]=FeatureExtraction(LogHandle,StackList,DataFolder,Seq_Par,Par_workers,All_overlaps,StackPositions_pixels,StackSizes_pixels,debug)
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

mkdir([DataFolder,'/tmp/']);
if LogHandle ~= 0
    AllGUI = findobj(Registrar);
    %     log_ctrl = findall(AllGUI,'Tag','listbox1');
    %     log_ctrl = log_ctrl(1);
    % for Stop and Pause Button
    tb11 = findobj(Registrar,'Tag', 'pushbutton10');
    set(tb11,'userdata',0);
    x = 0:0.1:100;
    tb9 = findobj(Registrar,'Tag', 'axes3');
%     cla(tb9,'reset')
%     tb9.XAxis.Visible='off';
%     tb9.YAxis.Visible='off';
    if ~isempty(tb9)
        tb9.XLim = [0 100];
        tb9.YTickLabel  = [];
    end
    patch('XData',[0,0,100,100],'YData',[0,20,20,0],'FaceColor','white','Parent',tb9);
    patch('XData',[0,0,x(40),x(40)],'YData',[0,20,20,0],'FaceColor','green','Parent',tb9);
    drawnow;
    hold on %hold('units','on')
    if ~isempty(tb9)
        tb9.Tag='axes3';
    end
else
    tb11 = 0;
end
stop = 0;


parfor_progress(size(StackList,1),100);
if Seq_Par > 1 % do parallel
    parpool(Par_workers)
    parfor i=1:size(StackList,1)
        parfor_progress;
        tifFile = StackList(i,1);
        overlap_ind=[i,find(All_overlaps(i,:)),find(All_overlaps(:,i))'];
        [ImportTime,FeatureExtractionTime,numberofFeatures,seedsFile,~] = FeatureExtractionFunc(LogHandle,tifFile,i,tb11,stop,debug,DataFolder,StackPositions_pixels(overlap_ind,:),StackSizes_pixels(overlap_ind,:));
        
    end
    parfor_progress(0);
    delete(gcp)
else % do sequential
    for i=1:size(StackList,1)
        TifFileExist = 1;%any(size(dir([Folders{i} '/*.tif' ]),1));
        if TifFileExist
            tifFile = StackList(i,1);
            
            if LogHandle ~= 0
                LogHandle.Children(2).String{end+1} = ['Extracting Features for ',char(tifFile)];
                LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);
                %                 listboxItems{v}  = ['Extracting Features for ',char(tifFile)];
                %                 v = v + 1;
                
                %                 log_ctrl.String = listboxItems;drawnow
                %                 tb.Value = v-1;drawnow
                
                
                if get(tb11,'userdata') || stop% stop condition
                    disp(num2str(tb11.UserData));
                    stop = 1;
                    LogHandle.Children(2).String{end+1} = 'Process Stopped';
                    LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);
                    %                     listboxItems{v}  = 'Process Stopped';
                    %                     v = v + 1;
                    %                 tb = findobj(NCT_Registration,'Tag', 'listbox1');
                    %                     log_ctrl.String = listboxItems;drawnow
                    %                 set(tb, 'String', listboxItems);drawnow
                    %                     tb.Value = v-1;drawnow
                    
                    break;
                end
            else
                tb11 = 0;
            end
            overlap_ind=[i,find(All_overlaps(i,:)),find(All_overlaps(:,i))'];
            [ImportTime,FeatureExtractionTime,numberofFeatures,seedsFile,stop] = FeatureExtractionFunc(LogHandle,tifFile,i,tb11,stop,debug,DataFolder,StackPositions_pixels(overlap_ind,:),StackSizes_pixels(overlap_ind,:));
            if LogHandle ~= 0
                q = round(i/size(StackList,1)*500);
                patch('XData',[0,0,x(min(q,size(x,2))),x(min(q,size(x,2)))],'YData',[0,20,20,0],'FaceColor','green','Parent',tb9);
                drawnow;
            end
            
        end
    end
end

if LogHandle ~= 0
    if ~isempty(tb9)
        tb9.Tag='axes3';
    end
end

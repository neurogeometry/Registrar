function [stop]=FeatureExtraction(handles,LogHandle,StackList,DataFolder,Seq_Par,Par_workers,All_overlaps,StackPositions_pixels,StackSizes_pixels)
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

% Create Features folder
mkdir([DataFolder,'/tmp/']);

% for Stop
handles.pushbutton10.UserData = 0;
stop = 0;

% Progressbar setting
x = 0:0.1:100;
if ~isempty(handles.axes3)
    handles.axes3.XLim = [0 100];
    handles.axes3.YTickLabel  = [];
end
patch('XData',[0,0,100,100],'YData',[0,20,20,0],'FaceColor','white','Parent',handles.axes3);
patch('XData',[0,0,x(40),x(40)],'YData',[0,20,20,0],'FaceColor','green','Parent',handles.axes3);
drawnow;
hold on %hold('units','on')
if ~isempty(handles.axes3)
    handles.axes3.Tag='axes3';
end

parfor_progress(size(StackList,1),100);
if Seq_Par > 1 % do parallel
    parpool(Par_workers)
    parfor i=1:size(StackList,1)
        parfor_progress;
        tifFile = StackList(i,1);
        overlap_ind=[i,find(All_overlaps(i,:)),find(All_overlaps(:,i))'];
        [~] = FeatureExtractionFunc(handles,LogHandle,tifFile,i,stop,DataFolder,StackPositions_pixels(overlap_ind,:),StackSizes_pixels(overlap_ind,:));
        
    end
    parfor_progress(0);
    delete(gcp)
else % do sequential
    for i=1:size(StackList,1)
            tifFile = StackList(i,1);
            
            % Update Log
            if handles.checkbox15.Value
                UpdateLog(LogHandle,['Extracting Features for ',char(tifFile)]);
            end

            % Stop Function
            if handles.pushbutton10.UserData || stop% stop condition
                handles.pushbutton10.UserData = 1;
                disp(num2str(tb11.UserData));
                stop = 1;
                if handles.checkbox15.Value
                    UpdateLog(LogHandle,'Process Stopped');
                end
                break;
            else
                handles.pushbutton10.UserData = 0;
            end
            
            % Find overlap indexes
            overlap_ind=[i,find(All_overlaps(i,:)),find(All_overlaps(:,i))'];
            
            % Run Feature Extraction for one Stack
            [stop] = FeatureExtractionFunc(handles,LogHandle,tifFile,i,stop,DataFolder,StackPositions_pixels(overlap_ind,:),StackSizes_pixels(overlap_ind,:));
            
            % Update Progressbar
            q = round(i/size(StackList,1)*500);
            patch('XData',[0,0,x(min(q,size(x,2))),x(min(q,size(x,2)))],'YData',[0,20,20,0],'FaceColor','green','Parent',handles.axes3);
            drawnow;
    end
end

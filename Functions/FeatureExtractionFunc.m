function [stop] = FeatureExtractionFunc(handles,LogHandle,File,stackID,stop,DataFolder,StackPositions_pixels,StackSizes_pixels)
% ============================== About ====================================
% -------------------------------------------------------------------------
% Purpose: Feature Extracion Function
% Input:
%   Folder                  : 1*1, The folder address of stack
%   File                    : 1*1, The file name of stack
%   featuresFolder          : 1*1, The address to save the features
%
%  Output:
%   ImportTime              : 1*1, elapsed time to import stack
%   FeatureExtractionTime   : 1*1, elapsed time to extract the features
%   numberofFeatures        : 1*1, the number of extracted features
% -------------------------------------------------------------------------
% Author: Seyed Mostafa Mousavi Kahaki, Armen Stepanyants
% Northeastern University, USA
% =========================================================================
% -------------------------------------------------------------------------

% Load Parameters
parameters;

% Import Stack
tic
IM_Original=ImportStack(char(File),StackSizes_pixels(1,:));
if paramsREuseHDF5
    hdf5write([DataFolder,'/tmp/temp_',num2str(stackID),'.h5'], '/dataset1', IM_Original);
end
ImportTime = toc;

% Update Log
if handles.checkbox15.Value
    UpdateLog(LogHandle,['End Reading the file ',char(File),' - Import Time:',num2str(ImportTime)]);
end

% Find seeds
tic
r_seed=Find_Seeds(IM_Original,StackPositions_pixels,StackSizes_pixels);
SeedsTime = toc;
% Update Log
if handles.checkbox15.Value
    UpdateLog(LogHandle,['End Finding the Seeds for ',char(File),' - Finding Seeds Time:',num2str(SeedsTime)]);
end

% Remove extra features
[X1,Y1,Z1]=size(IM_Original);
if Z1 < params.FE.windowsize(3)
    params.FE.windowsize(3) = 0;
end
rem_ind_1 = (r_seed(:,1) <=params.FE.windowsize(1) | r_seed(:,1) > (X1-params.FE.windowsize(1)) | r_seed(:,2) <=params.FE.windowsize(2) | r_seed(:,2) > (Y1-params.FE.windowsize(2)) | r_seed(:,3) <=params.FE.windowsize(3) | r_seed(:,3) > (Z1-params.FE.windowsize(3)));
r_seed(rem_ind_1,:) = [];

% Visualize features on Image
if handles.chkdebug.Value
    try
        DebugHandle=findobj(0,'Name','Debug');
        tb1 = DebugHandle.Children(1);
        set(tb1, 'visible', 'on');
        h_im=imshow(max(IM_Original,[],3),[0 max(IM_Original(:))],'Parent',tb1);hold(tb1,'on')
        tb1=h_im.Parent;
        tb1.Tag='axes1';
        tb1.XLim = [0.5 size(h_im.CData,2)+0.5];
        tb1.YLim = [0.5 size(h_im.CData,1)+0.5];
        h_plot=findobj(tb1,'Tag','Fpoints');
        if ~isempty(h_plot)
            delete(h_plot)
            drawnow;
        end
        title([num2str(size(r_seed,1)),' Features are Detected for: ',char(File)],'Parent',tb1);hold on
        plot(r_seed(:,2),r_seed(:,1),'r*','Parent',tb1,'Tag','Fpoints');hold on
        drawnow;
    catch
    end
end

% Update Log
if handles.checkbox15.Value
    UpdateLog(LogHandle,['End of Extracting Features for the file ',char(File),' - ',num2str(size(r_seed,1)), ' Features are Detected']);
end

% Get Feature vectors
FeatureVector=[];
tic
FeatureVector = zeros(size(r_seed,1),prod(2*params.FE.windowsize+1));
for i=1:size(r_seed,1)
    temp = double(IM_Original(r_seed(i,1)-params.FE.windowsize(1):r_seed(i,1)+params.FE.windowsize(1),r_seed(i,2)-params.FE.windowsize(2):r_seed(i,2)+params.FE.windowsize(2),r_seed(i,3)-params.FE.windowsize(3):r_seed(i,3)+params.FE.windowsize(3)));
    
    % Check for Stop
    if handles.pushbutton10.UserData || stop% stop condition
        handles.pushbutton10.UserData = 1;
        disp(num2str(handles.pushbutton10.UserData));
        stop = 1;
        
        % Update Log
        if handles.checkbox15.Value
            UpdateLog(LogHandle,'Process Stopped');
        end
        break;
    else
        handles.pushbutton10.UserData = 0;
    end
    FeatureVector(i,:) =  temp(:)';
end
NeighborExtractionTime = toc;
% Update Log
if handles.checkbox15.Value
    UpdateLog(LogHandle,['End of Reading the Features Neighbor of the file ',char(File),' - extracting neighbors time: ',num2str(NeighborExtractionTime) ]);
end

% Update Log
numberofFeatures = size(r_seed);
if handles.checkbox15.Value
    UpdateLog(LogHandle,[num2str(numberofFeatures(1)),' Features are extracted for the file ',char(File)]);
end

% Write feature locations and Feature vectors into HDF5 files
hdf5write([DataFolder,'/tmp/Feature_seeds',num2str(stackID),'.h5'], '/dataset1', r_seed);
hdf5write([DataFolder,'/tmp/Feature_vector',num2str(stackID),'.h5'], '/dataset1', FeatureVector);
seedsFile = [DataFolder,'/tmp/Feature_',num2str(stackID),'.h5'];

% Update Log
if handles.checkbox15.Value
    UpdateLog(LogHandle,['Seeds are saved as: ',seedsFile]);
end
end
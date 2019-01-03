function [ImportTime,FeatureExtractionTime,numberofFeatures,seedsFile,stop] = FeatureExtractionFunc(LogHandle,File,stackID,tb11,stop,debug,DataFolder,StackPositions_pixels,StackSizes_pixels)
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

ImportTime = 0;
numberofFeatures = 0;
FeatureExtractionTime =0;
parameters;
tic

tic

IM_Original=ImportStack(char(File),StackSizes_pixels(1,:));
if paramsREuseHDF5
    hdf5write([DataFolder,'/tmp/temp_',num2str(stackID),'.h5'], '/dataset1', IM_Original);
end

ImportTime = toc;
if LogHandle ~= 0
    LogHandle.Children(2).String{end+1} = ['End Reading the file ',char(File),' - Import Time:',num2str(ImportTime)];
    LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);
    
    %     listboxItems{v}  = ['End Reading the file ',char(File),' - Import Time:',num2str(ImportTime)];
    %     v = v + 1;
    %     set(tb, 'String', listboxItems);drawnow
    %     tb.Value = v-1;drawnow
end
tic
r_seed=Find_Seeds(IM_Original,StackPositions_pixels,StackSizes_pixels);
SeedsTime = toc;
if LogHandle ~= 0
    LogHandle.Children(2).String{end+1} = ['End Finding the Seeds for ',char(File),' - Finding Seeds Time:',num2str(SeedsTime)];
    LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);
    
    %     listboxItems{v}  = ['End Finding the Seeds for ',char(File),' - Finding Seeds Time:',num2str(SeedsTime)];
    %     v = v + 1;
    %     set(tb, 'String', listboxItems);drawnow
    %     tb.Value = v-1;drawnow
end

[X1,Y1,Z1]=size(IM_Original);
if Z1 < params.FE.windowsize(3)
    params.FE.windowsize(3) = 0;
end

rem_ind_1 = (r_seed(:,1) <=params.FE.windowsize(1) | r_seed(:,1) > (X1-params.FE.windowsize(1)) | r_seed(:,2) <=params.FE.windowsize(2) | r_seed(:,2) > (Y1-params.FE.windowsize(2)) | r_seed(:,3) <=params.FE.windowsize(3) | r_seed(:,3) > (Z1-params.FE.windowsize(3)));
r_seed(rem_ind_1,:) = [];

if debug
    tb1 = findobj(NCT_Registration,'Tag', 'axes1');
    
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
    
end
if LogHandle ~= 0
    LogHandle.Children(2).String{end+1} = ['End of Extracting Features for the file ',char(File),' - ',num2str(size(r_seed,1)), ' Features are Detected'];
    LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);
    
    %     listboxItems{v}  = ['End of Extracting Features for the file ',char(File),' - ',num2str(size(r_seed,1)), ' Features are Detected'];
    %     v = v + 1;
    %     set(tb, 'String', listboxItems);drawnow
    %     tb.Value = v-1;drawnow
end
FeatureVector=[];

tic
FeatureVector = zeros(size(r_seed,1),prod(2*params.FE.windowsize+1));
for i=1:size(r_seed,1)
    temp = double(IM_Original(r_seed(i,1)-params.FE.windowsize(1):r_seed(i,1)+params.FE.windowsize(1),r_seed(i,2)-params.FE.windowsize(2):r_seed(i,2)+params.FE.windowsize(2),r_seed(i,3)-params.FE.windowsize(3):r_seed(i,3)+params.FE.windowsize(3)));
    if LogHandle ~= 0
        if get(tb11,'userdata') || stop% stop condition
            disp(num2str(tb11.UserData));
            stop = 1;
            LogHandle.Children(2).String{end+1} = 'Process Stopped';
            LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);
            %             listboxItems{v}  = 'Process Stopped ';
            %             v = v + 1;
            %             tb = findobj(NCT_Registration,'Tag', 'listbox1');
            %             set(tb, 'String', listboxItems);drawnow
            %             tb.Value = v-1;drawnow
            break;
        end
    end
    FeatureVector(i,:) =  temp(:)';
end

NeighborExtractionTime = toc;
if LogHandle ~= 0
    LogHandle.Children(2).String{end+1} = ['End of Reading the Features Neighbor of the file ',char(File),' - extracting neighbors time: ',num2str(NeighborExtractionTime) ];
    LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);
    
    %     listboxItems{v}  = ['End of Reading the Features Neighbor of the file ',char(File),' - extracting neighbors time: ',num2str(NeighborExtractionTime) ];
    %     v = v + 1;
    %     set(tb, 'String', listboxItems);drawnow
    %     tb.Value = v-1;drawnow
end
FeatureExtractionTime = toc;
numberofFeatures = size(r_seed);
if LogHandle ~= 0
    LogHandle.Children(2).String{end+1} = [num2str(numberofFeatures(1)),' Features are extracted for the file ',char(File)];
    LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);
    %     listboxItems{v}  = [num2str(numberofFeatures(1)),' Features are extracted for the file ',char(File)];
    %     v = v + 1;
    %     set(tb, 'String', listboxItems);drawnow
    %     tb.Value = v-1;drawnow
end

hdf5write([DataFolder,'/tmp/Feature_seeds',num2str(stackID),'.h5'], '/dataset1', r_seed);
hdf5write([DataFolder,'/tmp/Feature_vector',num2str(stackID),'.h5'], '/dataset1', FeatureVector);
seedsFile = [DataFolder,'/tmp/Feature_',num2str(stackID),'.h5'];

if LogHandle ~= 0
    LogHandle.Children(2).String{end+1} = ['Seeds are saved as: ',seedsFile];
    LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);
    
    %     listboxItems{v}  = ['Seeds are saved as: ',seedsFile];
    %     v = v + 1;
    %     set(tb, 'String', listboxItems);drawnow
    %     tb.Value = v-1;drawnow
end
end
function [ImportTime,FeatureExtractionTime,numberofFeatures,seedsFile,listboxItems,v,stop] = FeatureExtractionFunc(v,File,stackID,listboxItems,tb11,stop,debug,DataFolder,StackPositions_pixels,StackSizes_pixels)
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
if v ~= 0
    tb = findobj(NCT_Registration,'Tag', 'listbox1');
end
tic


IM_Original=ImportStack(char(File),StackSizes_pixels(1,:));
if paramsREuseHDF5
    tic,
    hdf5write([DataFolder,'\tmp\temp_',num2str(stackID),'.h5'], '/dataset1', IM_Original);
    toc
    %
    %     tic,
    %     X = hdf5read([DataFolder,'\tmp\temp_',num2str(stackID),'.h5'], '/dataset1');
    %     toc
    %
    %     IM_Original2D = reshape(IM_Original,[11*512 7*512]);
    %     tic,
    %     imwrite(IM_Original2D,[DataFolder,'\test.jpg'],'jpg');
    %     toc
    %
    %     tic,imread([DataFolder,'\test.jpg']);toc
    
    
    
    
end

ImportTime = toc;
if v ~= 0
    listboxItems{v}  = ['End Reading the file ',char(File),' - Import Time:',num2str(ImportTime)];
    v = v + 1;
    set(tb, 'String', listboxItems);drawnow
    tb.Value = v-1;drawnow
end
tic
r_seed=Find_Seeds(IM_Original,StackPositions_pixels,StackSizes_pixels);
SeedsTime = toc;
if v ~= 0
    listboxItems{v}  = ['End Finding the Seeds for ',char(File),' - Finding Seeds Time:',num2str(SeedsTime)];
    v = v + 1;
    set(tb, 'String', listboxItems);drawnow
    tb.Value = v-1;drawnow
end

[X1,Y1,Z1]=size(IM_Original);
if Z1 < params.FE.windowsize(3)
    params.FE.windowsize(3) = 0;
end

rem_ind_1 = (r_seed(:,1) <=params.FE.windowsize(1) | r_seed(:,1) > (X1-params.FE.windowsize(1)) | r_seed(:,2) <=params.FE.windowsize(2) | r_seed(:,2) > (Y1-params.FE.windowsize(2)) | r_seed(:,3) <=params.FE.windowsize(3) | r_seed(:,3) > (Z1-params.FE.windowsize(3)));
r_seed(rem_ind_1,:) = [];

%     if ign oreCenter
%     T_eli = [X1 * el_perc,Y1 * el_perc,Z1 * el_perc];
%     Stack_center=size(IM_Original)./2;
%     rem_ind_2 = (r_seed(:,1) <=Stack_center(1)+T_eli(1) & r_seed(:,1) >=Stack_center(1)-T_eli(1) &...
%                  r_seed(:,2) <=Stack_center(2)+T_eli(2) & r_seed(:,2) >=Stack_center(2)-T_eli(2) &...
%                  r_seed(:,3) <=Stack_center(3)+T_eli(3) & r_seed(:,3) >=Stack_center(3)-T_eli(3));
%
% %     rem_ind_2 = sum((r_seed <=Stack_center+T_eli & r_seed >= Stack_center-T_eli),2);
% %     r_seed(~xor(rem_ind_2,2),:) = [];
%     r_seed(rem_ind_2,:) = [];
%     end


%         if ~strcmp(Dataset, 'Svoboda')
%
%             if strcmp(Dataset, 'Holtmaat_2_1')
%                 figure,imshow(max(IM_Original,[],3),[0 500])
%             else
%                 figure,imshow(max(IM_Original,[],3),[0 max(IM_Original(:))])
%             end
%             hold on
%             plot(r_seed(:,1),r_seed(:,2),'r*');
%         else
%             featuresFig = figure;imshow(max(im2double(IM_Original),[],3),[0 1]);
%             hold on
%             plot(r_seed(:,1),r_seed(:,2),'r*');
%         end
% Show Features
%         figure,imshow(max(IM_Original,[],3),[0 max(IM_Original(:))])
%         hold on;plot(r_seed(:,2),r_seed(:,1),'r*');
% -----------------------

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
if v ~= 0
    listboxItems{v}  = ['End of Extracting Features for the file ',char(File),' - ',num2str(size(r_seed,1)), ' Features are Detected'];
    v = v + 1;
    set(tb, 'String', listboxItems);drawnow
    tb.Value = v-1;drawnow
end
FeatureVector=[];%zeros(size(r_seed,1),(2*w+1)^3);

tic
FeatureVector = zeros(size(r_seed,1),prod(2*params.FE.windowsize+1));
for i=1:size(r_seed,1)
    temp = double(IM_Original(r_seed(i,1)-params.FE.windowsize(1):r_seed(i,1)+params.FE.windowsize(1),r_seed(i,2)-params.FE.windowsize(2):r_seed(i,2)+params.FE.windowsize(2),r_seed(i,3)-params.FE.windowsize(3):r_seed(i,3)+params.FE.windowsize(3)));
    if v ~= 0
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
    end
    %Hessian Matrix
    if paramsFMEigen
        w=params.FE.windowsize;
        Hxx=temp(w(1)+1,w(2),w(3))+temp(w(1)-1,w(2),w(3))-2.*temp(w(1),w(2),w(3));
        Hyy=temp(w(1),w(2)+1,w(3))+temp(w(1),w(2)-1,w(3))-2.*temp(w(1),w(2),w(3));
        Hzz=temp(w(1),w(2),w(3)+1)+temp(w(1),w(2),w(3)-1)-2.*temp(w(1),w(2),w(3));
        Hxy=(temp(w(1)+1,w(2)+1,w(3))+temp(w(1)-1,w(2)-1,w(3))-temp(w(1)+1,w(2)-1,w(3))-temp(w(1)-1,w(2)+1,w(3)))./4;
        Hxz=(temp(w(1)+1,w(2),w(3)+1)+temp(w(1)-1,w(2),w(3)-1)-temp(w(1)+1,w(2),w(3)-1)-temp(w(1)-1,w(2),w(3)+1))./4;
        Hyz=(temp(w(1),w(2)+1,w(3)+1)+temp(w(1),w(2)-1,w(3)-1)-temp(w(1),w(2)+1,w(3)-1)-temp(w(1),w(2)-1,w(3)+1))./4;
        H=[Hxx,Hxy,Hxz;Hxy,Hyy,Hyz;Hxz,Hyz,Hzz];
        [vec,l] = eig(H);
        FeatureVector =  [FeatureVector;vec(:)',l(:)'];
        %         [Lambda1i,Lambda2i,Lambda3i]=eig3volume(Hxx,Hxy,Hxz,Hyy,Hyz,Hzz);
        %             FeatureVector =  [FeatureVector;H(:)'];
        
    else
        FeatureVector(i,:) =  temp(:)';
    end
    
end

%     % GPU
%     I = gpuArray(IM_Original);
%     for i=1:size(r_seed,1)
%         temp = double(I(r_seed(i,2)-w:r_seed(i,2)+w,r_seed(i,1)-w:r_seed(i,1)+w,r_seed(i,3)-w:r_seed(i,3)+w));
%         FeatureVector =  [FeatureVector;temp(:)'];
%     end

NeighborExtractionTime = toc;
if v ~= 0
    listboxItems{v}  = ['End of Reading the Features Neighbor of the file ',char(File),' - extracting neighbors time: ',num2str(NeighborExtractionTime) ];
    v = v + 1;
    set(tb, 'String', listboxItems);drawnow
    tb.Value = v-1;drawnow
end
FeatureExtractionTime = toc;
numberofFeatures = size(r_seed);
if v ~= 0
    listboxItems{v}  = [num2str(numberofFeatures(1)),' Features are extracted for the file ',char(File)];
    v = v + 1;
    set(tb, 'String', listboxItems);drawnow
    tb.Value = v-1;drawnow
end

% save([featuresFolder,File{1}(1:5),'_seeds'],'r_seed','FeatureVector','-v7.3')
% [PathStr,seedsFileName]=fileparts(File{1});
% if isempty(seedsFileName)
%     seedsFileName = num2str(stackID);
% end
% seedsFileName = regexp(File{1},'\\([^\\]*)$','tokens','once');
hdf5write([DataFolder,'\tmp\Feature_seeds',num2str(stackID),'.h5'], '/dataset1', r_seed);
hdf5write([DataFolder,'\tmp\Feature_vector',num2str(stackID),'.h5'], '/dataset1', FeatureVector);
seedsFile = [DataFolder,'\tmp\Feature_',num2str(stackID),'.h5'];
% save([featuresFolder,char(seedsFileName),'_seeds.mat'],'r_seed','FeatureVector','-v7.3')
% seedsFile = [featuresFolder,char(seedsFileName),'_seeds.mat'];
if v ~= 0
    listboxItems{v}  = ['Seeds are saved as: ',seedsFile];
    v = v + 1;
    set(tb, 'String', listboxItems);drawnow
    tb.Value = v-1;drawnow
end
end
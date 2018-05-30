function [Registrationtime,MatchLocations,listboxItems,v,Transformation_T,b,stop] = Stitching_3D_Func(SourceStackSize,TargetStackSize,listboxItems,v,StackList,SourceID,TargetID,Source_seed,SourceFeatures,Target_seed,TargetFeatures,Source_StackPositions,Target_StackPositions,TransformationValue,Seq_Par,tb11,stop,debug,DataFolder)
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
addpath('Functions');
parameters;
tic

Displacement = Source_StackPositions-Target_StackPositions;

if debug ==1 && Seq_Par ~= 2
    Source_Stack_File = StackList(SourceID,1);
    Target_Stack_File = StackList(TargetID,1);
    
    IM_Source=ImportStack(char(Source_Stack_File),SourceStackSize);
    IM_Target=ImportStack(char(Target_Stack_File),TargetStackSize);
    
    
    IM_source_max=max(IM_Source,[],3);
    IM_target_max=max(IM_Target,[],3);
    M_source=max(IM_Source(:));
    
    if abs(Displacement(2)) > abs(Displacement(1))
        Direction='horizontal';
    elseif abs(Displacement(1)) > abs(Displacement(2))
        Direction='vertical';
    else
        Direction='vertical';
    end
    %
    [X1,Y1,~]=size(IM_Source);
    [X2,Y2,~]=size(IM_Target);
    
end








hungInput = zeros(size(Source_seed,1),size(Target_seed,1));

% ZMAD = ZeroMeanAbsoluteDifferences
% MSE: MeanSquareError
% AD = AverageDifferences
% NAE = NormalizedAbsoluteError
% NCC= NormalizedCrossCorrelation
% ZNCC= Zero mean NormalizedCrossCorrelation
% SAD = Sum of Absolute Differences
% SAD =  Sum of Squared Differences
% RMSD = Root Mean Square Distance
% SSIM = StructuralSimilarity
% SSD = Sum of Square Differences (SSD)
if strcmp(paramsFMMetric,'ZMAD')
    nTargetFeatures=bsxfun(@rdivide,TargetFeatures,mean(TargetFeatures,2));
    nSourceFeatures=bsxfun(@rdivide,SourceFeatures,mean(SourceFeatures,2));
    for i=1:size(SourceFeatures,2)
        hungInput=hungInput+abs(bsxfun(@minus,nTargetFeatures(:,i)',nSourceFeatures(:,i)));
    end
    hungInput=hungInput./size(SourceFeatures,2);
    
    TempDisp=zeros(size(SourceFeatures,1),size(TargetFeatures,1));
    for i=1:size(Target_seed,2)
        TempDisp=TempDisp+(bsxfun(@minus,Target_seed(:,i)',Source_seed(:,i))).^2;
    end
    Dthr=(sum(Displacement.^2))^0.5;
    hungInput(isnan(hungInput)) = 10^12;
    hungInput(abs(TempDisp.^0.5-Dthr)>paramsFMDT)=10^12;
else
    
    for i = 1 : size(Source_seed,1)
        for j = 1 : size(Target_seed,1)
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
            TempDisp = sum((Target_seed(j,:)-Source_seed(i,:))-Displacement);
            
            if  TempDisp > paramsFMDT
                hungInput(i,j) = 10^12;
            else
                switch  paramsFMMetric
                    case 'ZMAD'
                        hungInput(i,j) = mean(abs(TargetFeatures(j,:)./mean(TargetFeatures(j,:))-SourceFeatures(i,:)./mean(SourceFeatures(i,:))));
                    case 'MSE'
                        %                     error = TargetFeatures(j,:) - SourceFeatures(i,:);
                        %                     hungInput(i,j) = sum(error .* error) / (size(TargetFeatures(j,:),1) * size(TargetFeatures(j,:),2));
                        hungInput(i,j) = immse(TargetFeatures(j,:),SourceFeatures(i,:));
                    case 'AD'
                        hungInput(i,j) = (sum(TargetFeatures(j,:) - SourceFeatures(i,:)))/(size(TargetFeatures(j,:),1) * size(TargetFeatures(j,:),2));
                    case 'NAE'
                        error = TargetFeatures(j,:) - SourceFeatures(i,:);
                        hungInput(i,j) = sum(abs(error))/ sum(TargetFeatures(j,:));
                    case 'NCC' % https://en.wikipedia.org/wiki/Cross-correlation
                        hungInput(i,j) = (sum(TargetFeatures(j,:) .* SourceFeatures(i,:)) / sum(std(TargetFeatures(j,:)) .* std(TargetFeatures(j,:))))/size(TargetFeatures(j,:),2);
                        %                     hungInput(i,j) = normxcorr2(TargetFeatures(j,:),SourceFeatures(i,:))
                    case 'ZNCC' % https://en.wikipedia.org/wiki/Cross-correlation
                        hungInput(i,j) = (sum((TargetFeatures(j,:)-mean(TargetFeatures(j,:))) .* (SourceFeatures(i,:)-mean(SourceFeatures(i,:)))) / sum(std(TargetFeatures(j,:)) .* std(TargetFeatures(j,:))))/size(TargetFeatures(j,:),2);
                    case 'SAD'
                        hungInput(i,j) = sum(abs(TargetFeatures(j,:)-SourceFeatures(i,:)));
                    case 'ZSAD'
                        hungInput(i,j) = sum(abs((TargetFeatures(j,:)-mean(TargetFeatures(j,:)))-(SourceFeatures(i,:)-mean(SourceFeatures(i,:)))));
                    otherwise % ZMAD
                        hungInput(i,j) = mean(abs(TargetFeatures(j,:)./mean(TargetFeatures(j,:))-SourceFeatures(i,:)./mean(SourceFeatures(i,:))));
                end
                
            end
        end
    end
    
end

Am = Hungarian_fast(hungInput);
Am(hungInput==10^12)=0;
[idx1,idx2]=find(Am);


if Seq_Par ~= 2 && v ~= 0
    tb = findobj(NCT_Registration,'Tag', 'listbox1');
    listboxItems{v}  = ['Number of Hungarian Matches:',num2str(size(idx1,1))];
    set(tb, 'String', listboxItems);drawnow
    v = v + 1;
    tb.Value = v-1;drawnow
end
x_Source = Source_seed(idx1,1);
x_Target = Target_seed(idx2,1);
y_Source = Source_seed(idx1,2);
y_Target = Target_seed(idx2,2);
z_Source = Source_seed(idx1,3);
z_Target = Target_seed(idx2,3);

if debug ==1 && Seq_Par ~= 2 && v ~= 0
    if strcmp(Direction,'horizontal')
        if Displacement(2)<0
            stitched = appendimages(IM_source_max,IM_target_max,Direction);
        else
            stitched = appendimages(IM_target_max,IM_source_max,Direction);
        end
    else
        if Displacement(1)<0
            stitched = appendimages(IM_source_max,IM_target_max,Direction);
        else
            stitched = appendimages(IM_target_max,IM_source_max,Direction);
        end
    end
end

matchLoc_Target = [x_Target,y_Target,z_Target];
matchLoc_Source = [x_Source,y_Source,z_Source];

i=1;
Match_Indexes=[];
Transformation_T = [];
b = [];
Global_Matched_Source = matchLoc_Source'+(ones(size(matchLoc_Source,1),1)*Source_StackPositions)'-1;
Global_Matched_Target = matchLoc_Target'+(ones(size(matchLoc_Source,1),1)*Target_StackPositions)'-1;
while numel(Match_Indexes)<paramsFMminmatches && i<paramsFMmaxiter
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
    Match_Indexes_temp = RANSAC(Global_Matched_Source,Global_Matched_Target,TransformationValue);
    if numel(Match_Indexes_temp)<paramsFMminmatches
        disp('Trying again iteration ');
        i=i+1;
    end
    if numel(Match_Indexes)<numel(Match_Indexes_temp)
        Match_Indexes=Match_Indexes_temp;
    end
end

toc

if Seq_Par ~= 2 && v ~= 0
    tb = findobj(NCT_Registration,'Tag', 'listbox1');
    listboxItems{v}  = ['Number of Final Matches:',num2str(size(Match_Indexes,2))];
    set(tb, 'String', listboxItems);drawnow
    v = v + 1;
    tb.Value = v-1;drawnow
end
if debug ==1 && Seq_Par ~= 2 && v ~= 0
    
    tb2 = findobj(NCT_Registration,'Tag', 'axes1');
    
    set(tb2, 'visible', 'on');
    h_plot=findobj(tb2,'Tag','Fpoints');
    h_plot3=findobj(tb2,'Tag','Flines');
    if ~isempty(h_plot)
        delete(h_plot)
        delete(h_plot3)
        drawnow;
    end
    h_im=imshow(stitched,[0 M_source],'Parent',tb2);hold(tb2,'on')
    tb2=h_im.Parent;
    tb2.Tag='axes1';
    tb2.XLim = [0.5 size(h_im.CData,2)+0.5];
    tb2.YLim = [0.5 size(h_im.CData,1)+0.5];
    title('Correct Matches','Parent',tb2);
    
    if strcmp(Direction,'horizontal')
        if Displacement(2)<0
            plot(matchLoc_Source(Match_Indexes,2),matchLoc_Source(Match_Indexes,1),'*','Parent',tb2,'Tag','Fpoints');
        else
            plot(matchLoc_Source(Match_Indexes,2)+Y1,matchLoc_Source(Match_Indexes,1),'*','Parent',tb2,'Tag','Fpoints');
        end
        if Displacement(2)<0
            plot(matchLoc_Target(Match_Indexes,2)+Y2,matchLoc_Target(Match_Indexes,1),'*','Parent',tb2,'Tag','Fpoints');
        else
            plot(matchLoc_Target(Match_Indexes,2),matchLoc_Target(Match_Indexes,1),'*','Parent',tb2,'Tag','Fpoints');
        end
        
    else
        if Displacement(1)<0
            plot(matchLoc_Source(Match_Indexes,2),matchLoc_Source(Match_Indexes,1),'*','Parent',tb2,'Tag','Fpoints');
        else
            plot(matchLoc_Source(Match_Indexes,2),matchLoc_Source(Match_Indexes,1)+X1,'*','Parent',tb2,'Tag','Fpoints');
        end
        if Displacement(1)<0
            plot(matchLoc_Target(Match_Indexes,2),matchLoc_Target(Match_Indexes,1)+X2,'*','Parent',tb2,'Tag','Fpoints');
        else
            plot(matchLoc_Target(Match_Indexes,2),matchLoc_Target(Match_Indexes,1),'*','Parent',tb2,'Tag','Fpoints');
        end
    end
    
    for i = 1: length(Match_Indexes)
        if strcmp(Direction,'horizontal')
            if Displacement(2)<0
                line([matchLoc_Source(Match_Indexes(i),2) matchLoc_Target(Match_Indexes(i),2)+Y1], ...
                    [matchLoc_Source(Match_Indexes(i),1) matchLoc_Target(Match_Indexes(i),1)], 'Color', rand(1,3),'Parent',tb2,'Tag','Flines');
            else
                line([matchLoc_Source(Match_Indexes(i),2)+Y1 matchLoc_Target(Match_Indexes(i),2)], ...
                    [matchLoc_Source(Match_Indexes(i),1) matchLoc_Target(Match_Indexes(i),1)], 'Color', rand(1,3),'Parent',tb2,'Tag','Flines');
            end
        else
            if Displacement(1)<0
                line([matchLoc_Source(Match_Indexes(i),2) matchLoc_Target(Match_Indexes(i),2)], ...
                    [matchLoc_Source(Match_Indexes(i),1) matchLoc_Target(Match_Indexes(i),1)+X1], 'Color', rand(1,3),'Parent',tb2,'Tag','Flines');
            else
                line([matchLoc_Source(Match_Indexes(i),2) matchLoc_Target(Match_Indexes(i),2)], ...
                    [matchLoc_Source(Match_Indexes(i),1)+X1 matchLoc_Target(Match_Indexes(i),1)], 'Color', rand(1,3),'Parent',tb2,'Tag','Flines');
            end
        end
        
    end
    
end

MatchLocations=NaN(length(Match_Indexes),6);
if ~isempty(Match_Indexes)
    for i = 1:length(Match_Indexes)
        MatchLocations(i,:) = [matchLoc_Source(Match_Indexes(i),:),matchLoc_Target(Match_Indexes(i),:)];
    end
end

Registrationtime=toc;
if Seq_Par ~= 2 && v ~= 0
    tb = findobj(NCT_Registration,'Tag', 'listbox1');
    listboxItems{v}  = ['Correspondence Finding Time:',num2str(Registrationtime)];
    set(tb, 'String', listboxItems);drawnow
    v = v + 1;
    tb.Value = v-1;drawnow
end

end
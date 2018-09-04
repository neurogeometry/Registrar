r_source=[];SourcePoints=[];AM_source=[];
            r_target=[];TargetPoints=[];AM_target=[];
            r_source = r_targetAll{sourceID-1};
            SourcePoints = TargetPointsAll{sourceID-1};
            AM_source = AM_targetAll{sourceID-1};
            r_target = r_targetAll{sourceID};
            TargetPoints = TargetPointsAll{sourceID};
            AM_target = AM_targetAll{sourceID};
            
        HC = [];
        for i=1:size(SourcePoints,1)
            for j=1:size(TargetPoints,1)
                HC(i,j) = mean(sum((SourcePoints(i,:)-TargetPoints(j,:)).^2,2).^0.5);
            end
        end
        
        
        
        
        
        
%%                     ----------------- Fig 4.B  - Points Correct - Lines Wrong
L=[];b=[];Cxyz=[];Grid_start=[];Nxyz=[];
IMAll = uint8(zeros(1024,1024));
TR = 0;
AddSpace = 40;
for ID = 1: size(StackList,1)
    sourceID = ID;
    targetID = ID + 1;
    
    Global_Matched_Source = Matched{sourceID,targetID}(:,4:6)';
    Global_Matched_Target = Matched{sourceID,targetID}(:,1:3)';
    [~,L{sourceID},b{sourceID},Cxyz{sourceID},Nxyz{sourceID},nxyz,Grid_start{sourceID}]=Optimal_Bspline_Transform(Global_Matched_Source,Global_Matched_Target,nxyz,affine,mu);
    
    fname = dir([GTpath,'Matches\Traces\DL083',T_Names{sourceID},'001-A0*']);
    fname={fname.name}';
    Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
    [AM_source,r_source,R]=swc2AM(Path);
    [AM_source,r_source,~] = AdjustPPM(AM_source,r_source,R,ppm);
    
    %         Stack_File = char(StackList(sourceID,1));
    %         IM=ImportStack(char(Stack_File));
    %         IM = uint8(double(IM)./double(max(IM(:))).*255);
    %         IM_max=max(IM,[],3);
    %
    %         [KT]=FastMarchingTube(size(IM),r_source,3,[1,1,1]);
    %         IMmax =  max(uint8(KT).*IM,[],3);
    %         IMmax = imtranslate(IMmax,[0, TR]);
    %         IMAll = max(IMAll,IMmax);
    %         figure(N);imshow(IMAll,[0 20])
    
    if ID < size(StackList,1)
        
        fname = dir([GTpath,'Matches\Traces\DL083',T_Names{targetID},'001-A0*']);
        fname={fname.name}';
        Path = [GTpath,'Matches\Traces\',fname{TraceNum}];
        [AM_target,r_target,R]=swc2AM(Path);
        [AM_target,r_target,~] = AdjustPPM(AM_target,r_target,R,ppm);
        
        B1 = Boutons{sourceID,1};
        SourcePoints = B1.r1;
        TargetPoints = B1.r2;
        
        if ID > 1
            %         sourcePoints_NR_temp = SourcePoints';
            %         r_source_NR_temp = r_source';
            %         for j=size(b,2):-1:1
            %             [sourcePoints_NR_temp,~]=Perform_Bspline_Transform(sourcePoints_NR_temp,[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
            %             [r_source_NR_temp,~]=Perform_Bspline_Transform(r_source_NR_temp,[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
            %         end
            %         SourcePoints = sourcePoints_NR_temp';
            %         r_source = r_source_NR_temp';
            r_source=[];SourcePoints=[];AM_source=[];
            r_source = r_targetAll{sourceID-1};
            SourcePoints = TargetPointsAll{sourceID-1};
            AM_source = AM_targetAll{sourceID-1};
        end
        
        targetPoints_NR_temp = TargetPoints';
        r_target_NR_temp = r_target';
        for j=size(b,2):-1:1
            [targetPoints_NR_temp,~]=Perform_Bspline_Transform(targetPoints_NR_temp,[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
            [r_target_NR_temp,~]=Perform_Bspline_Transform(r_target_NR_temp,[],L{j},b{j},Cxyz{j},Nxyz{j},nxyz,Grid_start{j},affine);
        end
        TargetPoints = targetPoints_NR_temp';
        r_target = r_target_NR_temp';
        r_targetAll{sourceID} = r_target;
        AM_targetAll{sourceID}=AM_target;
        TargetPointsAll{sourceID} = TargetPoints;
        
        [k,f] = dsearchn(r_source,SourcePoints);
        
        s = 1;
        SourcePoints_onTrace=[];
        for j=1:size(k)
            if f(j)< 3
                SourcePoints_onTrace(s,:) = SourcePoints(j,:);
                s = s +1;
            end
        end
        figure(N); hold on;
        %     if ID==1
        %     figure(N); hold on; PlotAM(AM_source,r_source,'c')
        %     end
        
        plot3(SourcePoints_onTrace(:,2),SourcePoints_onTrace(:,1)+TR,ones(1,length(SourcePoints_onTrace(:,2))),'or')
        
        TR = TR - AddSpace;
        [k,f] = dsearchn(r_target,TargetPoints);
        TargetPoints(:,1) = TargetPoints(:,1)+TR;
        r_target(:,1)=r_target(:,1)+TR;
        %     PlotAM(AM_target,r_target,'c')
        s = 1;
        TargetPoints_onTrace=[];
        for j=1:size(k)
            if f(j)< 3
                TargetPoints_onTrace(s,:) = TargetPoints(j,:);
                s = s +1;
            end
        end
        plot3(TargetPoints_onTrace(:,2),TargetPoints_onTrace(:,1),ones(1,length(TargetPoints_onTrace(:,2))),'or')
        
        %         if ID==1
        %             line([SourcePoints_onTrace(:,2) TargetPoints_onTrace(:,2)]', ...
        %                 [SourcePoints_onTrace(:,1) TargetPoints_onTrace(:,1)]',[1;1]*ones(1,length(TargetPoints_onTrace(:,1))), 'Color', 'c');
        %         else
        %             line([SourcePoints_onTrace(:,2) TargetPoints_onTrace(:,2)]', ...
        %                 [SourcePoints_onTrace(:,1)+TR+AddSpace TargetPoints_onTrace(:,1)]',[1;1]*ones(1,length(TargetPoints_onTrace(:,1))), 'Color', 'c');
        %         end
        SourcePoints_onTrace;
        TargetPoints_onTrace;
        HC = [];
        for i=1:size(SourcePoints_onTrace,1)
            for j=1:size(TargetPoints_onTrace,1)
                HC(i,j) = mean(sum((SourcePoints_onTrace(i,:)-TargetPoints_onTrace(j,:)).^2,2).^0.5);
            end
        end
        
        Am = Hungarian_fast(HC);
        [idx1,idx2]=find(Am);
        x_Source = SourcePoints_onTrace(idx1,1);
        x_Target = TargetPoints_onTrace(idx2,1);
        y_Source = SourcePoints_onTrace(idx1,2);
        y_Target = TargetPoints_onTrace(idx2,2);
        z_Source = SourcePoints_onTrace(idx1,3);
        z_Target = TargetPoints_onTrace(idx2,3);
        matchLoc_Source = [x_Target,y_Target,z_Target];
        matchLoc_Target = [x_Source,y_Source,z_Source];
        
        if ID==1
        line([matchLoc_Source(:,2) matchLoc_Target(:,2)]', ...
                [matchLoc_Source(:,1) matchLoc_Target(:,1)]',...
                [1;1]*ones(1,length(matchLoc_Source(:,1))),'Color', rand(1,3));
        else
         line([matchLoc_Source(:,2) matchLoc_Target(:,2)]', ...
                [matchLoc_Source(:,1) matchLoc_Target(:,1)+TR+AddSpace]',...
                [1;1]*ones(1,length(matchLoc_Source(:,1))),'Color', rand(1,3));  
        end
        axis equal
        
        
        
        
        
        
        
        drawnow
    end
    
end
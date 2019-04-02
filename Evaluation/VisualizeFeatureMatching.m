function VisualizeFeatureMatching(Matches)%(PathSource,PathTarget,Matches)
PathSource = '/Users/chuhan/Documents/study/2019Spring/neurogeometry/NCTWeb/NCT-Web/src/main/data/registration/Images/01/00001.tif';
PathTarget = '/Users/chuhan/Documents/study/2019Spring/neurogeometry/NCTWeb/NCT-Web/src/main/data/registration/Images/02/00002.tif';
Source_StackPositions = [1851,749,-22];
Target_StackPositions = [1857,1192,-7];
matchLoc_Source = Matches(:,1:3);
matchLoc_Target = Matches(:,4:6);
Match_Indexes = 1:size(Matches,1);
IMSource=ImportStack(PathSource);
IMTarget=ImportStack(PathTarget);
IM_source_max=max(IMSource,[],3);
IM_target_max=max(IMTarget,[],3);
[X1,Y1,~]=size(IMSource);
[X2,Y2,~]=size(IM_target_max);
M_source=max(IMSource(:));

Displacement = Source_StackPositions-Target_StackPositions;
Direction='vertical';

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
figure, imshow(stitched,[0 M_source]);hold on;
if strcmp(Direction,'horizontal')
    if Displacement(2)<0
        plot(matchLoc_Source(:,2),matchLoc_Source(Match_Indexes,1),'*');
    else
        plot(matchLoc_Source(Match_Indexes,2)+Y1,matchLoc_Source(Match_Indexes,1),'*');
    end
    if Displacement(2)<0
        plot(matchLoc_Target(Match_Indexes,2)+Y2,matchLoc_Target(Match_Indexes,1),'*');
    else
        plot(matchLoc_Target(Match_Indexes,2),matchLoc_Target(Match_Indexes,1),'*');
    end
    
else
    if Displacement(1)<0
        plot(matchLoc_Source(Match_Indexes,2),matchLoc_Source(Match_Indexes,1),'*');
    else
        plot(matchLoc_Source(Match_Indexes,2),matchLoc_Source(Match_Indexes,1)+X1,'*');
    end
    if Displacement(1)<0
        plot(matchLoc_Target(Match_Indexes,2),matchLoc_Target(Match_Indexes,1)+X2,'*');
    else
        plot(matchLoc_Target(Match_Indexes,2),matchLoc_Target(Match_Indexes,1),'*');
    end
end

for i = 1: length(Match_Indexes)
    if strcmp(Direction,'horizontal')
        if Displacement(2)<0
            line([matchLoc_Source(Match_Indexes(i),2) matchLoc_Target(Match_Indexes(i),2)+Y1], ...
                [matchLoc_Source(Match_Indexes(i),1) matchLoc_Target(Match_Indexes(i),1)], 'Color', rand(1,3));
        else
            line([matchLoc_Source(Match_Indexes(i),2)+Y1 matchLoc_Target(Match_Indexes(i),2)], ...
                [matchLoc_Source(Match_Indexes(i),1) matchLoc_Target(Match_Indexes(i),1)], 'Color', rand(1,3));
        end
    else
        if Displacement(1)<0
            line([matchLoc_Source(Match_Indexes(i),2) matchLoc_Target(Match_Indexes(i),2)], ...
                [matchLoc_Source(Match_Indexes(i),1) matchLoc_Target(Match_Indexes(i),1)+X1], 'Color', rand(1,3));
        else
            line([matchLoc_Source(Match_Indexes(i),2) matchLoc_Target(Match_Indexes(i),2)], ...
                [matchLoc_Source(Match_Indexes(i),1)+X1 matchLoc_Target(Match_Indexes(i),1)], 'Color', rand(1,3));
        end
    end
    
end





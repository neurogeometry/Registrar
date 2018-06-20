function RegistrationVisualization(Sourcepath,Targetpath,Matched_Points,Direction)

IM_Source=ImportStack(Sourcepath);

IM_Target=ImportStack(Targetpath);

M_source=max(IM_Source(:));
IM_source_max=max(IM_Source,[],3);
IM_target_max=max(IM_Target,[],3);

[Y1,X1,Z1]=size(IM_Source);
[Y2,X2,Z2]=size(IM_Target);

stitched = appendimages(IM_source_max,IM_target_max,Direction);
    figure(20);imshow(stitched,[0 M_source]);
    hold on
    plot(Matched_Points(:,1),Matched_Points(:,2),'*');
    if strcmp(Direction,'horizontal')
        plot(Matched_Points(:,4)+X1,Matched_Points(:,5),'*');
    else
        plot(Matched_Points(:,4),Matched_Points(:,5)+Y1,'*');
    end
    figure(20);
    hold on
    for i = 1:size(Matched_Points,1)
        if strcmp(Direction,'horizontal')
            line([Matched_Points(i,1) Matched_Points(i,4)+X1],[Matched_Points(i,2) Matched_Points(i,5)], 'Color', rand(1,3));
        else
            line([Matched_Points(i,1) Matched_Points(i,4)],[Matched_Points(i,2) Matched_Points(i,5)+Y1], 'Color', rand(1,3));
        end
    end

end
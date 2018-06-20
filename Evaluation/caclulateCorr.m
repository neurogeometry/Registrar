function correl = caclulateCorr(StackPosition_Source,StackPosition_Target,IM_Source,IM_Target,pad,showIM)



temp_Source=IM_Source(pad(1)+max(1,StackPosition_Target(1)-StackPosition_Source(1)+1):-pad(1)+min(StackPosition_Target(1)-StackPosition_Source(1)+size(IM_Target,1),size(IM_Source,1)),...
        pad(2)+max(1,StackPosition_Target(2)-StackPosition_Source(2)+1):-pad(2)+min(StackPosition_Target(2)-StackPosition_Source(2)+size(IM_Target,2),size(IM_Source,2)),...
        pad(3)+max(1,StackPosition_Target(3)-StackPosition_Source(3)+1):-pad(3)+min(StackPosition_Target(3)-StackPosition_Source(3)+size(IM_Target,3),size(IM_Source,3)));
    
    
    temp_Target=IM_Target(pad(1)+max(1,StackPosition_Source(1)-StackPosition_Target(1)+1):-pad(1)+min(StackPosition_Source(1)-StackPosition_Target(1)+size(IM_Source,1),size(IM_Target,1)),...
        pad(2)+max(1,StackPosition_Source(2)-StackPosition_Target(2)+1):-pad(2)+min(StackPosition_Source(2)-StackPosition_Target(2)+size(IM_Source,2),size(IM_Target,2)),...
        pad(3)+max(1,StackPosition_Source(3)-StackPosition_Target(3)+1):-pad(3)+min(StackPosition_Source(3)-StackPosition_Target(3)+size(IM_Source,3),size(IM_Target,3)));
    
    if showIM
        IM_Source_max=max(IM_Source,[],3);
        IM_Target_max=max(IM_Target,[],3);
        figure;imshow(max(temp_Source,[],3),[0 max(IM_Source_max(:))]);
        figure;imshow(max(temp_Target,[],3),[0 max(IM_Target_max(:))]);
    end

    cor3 = corrcoef(double(temp_Target),double(temp_Source));





correl = cor3(1,2);




end
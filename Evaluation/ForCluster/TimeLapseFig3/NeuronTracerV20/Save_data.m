% Function needed for MatLab - Java data format conversion.

function Save_data(filepath,Original,IM,AM,r,R,reduction_x,reduction_y,reduction_z)

if ~isempty(AM)
    rem_ind=(sum(AM)==0);
    AM(rem_ind,:)=[];
    AM(:,rem_ind)=[];
    r(rem_ind,:)=[];
    R(rem_ind)=[];
else
    AM=[];
    r=[];
    R=[];
end

save(filepath,'Original','IM','AM','r','R','reduction_x','reduction_y','reduction_z')
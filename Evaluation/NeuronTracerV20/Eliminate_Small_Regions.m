% This function eliminates small (by number of voxels) regions from Im

function Im=Eliminate_Small_Regions(Im,size_thr)

disp('Small region elimination started.')

[L, NUM] = bwlabeln(Im>0, 26);
[RegSizes, ~]=hist(L(L>0),1:NUM);
temp=Im(L>0);
temp(RegSizes(L(L>0))<size_thr)=0;
Im(L>0)=temp;

disp('Small region elimination is complete.')
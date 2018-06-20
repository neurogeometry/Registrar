% This function thresholds the stack.

function Im=Threshold(Im,thr_low)

disp('Thresholding started.')

rem_ind=(Im(:)<thr_low);
Im(rem_ind)=0;

disp('Thresholding is complete.')

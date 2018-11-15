% This function thresholds the stack.

function Im=Threshold_Background(Im,thr_low)

disp('Backgroung thresholding started.')

rem_ind=(Im(:)<thr_low);
Im(rem_ind)=0;

disp('Backgroung thresholding is complete.')

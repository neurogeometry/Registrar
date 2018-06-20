% This function calculates correlation coefficient for voxels in the overlap 
% region of Stack1 and Stack2. StackPosition1 and StackPosition2 are the global 
% positions of the stacks. NaN's do not contribute to CC.  

function CC=Stack_Correlation(Stack1,Stack2,StackPosition1,StackPosition2)

StackPosition1=round(StackPosition1);
StackPosition2=round(StackPosition2);

StackSize1=size(Stack1);
if length(StackSize1)==2
    StackSize1(3)=1;
end
StackSize2=size(Stack2);
if length(StackSize2)==2
    StackSize2(3)=1;
end

StackStart1=max(StackPosition2-StackPosition1+1,1);
StackEnd1=min(StackPosition2-StackPosition1+StackSize2,StackSize1);

StackStart2=max(StackPosition1-StackPosition2+1,1);
StackEnd2=min(StackPosition1-StackPosition2+StackSize1,StackSize2);

if nnz(StackStart1<=StackEnd1)==3
    Temp1=Stack1(StackStart1(1):StackEnd1(1),StackStart1(2):StackEnd1(2),StackStart1(3):StackEnd1(3));
    Temp2=Stack2(StackStart2(1):StackEnd2(1),StackStart2(2):StackEnd2(2),StackStart2(3):StackEnd2(3));
    
    ind=(Temp1>0 & Temp2>0);
    
    CC=corrcoef(double(Temp1(ind)),double(Temp2(ind)));
    CC=CC(1,2);
else
    CC=NaN;
end

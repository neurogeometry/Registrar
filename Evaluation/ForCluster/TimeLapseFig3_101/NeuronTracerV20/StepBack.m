% This function determines branch tip vertexes [tip1V, tip2V] of 
% every branch contained in AMlbl by steping back from the tips stepsback micrometers

function [tip1V, tip2V]=StepBack(AMlbl,stepsback)

stepsback=fix(stepsback)+1;
AM=spones(AMlbl);
L=unique(AMlbl(AMlbl>0));

tip1V=cell(length(L),1);
tip2V=cell(length(L),1);
for i=1:length(L)
    [v1 ~]=find(AMlbl==L(i));
    tip12=sort(v1(sum(AM(:,v1))>=3 | sum(AM(:,v1))==1));
    
    if ~isempty(tip12)
        tempAM=(AMlbl==L(i));
        tempV=zeros(1,stepsback+1);
        nextV=tip12(1);
        count=1;
        while count<=stepsback+1 && ~isempty(nextV)
            tempV(count)=nextV;
            nextV=find(tempAM(:,nextV),1,'first'); % in case branch makes a loop
            tempAM(nextV,tempV(count))=0;
            tempAM(tempV(count),nextV)=0;
            count=count+1;
        end
        tempV(tempV==0)=[];
        tip1V{L(i)}=tempV;
        
        tempAM=(AMlbl==L(i));
        tempV=zeros(1,stepsback+1);
        nextV=tip12(2);
        count=1;
        while count<=stepsback+1 && ~isempty(nextV)
            tempV(count)=nextV;
            nextV=find(tempAM(:,nextV),1,'last'); % in case branch makes a loop
            tempAM(nextV,tempV(count))=0;
            tempAM(tempV(count),nextV)=0;
            count=count+1;
        end
        tempV(tempV==0)=[];
        tip2V{L(i)}=tempV;
    else
        tip1V{L(i)}=[];
        tip2V{L(i)}=[];
    end   
end

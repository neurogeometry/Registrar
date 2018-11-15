% This function determines branch tip vertexes, tipV, for 
% every tree contained in AMlbl by steping back from the tips stepsback
% steps

function tipV=StepBackTips(AMlbl,stepsback)

stepsback=ceil(stepsback);
AM=spones(AMlbl);
tips=find(sum(AM)==1);

tipV=cell(length(tips),1);
for i=1:length(tips)
    tempAM=AM;
    tempV=zeros(1,stepsback+1);
    nextV=tips(i);
    count=1;
    while count<=stepsback && length(nextV)==1
        tempV(count)=nextV;
        nextV=find(tempAM(:,nextV));
        tempAM(nextV,tempV(count))=0;
        tempAM(tempV(count),nextV)=0;
        count=count+1;
    end
    tempV(tempV==0)=[];
    tipV{i}=tempV;
end
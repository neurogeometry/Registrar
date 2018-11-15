% Function needed for MatLab - Java data format conversion.
% Training history contained in the ClustersStr is appended to OldTrainingHistory 
% w is updated

function Save_TrainingHistory(filepath,ClustersStr,OldTrainingHistory)

if ~isempty(ClustersStr(1).w)
    TrainingHistory=[ClustersStr,OldTrainingHistory];
else
    TrainingHistory=OldTrainingHistory;
end

save(filepath,'TrainingHistory')
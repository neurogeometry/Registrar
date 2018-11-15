% Function needed for MatLab - Java data format conversion.

function TrainingHistory = Load_TrainingHistory(filepath)

arr = load(filepath);

if isfield(arr,'TrainingHistory')
    TrainingHistory = arr.TrainingHistory;
else
    TrainingHistory = [];
    disp('Incorrect Training History file format.')
end
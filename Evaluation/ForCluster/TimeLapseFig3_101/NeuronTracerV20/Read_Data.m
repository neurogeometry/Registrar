% Function needed for MatLab - Java data format conversion.

function [Original, IM, AM, r, R, reduction_x, reduction_y, reduction_z] = Read_Data(filepath)

arr = load(filepath);

if isfield(arr,'Original') && isfield(arr,'IM') && isfield(arr,'AM') && isfield(arr,'r') && isfield(arr,'R')
    Original = arr.Original;
    IM = arr.IM;
    AM = arr.AM;
    r = arr.r;
    R = arr.R;
    if isfield(arr,'reduction_xy')
        reduction_x= arr.reduction_xy;
        reduction_y= arr.reduction_xy;
        
    elseif isfield(arr,'reduction_x') && isfield(arr,'reduction_y')
        reduction_x= arr.reduction_x;
        reduction_y= arr.reduction_y;
    else
        reduction_x= 1;
        reduction_y= 1;
    end
    
    if isfield(arr,'reduction_z')
        reduction_z= arr.reduction_z;
    else
        reduction_z= 1;
    end
else
    Original = [];
    IM = [];
    AM = [];
    r = [];
    R = [];
    reduction_x= [];
    reduction_y= [];
    reduction_z = [];
    disp('Incorrect Library file format.')
end
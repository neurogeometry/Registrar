function UpdateLog(LogHandle,String)
try
    if isempty(LogHandle)
        Log();
        LogHandle=findobj(0,'Name','Log');
        LogHandle.Children(2).String = {};
    end
    LogHandle.Children(2).String{end+1} = String;
    LogHandle.Children(2).Value = size(LogHandle.Children(2).String,1);
catch
end


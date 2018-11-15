function Parameters=Read_Parameter_File(file_path)

Parameters=[];
fid = fopen(file_path);

tline = fgetl(fid);
while ischar(tline)
    tline = strtrim(tline);
    if ~isempty(tline) && ~strcmp(tline(1),'%')
        f1=find(tline=='.');
        f2=find(tline=='=');
        Str1=strtrim(tline(1:f1-1));
        Str2=strtrim(tline(f1+1:f2-1));
        Val=str2num(strtrim(tline(f2+1:end)));
        
        if ~isfield(Parameters,Str1)
            temp=[];
            temp=setfield(temp,Str2,Val);
            Parameters=setfield(Parameters,Str1,temp);
        else
            temp=getfield(Parameters,Str1);
            temp=setfield(temp,Str2,Val);
            Parameters=setfield(Parameters,Str1,temp);
        end
    end
    tline = fgetl(fid);
end

fclose(fid);
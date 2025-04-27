function fileContents = logger(path,log,opt,n)
switch opt  
    case 'read' 
        fileID = fopen(path,'r');
        ret = textscan(fileID, '%s', [1, 4]);  % 使用换行符作为分隔符 
        fileContents = ret{1}{1};
        fclose(fileID); 
    case 'write' 
        fileID = fopen(path, 'a');  
        fprintf(fileID, '%s\n', log);  
        fclose(fileID);
        fileContents = log;
    case 'replace' 
        fileID = fopen(path, 'r');
        ret = textscan(fileID, '%s', [1, 4]);  % 使用换行符作为分隔符 
        fileContents = ret{1}{1};
        fileContents(n)=log;
        fclose(fileID); 
        fileID = fopen(path, 'w');  
        fprintf(fileID, '%s', fileContents);  
        fclose(fileID);
        fileContents = log;
    otherwise  
        disp('The value is something else');  
end
end
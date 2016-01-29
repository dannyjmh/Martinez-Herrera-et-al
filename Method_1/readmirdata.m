function [num txt] = readmirdata(FileName, columns)
    fileID = fopen(FileName);
    C = textscan(fileID, '%s', columns,'Delimiter','\t');
    txt = C{1};
    fclose(fileID);
    num = dlmread(FileName, '\t',1,1);
end 
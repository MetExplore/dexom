function saveArrayStrings(file, cellArray)
    filePh = fopen(file,'w');
    fprintf(filePh,'%s\n',cellArray{:});
    fclose(filePh);
end


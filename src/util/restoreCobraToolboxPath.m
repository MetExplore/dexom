function restoreCobraToolboxPath()
    tempFile = [fileparts(which('dexomInit')) filesep 'cobrapath.old'];
    if isempty(tempFile)
        error('Missing cobrapath.old file');
    end
    fid = fopen(tempFile);
    while ~feof(fid)
        fpath = fgetl(fid);
        addpath(fpath);
    end
    fclose(fid);
    removeCobraFromPath('initializeCobraToolboxLib');
end


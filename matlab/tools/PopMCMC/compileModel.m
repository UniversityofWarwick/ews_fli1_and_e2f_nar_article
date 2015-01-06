function [name] = compileModel(modelVf, outPath)

name = [];

if ~exist(modelVf, 'file')
    fprintf('%s\n', modelVf);
    error('compileModel: FileNotFound', 'vf file not found');
end

[inPath modelname ext] = fileparts(modelVf);

if exist(fullfile(inPath, [modelname '_cvs.c.so']),'file')
    name = fullfile(inPath, [modelname '_cvs.c.so']);
    return 
end

user_pwd = pwd;

if ~isempty(inPath)
    cd(inPath)
else
    inPath = user_pwd;
end

[status, result] = system(['vfgen cvodes:sens=yes ' modelVf]);

if status == 0
    cfile = fullfile(inPath, [modelname '_cvs.c']);
    if ~exist(cfile, 'file')
        fprintf('%s\n', cfile);
        error('compileModel: FileNotFound', 'C file not found. Make sure the vf file name matches the model name in the coresponding xml element.');
    end

    [status, result] = system(['ginac-excompiler ' cfile]);
    
    if status == 0    
       name = fullfile(inPath, [modelname '_cvs.c.so']);
       if ~exist(name, 'file')
           fprintf('%s\n', name);
           error('compileModel: FileNotFound', 'Compiled file not found. Check file names.');
       end
    else
        fprintf('%s\n%s', status, result);
        error('compileModel: ginacError', 'Error compiling the C file.');
    end
else
    fprintf('%s\n%s', status,result);
    error('compileModel: vfgenError', 'Error executing vfgen.');
end

if nargin == 2
  movefile(name, outPath);
  movefile(cfile, outPath);
  movefile(fullfile(inPath, [modelname '_cvs.h']), outPath);
end
  
cd(user_pwd);

end

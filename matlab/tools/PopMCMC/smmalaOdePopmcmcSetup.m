function [xPosterior, BurnInPath, LL] = ...
  smmalaOdePopmcmcSetup(numMcmcIters, numMcmcBurnin, ...
  thinBurnin, verboseStep, ...
  logFile, modelVfFile, dataMatFile, initX, SOmodelFile,remoteOutFile)


% SET PATHS
%addpath(genpath('/home/ucakvst/matlab_tools/'));

%sys_CINCL = getenv('C_INCLUDE_PATH');
%setenv('C_INCLUDE_PATH',[sys_CINCL  ':/home/ucakvst/local/include']);

%sys_LDLIB = getenv('LD_LIBRARY_PATH');
%setenv('LD_LIBRARY_PATH',['/home/ucakvst/local/lib']);

%sys_PATH = getenv('PATH');
%setenv('PATH',[sys_PATH ':/home/ucakvst/local/bin']);

addpath(genpath('/home/ucaktpa/opt/matlab/tools/odeSolverMex'));

sys_CINCL = getenv('C_INCLUDE_PATH');
setenv('C_INCLUDE_PATH',[sys_CINCL  ':/home/ucaktpa/opt/sundials/2.5.0/include']);

sys_LDLIB = getenv('LD_LIBRARY_PATH');
setenv('LD_LIBRARY_PATH',[sys_LDLIB  ':/home/ucaktpa/opt/sundials/2.5.0/lib']);

sys_PATH = getenv('PATH');
setenv('PATH',[sys_PATH ':/home/ucaktpa/opt/vfgen/2.4.1/bin/:/home/ucaktpa/opt/ginac/1.6.2/bin']);


randStream = RandStream('mt19937ar','Seed',labindex);
RandStream.setGlobalStream(randStream);



%fprintf(tmpLogFile01, 'Random stream set\n')

data = load(dataMatFile);
%fprintf(tmpLogFile01, 'Data loaded\n')

%dlmwrite('/home/ucaktpa/tmp/tmpDataFile01.txt', data.Observations);
%fprintf(tmpLogFile01, 'Wrote Observations to file\n')
%dlmwrite('/home/ucaktpa/tmp/tmpDataFile02.txt', data.TimePoints);
%fprintf(tmpLogFile01, 'Wrote TimePoints to file\n')

%if exist(modelVfFile, 'file')
%  fprintf(tmpLogFile01, 'Model VF file exists\n')
%end
%if exist(['/home/ucaktpa/opt/lib/' modelSoLibFile], 'file')
%  fprintf(tmpLogFile01, 'Model dynamic library file exists\n')
%end
odeModel = odeModelDist.odeModelDistFromData(data, modelVfFile, SOmodelFile);
%fprintf(tmpLogFile01, 'odeModel instantiated\n')

%dlmwrite('/home/ucaktpa/tmp/tmpDataFile03.txt', odeModel.Data.Observations)
%fprintf(tmpLogFile01, 'Wrote Observations from odeModel to file\n')
%dlmwrite('/home/ucaktpa/tmp/tmpDataFile04.txt', odeModel.Data.TimePoints)
%fprintf(tmpLogFile01, 'Wrote TimePoints from odeModel to file\n')

%fclose(tmpLogFile01)


[xPosterior, BurnInPath, LL] = ...
  smmalaOdePopmcmc(odeModel, ...
  'Iterations', numMcmcIters, ...
  'BurnIn', numMcmcBurnin, ...
  'Verbose', verboseStep, ...
  'LogFile', logFile, ...
  'ThinBurnIn', thinBurnin,...
  'Init', initX(:,labindex) );
remoteOutFile = [remoteOutFile '_' num2str(labindex) '.mat'];
save(remoteOutFile,'xPosterior','BurnInPath','LL');
end

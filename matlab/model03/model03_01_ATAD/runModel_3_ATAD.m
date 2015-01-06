%function [pjob, outFile] = runModel_3_ATAD()
%=== SET EXPERIMENT PARAMETERS 
logFile = '/home/ucaktpa/model3_ATAD/log/MPI_POP_SMMALA_job_Model_3_ATAD.txt';
modelFile = 'ewsfli1_Model_3_ATAD.vf';
somodelFile = '/home/ucaktpa/opt/lib/ewsfli1_kd_model03_01_02_cvs.c.so';
%somodelFile = '/home/ucakvst/e2f3_so/ewsfli1_kd_model01_01_02_cvs.c.so';
remoteOutFile = '/home/ucaktpa/model3_ATAD/posteriorModel_3_ATAD';
dataFile = 'dataATAD.mat';
out_m3_atad = 'posteriorModel_3_ATAD.mat';


%=== GET INITIAL CONDITIONS FROM PREVIOUS SIMULATION
pPostx = load('prevSimModel3.mat');
% take last row from each matrix
initX = cellfun(@(x) x(end,:), pPostx.results, 'UniformOutput', false);
% first cell column is the posterior
initX = cell2mat(initX(:,1));
% use the common variance init for both rna and protein variances
initX = [initX initX(:,end)]';


sched = findResource('scheduler', 'configuration', 'genericconfig1');

job_m3_atad = createParallelJob(sched);

set(job_m3_atad, 'MaximumNumberOfWorkers', 50);
set(job_m3_atad, 'MinimumNumberOfWorkers', 50);

Dependencies = {'odeLogLikelihood.m', ...
  'smmalaOdePopmcmcSetup.m', ...
  'smmalaOdePopmcmc.m', ...
  'compileModel.m',...
  'odeModelDist.m', ...
  modelFile, ...
  dataFile };

set(job_m3_atad, 'FileDependencies', Dependencies);

%=== MCMC PARAMETERS
Iterations = 25000;
Burnin = 5000;
ThinBurnin = 10;
Adapt = 100;

%=== MCMC PARAMETERS
%Iterations = 10;
%Burnin = 5;
%ThinBurnin = 1;
%Adapt = 2;

t_m3_atad = createTask(job_m3_atad, @() smmalaOdePopmcmcSetup(Iterations, Burnin, ...
  ThinBurnin, Adapt, ...
  logFile, modelFile, dataFile, initX, somodelFile,remoteOutFile), 3, {});

submit(job_m3_atad);

waitForState(job_m3_atad);

results = getAllOutputArguments(job_m3_atad);

save(out_m3_atad, 'results');

%destroy(pjob);

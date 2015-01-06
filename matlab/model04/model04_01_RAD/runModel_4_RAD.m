%function [job_m4_rad, outFile] = runModel_4_RAD()
%=== SET EXPERIMENT PARAMETERS 
logFile = '/home/ucaktpa/model4_RAD/log/MPI_POP_SMMALA_job_Model_4_RAD.txt';
modelFile = 'ewsfli1_Model_4_RAD.vf';
somodelFile = '/home/ucaktpa/opt/lib/ewsfli1_kd_model04_01_02_cvs.c.so';
%somodelFile = '/home/ucakvst/e2f3_so/ewsfli1_kd_model01_01_02_cvs.c.so';
remoteOutFile = '/home/ucaktpa/model4_RAD/posteriorModel_4_RAD';
dataFile = 'dataRAD.mat';
out_m4_rad = 'posteriorModel_4_RAD.mat';


%=== GET INITIAL CONDITIONS FROM PREVIOUS SIMULATION
pPostx = load('prevSimModel4.mat');
% take last row from each matrix
initX = cellfun(@(x) x(end,:), pPostx.results, 'UniformOutput', false);
% first cell column is the posterior
initX = cell2mat(initX(:,1));
% use the common variance init for both rna and protein variances
initX = [initX initX(:,end)]';


sched = findResource('scheduler', 'configuration', 'genericconfig1');

job_m4_rad = createParallelJob(sched);

set(job_m4_rad, 'MaximumNumberOfWorkers', 50);
set(job_m4_rad, 'MinimumNumberOfWorkers', 50);

Dependencies = {'odeLogLikelihood.m', ...
  'smmalaOdePopmcmcSetup.m', ...
  'smmalaOdePopmcmc.m', ...
  'compileModel.m',...
  'odeModelDist.m', ...
  modelFile, ...
  dataFile };

set(job_m4_rad, 'FileDependencies', Dependencies);

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

t_m4_rad = createTask(job_m4_rad, @() smmalaOdePopmcmcSetup(Iterations, Burnin, ...
  ThinBurnin, Adapt, ...
  logFile, modelFile, dataFile, initX, somodelFile,remoteOutFile), 3, {});

submit(job_m4_rad);

waitForState(job_m4_rad);

results = getAllOutputArguments(job_m4_rad);

save(out_m4_rad, 'results');


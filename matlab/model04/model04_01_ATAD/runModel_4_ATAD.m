%function [pjob, outFile] = runModel_4_ATAD()
%=== SET EXPERIMENT PARAMETERS 
logFile = '/home/ucaktpa/model4_ATAD/log/MPI_POP_SMMALA_job_Model_4_ATAD.txt';
modelFile = 'ewsfli1_Model_4_ATAD.vf';
somodelFile = '/home/ucaktpa/opt/lib/ewsfli1_kd_model04_01_02_cvs.c.so';
%somodelFile = '/home/ucakvst/e2f3_so/ewsfli1_kd_model01_01_02_cvs.c.so';
remoteOutFile = '/home/ucaktpa/model4_ATAD/posteriorModel_4_ATAD';
dataFile = 'dataATAD.mat';
out_m4_atad = 'posteriorModel_4_ATAD.mat';


%=== GET INITIAL CONDITIONS FROM PREVIOUS SIMULATION
pPostx = load('prevSimModel4.mat');
% take last row from each matrix
initX = cellfun(@(x) x(end,:), pPostx.results, 'UniformOutput', false);
% first cell column is the posterior
initX = cell2mat(initX(:,1));
% use the common variance init for both rna and protein variances
initX = [initX initX(:,end)]';


sched = findResource('scheduler', 'configuration', 'genericconfig1');

job_m4_atad = createParallelJob(sched);

set(job_m4_atad, 'MaximumNumberOfWorkers', 50);
set(job_m4_atad, 'MinimumNumberOfWorkers', 50);

Dependencies = {'odeLogLikelihood.m', ...
  'smmalaOdePopmcmcSetup.m', ...
  'smmalaOdePopmcmc.m', ...
  'compileModel.m',...
  'odeModelDist.m', ...
  modelFile, ...
  dataFile };

set(job_m4_atad, 'FileDependencies', Dependencies);

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

t_m4_atad = createTask(job_m4_atad, @() smmalaOdePopmcmcSetup(Iterations, Burnin, ...
  ThinBurnin, Adapt, ...
  logFile, modelFile, dataFile, initX, somodelFile,remoteOutFile), 3, {});

submit(job_m4_atad);

waitForState(job_m4_atad);

results = getAllOutputArguments(job_m4_atad);

save(out_m4_atad, 'results');

%destroy(pjob);

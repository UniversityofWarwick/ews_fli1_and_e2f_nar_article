%function [pjob, outFile] = runModel_1_RAD()
%=== SET EXPERIMENT PARAMETERS 
logFile = '/home/ucaktpa/model1_RAD/log/MPI_POP_SMMALA_job_Model_1_RAD.txt';
modelFile = 'ewsfli1_Model_1_RAD.vf';
somodelFile = '/home/ucaktpa/opt/lib/ewsfli1_kd_model01_01_02_cvs.c.so';
%somodelFile = '/home/ucakvst/e2f3_so/ewsfli1_kd_model01_01_02_cvs.c.so';
remoteOutFile = '/home/ucaktpa/model1_RAD/posteriorModel_1_RAD';
dataFile = 'dataRAD.mat';
out_m1_rad = 'posteriorModel_1_RAD.mat';


%=== GET INITIAL CONDITIONS FROM PREVIOUS SIMULATION
pPostx = load('prevSimModel1.mat');
% take last row from each matrix
initX = cellfun(@(x) x(end,:), pPostx.results, 'UniformOutput', false);
% first cell column is the posterior
initX = cell2mat(initX(:,1));
% use the common variance init for both rna and protein variances
initX = [initX initX(:,end)]';


sched = findResource('scheduler', 'configuration', 'genericconfig1');

job_m1_rad = createParallelJob(sched);

set(job_m1_rad, 'MaximumNumberOfWorkers', 50);
set(job_m1_rad, 'MinimumNumberOfWorkers', 50);

Dependencies = {'odeLogLikelihood.m', ...
  'smmalaOdePopmcmcSetup.m', ...
  'smmalaOdePopmcmc.m', ...
  'compileModel.m',...
  'odeModelDist.m', ...
  modelFile, ...
  dataFile };

set(job_m1_rad, 'FileDependencies', Dependencies);

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

t_m1_rad = createTask(job_m1_rad, @() smmalaOdePopmcmcSetup(Iterations, Burnin, ...
  ThinBurnin, Adapt, ...
  logFile, modelFile, dataFile, initX, somodelFile,remoteOutFile), 3, {});

submit(job_m1_rad);

waitForState(job_m1_rad);

results = getAllOutputArguments(job_m1_rad);

save(out_m1_rad, 'results');

%destroy(pjob);

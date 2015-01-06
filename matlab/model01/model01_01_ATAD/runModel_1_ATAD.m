clear;
%=== SET EXPERIMENT PARAMETERS 
logFile = '/home/ucaktpa/model1_ATAD/log/MPI_POP_SMMALA_job_Model_1_ATAD.txt';
modelFile = 'ewsfli1_Model_1_ATAD.vf';
somodelFile = '/home/ucaktpa/opt/lib/ewsfli1_kd_model01_01_02_cvs.c.so';
%somodelFile = '/home/ucakvst/e2f3_so/ewsfli1_kd_model01_01_02_cvs.c.so';
remoteOutFile = '/home/ucaktpa/model1_ATAD/posteriorModel_1_ATAD';
dataFile = 'dataATAD.mat';
outFile = 'posteriorModel_1_ATAD.mat';


%=== GET INITIAL CONDITIONS FROM PREVIOUS SIMULATION
pPostx = load('prevSimModel1.mat');
% take last row from each matrix
initX = cellfun(@(x) x(end,:), pPostx.results, 'UniformOutput', false);
% first cell column is the posterior
initX = cell2mat(initX(:,1));
% use the common variance init for both rna and protein variances
initX = [initX initX(:,end)]';


sched = findResource('scheduler', 'configuration', 'genericconfig1');

pjob = createParallelJob(sched);

set(pjob, 'MaximumNumberOfWorkers', 50);
set(pjob, 'MinimumNumberOfWorkers', 50);

Dependencies = {'odeLogLikelihood.m', ...
  'smmalaOdePopmcmcSetup.m', ...
  'smmalaOdePopmcmc.m', ...
  'compileModel.m',...
  'odeModelDist.m', ...
  modelFile, ...
  dataFile };

set(pjob, 'FileDependencies', Dependencies);

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

t = createTask(pjob, @() smmalaOdePopmcmcSetup(Iterations, Burnin, ...
  ThinBurnin, Adapt, ...
  logFile, modelFile, dataFile, initX, somodelFile,remoteOutFile), 3, {});

submit(pjob);

waitForState(pjob);

results = getAllOutputArguments(pjob);

save(outFile, 'results');

%destroy(pjob);

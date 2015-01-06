clear;
%ATAD data
TimePoints = [0 3 6 10 12:3:42 44 48];

Observations = csvread('ewsfli1_kd_2013_08_23_experiment04_hrs00_48.csv');
Observations(Observations==-1) = nan;
Observations = Observations(:, [1:4, 6])';

TransformationMatrix = [];
StateIdx{1} = [1 3 5];  %RNA states
StateIdx{2} = [ 2 4 ];  %Protein states

outfile = './dataATAD.mat';
save(outfile, 'TimePoints', 'Observations', 'TransformationMatrix', 'StateIdx');


clear;
%RAD data

TimePoints = [0 3 6 10 12:3:42 44 48];

Observations = csvread('ewsfli1_kd_2013_08_23_experiment04_hrs00_48.csv');
Observations(Observations==-1) = nan;
Observations = Observations(:, 1:5)';

TransformationMatrix = [];
StateIdx{1} = [1 3 5];  %RNA states
StateIdx{2} = [ 2 4 ];  %Protein states

outfile = './dataRAD.mat';
save(outfile, 'TimePoints', 'Observations', 'TransformationMatrix','StateIdx');
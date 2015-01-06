function [logModelMarginal] = computeMarginal(resultsFile)

post = load(resultsFile);


a = 0:0.02:1;
t = a.^4;

LLCmpValues = cell2mat(post.results(:,3)')';
[nPowers, nSamples] = size(LLCmpValues);

LLCmpMeans = mean(LLCmpValues, 2);

logModelMarginal = 0;
for i = 1:(nPowers-1)
  logModelMarginal = logModelMarginal + (t(i+2)-t(i+1))*(LLCmpMeans(i+1)+LLCmpMeans(i));
end
logModelMarginal = logModelMarginal/2;



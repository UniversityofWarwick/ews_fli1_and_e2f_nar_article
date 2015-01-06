function posModFig = responsePosteriorPlot(odeSystem, dataFile, xPostFile, varargin)
%function posModFig = responsePosteriorPlot(odeSystem, Ny, Np, xPost, timePoints, varargin)

%========== Parsing Parameters ===========
argParser = inputParser;
argParser.addParamValue('Samples', 500, @(x) isscalar(x) && x >= 0 );
argParser.addParamValue('Remove', 1, @(x) isscalar(x) && x >= 0 );
argParser.parse( varargin{:} );
Samples =  argParser.Results.Samples;
Remove = argParser.Results.Remove;
clear argParser;

%======================================%
%%% PLOT Posterior responses         %%%
%======================================%
xpost = load(xPostFile);
x = xpost.results{end,1};
x = x(Remove:end,:);

odeSys = compileModel(which(odeSystem));
data = load(dataFile);
odeModel = odeModelDist.odeModelDistFromData(data,odeSystem,odeSys);


idx = 1:floor(size(x,1)/Samples):size(x,1);

Ny = length(odeModel.States);

plotR = 2;
plorC = 3;
%plotR = ceil( sqrt(Ny) );
%plorC = plotR;

posModFig = figure;
for i = 1:Ny
    subplot(plotR,plorC,i)
    ylabel(odeModel.States(i).Name,'FontSize',14)
    xlabel('time','FontSize',14)
    set(gca,'FontSize',14)
    hold on
end

postMean = zeros(Ny, length(data.TimePoints),Samples);

for i = 1:Samples
    params = x(idx(i),:);
    Yp = odeModel.simulateODE(data.TimePoints, params');
    
    %epsilon = sqrt(10.^params(end-1:end));
    %Yp(data.StateIdx{1},:) =  Yp(data.StateIdx{1},:)+normrnd(0,epsilon(1),length(data.StateIdx{1}),length(data.TimePoints));
    %Yp(data.StateIdx{2},:) =  Yp(data.StateIdx{2},:)+normrnd(0,epsilon(2),length(data.StateIdx{2}),length(data.TimePoints));
    postMean(:,:,i) = Yp;
    
    for j = 1:Ny
        subplot(plotR,plorC, j)
        patchline(data.TimePoints,Yp(j,:),'linestyle','-','edgecolor','b','linewidth',1,'edgealpha',0.05);
    end
end
postMean = median(postMean,3);

for i = 1:Ny
    subplot(plotR,plorC,i)
    line(data.TimePoints,postMean(i,:),ones(1,length(data.TimePoints)), 'Color','green');
    scatter3(data.TimePoints, data.Observations(i,:), ones(1,length(data.TimePoints)),'MarkerEdgeColor','red');
%     
%     if ~isempty(Data)
%         if ~isempty(TransormationMatrix)
%             dataIdx = find(TransormationMatrix(:,i));
%         else
%             dataIdx = i;
%         end
%         if ~isempty(dataIdx)
%            scatter3(timePoints, Data(dataIdx,:), ones(1,length(timePoints)),'MarkerEdgeColor','red');
%            %line(timePoints, Data(dataIdx,:), ones(1,length(timePoints)), 'Color','red','LineWidth',1.0);
%         else
%             ylim( [min(postMean(i,:)) 2*max(postMean(i,:))] );
%         end
%     end
%     if i == 5
%         ylim( [min(postMean(i,:)) 2*max(postMean(i,:))] );
%     end
    %ylim( [min(postMean(i,:)) 2*max(postMean(i,:))] );
    xlim( [0 max(data.TimePoints)]);
    hold off
end

end

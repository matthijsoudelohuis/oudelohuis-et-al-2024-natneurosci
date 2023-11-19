function temptrialData = MOL_StratLicks(temptrialData)

params.createprobetrials    = 1;
params.binres               = 50e3;
params.minresplat           = 100e3;
params.maxresplat           = 700e3;    %
params.tolerance            = 10;       %

trialtypes          = {'Y' 'X' 'C' 'P'};
colorstrialType     = {[1 0.2 0.2] [0.2 0.2 1] [0.9 0 0.9] [0.5 0.5 0.5]};

edges = params.minresplat:params.binres:params.maxresplat;
resplat_dist = NaN(length(trialtypes),length(edges)-1);
for iTr = 1:length(trialtypes)
     resplat_dist(iTr,:) = histcounts(temptrialData.responseLatency(strcmp(temptrialData.trialType,trialtypes{iTr})),edges);
     resplat_dist(iTr,:) = histcounts(temptrialData.responseLatency(strcmp(temptrialData.trialType,trialtypes{iTr}) & temptrialData.correctResponse==1),edges);

end

figure; hold all;
for iTr = 1:length(trialtypes)
    plot(edges(1:end-1),resplat_dist(iTr,:),'color',colorstrialType{iTr},'LineWidth',2)
end





end

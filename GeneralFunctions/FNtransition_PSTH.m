function [PSTH,FNratio,ts] = FNtransition_PSTH(metaDS,signal)
% Plot calcium traces aligned to the familiar to novel transition for
% individual cells or boutons.
% 
% Thomas Hainmueller, 2023

if nargin<3
    signal = 'dFoT';
end

fr = 16; % Hz, standard for axon imaging datasets
interval = 60; % In seconds, symmetrical around zero
ts = -interval:1/fr:interval-1/fr;
colors = {[1 .5 .5],[1 .7 .7]};
minmaxnorm = true;

for ds = 1:length(metaDS)
    
    data = metaDS{ds};
    PSTH{ds} = NaN(2*interval*fr,length(data.cells));
    
    %fr = data.metadata.categories{1}.acquisition_rate{1};
    nruns = min(length(data.metadata.categories{1}.y),...
        length(data.metadata.categories{2}.y));
    
    for n = 1:length(data.cells)
        for r = 1:nruns
            Ftr = data.cells{n}.categories{1}.(signal){r};
            Ntr = data.cells{n}.categories{2}.(signal){r};
            
            thisPSTH(1:interval*fr,r) = Ftr(end-interval*fr+1:end);
            thisPSTH(interval*fr+1:2*interval*fr,r) = Ntr(1:interval*fr);
            %trace = cat(2,Ftr,Ntr);
            %ts = (-length(Ftr):length(Ntr)-1)/fr;
        end
        PSTH{ds}(:,n) = nanmean(thisPSTH,2);
        
        if minmaxnorm
            PSTH{ds}(:,n) = PSTH{ds}(:,n)-min(PSTH{ds}(:,n));
            PSTH{ds}(:,n) = PSTH{ds}(:,n)/max(PSTH{ds}(:,n));
        end
        
        FNratio{ds}(n) = nanmean(PSTH{ds}(1:interval*fr,n)) / nanmean(PSTH{ds}(interval*fr+1:end,n));
    end
end

PSTH = cat(2,PSTH{:});
FNratio = cat(2,FNratio{:});

end
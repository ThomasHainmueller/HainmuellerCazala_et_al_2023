function data = transientrates(data,moving)
% Define the rate of Calcium signals for all cells in all categories of a
% dataset. Store the average rate for each category. Also stores the
% the acquisition rate with each trace in the metadata. 
% TODO=Swap this step to the dataset conversion routine.
if nargin < 2
    % Use only periods where the animal was moving.
    moving = true;
end

for n = 1:length(data.cells)
    for c = 1:length(data.cells{n}.categories)
        transients = deal(data.cells{n}.categories{c}.transientmask);
        transients = cat(2,transients{:}); % 170526 adapted for recent ds
        
        signals = deal(data.cells{n}.categories{c}.dFoT);
        signals = cat(2,signals{:});
        
        if moving
            lenmoving = deal(data.metadata.categories{c}.moving);
            lenmoving = cat(2,lenmoving{:});
            transients = transients.*lenmoving;
            lenmoving = length(find(lenmoving>0));
        else
            lenmoving = length(transients);
        end
        % Mask the calciumsignals to the periods of transients and moving
        signals =signals.*transients;
        signalAUC = sum(signals);
        
        ntransients = length(find(diff(transients)>0));

        %framerate = deal(data.metadata.categories{c}.acquisition_rate);
        %framerate = mean(cat(2,framerate{:})); % frames/second
        framerate = mean(data.metadata.categories{c}.acquisition_rate);
        
        % Calcium transients per second
        data.cells{n}.transientrate(c)=ntransients/lenmoving*framerate;
        % AUC of Calcium signal per second.
        data.cells{n}.AUCrate(c)=signalAUC/lenmoving*framerate;
    end
end
end
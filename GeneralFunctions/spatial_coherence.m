function data = spatial_coherence(data)
% spatial_coherence(data). Calculate the spatial coherence based on the
% mean dFoY data stored with the dataset as the mean cross-correlation of
% each pixel with its two nearest neighbours, excluding first and last
% element.

for c = 1:length(data.metadata.categories)
    for n = 1:length(data.cells)
        % Get average Fluorescence over track distance (Y)
        thisdFoY = deal(data.cells{n}.categories{c}.dFoY);
        thisdFoY = cat(2,thisdFoY{:});
        thisdFoY = nanmean(thisdFoY,2);
        
        % Calculate this for z-scored traces, too.
%         thisdFoYz = (thisdFoY-nanmean(thisdFoY))./nanstd(thisdFoY);
        
        %figure; hold on; 
        % Calculate the trace for neighbours
        for f = length(thisdFoY)-1:-1:2
            thisdFoYsm(f-1) = nanmean([thisdFoY(f+1),thisdFoY(f-1)]);
%             thisdFoYsmz(f-1) = nanmean([thisdFoYz(f+1),thisdFoYz(f-1)]);
        end
        
        %plot(thisdFoYsm); plot(thisdFoY(2:end-1));
        
        % Calculate correlation
        thiscorr = corrcoef(thisdFoY(2:end-1),thisdFoYsm);
        data.cells{n}.spatial_coherence(c)= thiscorr(1,2);
        
%         thiscorrz = corrcoef(thisdFoYz(2:end-1),thisdFoYsmz);
%         data.cells{n}.spatial_coherence_z(c)= thiscorrz(1,2);
    end
end
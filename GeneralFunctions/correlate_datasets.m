function PVdata = correlate_datasets(PVdata,PCdata,window,filter)
% Take a dataset with PVcells and append a matrix with the correlations
% with all cells in PCdata, together with the inter-somatic distances.
% Default for significant traces of PCs and the unmasked PV cell signal
% restricted to movement periods.
if nargin<4
    filter=fspecial('gaussian',[9,1],1); % input =1 turns filtering of;
end
if nargin<3
    window = 1; % Window for xcorr, default off, i.e. 200 ms bins for 5Hz (=sampling rate).
end

for PVno=1:length(PVdata.cells)
    % PART1: Correlate PV-cell and PCs
    % Get distances between this PVcell and all PC's
    for PCno=length(PCdata.cells):-1:1
        thisdist=pdist([PVdata.cells{PVno}.roi.center;...
            PCdata.cells{PCno}.roi.center],'euclidean');
        for c=1:length(PVdata.cells{PVno}.categories)
            thisPVtr=[];
            thisPCtr=[];
            for run=length(PVdata.cells{PVno}.categories{c}.dFoT):-1:1
                thismov=logical(PVdata.metadata.categories{c}.moving{run});
                thisPVtr=cat(2,PVdata.cells{PVno}.categories{c}.dFoT{run}(thismov),...
                    thisPVtr);
                thisPCtr=cat(2,PCdata.cells{PCno}.categories{c}.dFoT{run}(thismov).*...
                    PCdata.cells{PCno}.categories{c}.transientmask{run}(thismov),thisPCtr);
            end
            % Apply bins for the signals
            thisPVtr=resample(thisPVtr,1,window);
            thisPCtr=resample(thisPCtr,1,window);
            % Apply filter to data
            thisPVtr=conv(thisPVtr,filter,'same');
            thisPCtr=conv(thisPCtr,filter,'same');
            % Subtract mean from signals to get zero centered correlations
            thisPVtr=thisPVtr-mean(thisPVtr);
            thisPCtr=thisPCtr-mean(thisPCtr);
            %PVdata.cells{PVno}.categories{c}.PCcorr(PCno,1:2)=...
                %[corr(thisPVtr(:),thisPCtr(:)),thisdist];
            PVdata.cells{PVno}.categories{c}.PCcorr(PCno,1:2)=...
               [xcorr(thisPVtr(:),thisPCtr(:),0,'coef'),thisdist];
        end
    end
    
    % PART2 correlate PV-cell (PV) with all other PV cells(CPV)
    for CPVno=length(PVdata.cells):-1:1
        thisdist=pdist([PVdata.cells{PVno}.roi.center;...
            PVdata.cells{CPVno}.roi.center],'euclidean');
        for c=1:length(PVdata.cells{PVno}.categories)
            thisPVtr=[];
            thisCPVtr=[];
            for run=length(PVdata.cells{PVno}.categories{c}.dFoT):-1:1
                thismov=logical(PVdata.metadata.categories{c}.moving{run});
                thisPVtr=cat(2,PVdata.cells{PVno}.categories{c}.dFoT{run}(thismov),...
                    thisPVtr);
                thisCPVtr=cat(2,PVdata.cells{CPVno}.categories{c}.dFoT{run}(thismov),...
                    thisCPVtr);
            end
            thisPVtr=resample(thisPVtr,1,window);
            thisCPVtr=resample(thisCPVtr,1,window);
            % Apply filter to data
            thisPVtr=conv(thisPVtr,filter,'same');
            thisCPVtr=conv(thisCPVtr,filter,'same');
            % Mean subtraction for zero-centered correlations
            thisPVtr=thisPVtr-mean(thisPVtr);
            thisCPVtr=thisCPVtr-mean(thisCPVtr);
            PVdata.cells{PVno}.categories{c}.PVcorr(CPVno,1:2)=...
                [xcorr(thisPVtr(:),thisCPVtr(:),0,'coef'),thisdist]; 
        end    
    end 
    warning([num2str(PVno),' cells processed.']);
end
end
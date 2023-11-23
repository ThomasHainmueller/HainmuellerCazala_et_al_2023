function data = placefield_COM(data, sigma)
% Calculate the center of mass for each placefield of a cell. If a cell has
% two placefields on one track, the COM would be located in the middle
% between them.
% simga: filter dFoY with gaussian of width 'sigma' bins before calculation

if nargin<2
    sigma=1; % Gaussian function with sigma = 1 bin
end

if sigma
    filter = fspecial('gaussian',[6*sigma,1],sigma);
end

for n=length(data.cells):-1:1
    for c=length(data.cells{n}.categories):-1:1
        if data.cells{n}.Placefield_P(c) < .05;
            % Get the placefield (logical) as predetermined
            thisfield = data.cells{n}.categories{c}.placefield;
            
            % Get mean dFoY trace from a random cell in this category
            thisdFoY=deal(data.cells{n}.categories{c}.dFoY);
            thisdFoY=cat(2,thisdFoY{:});
            thisdFoY=nanmean(thisdFoY,2);

            if sigma
                nanvals = isnan(thisdFoY);
                thisdFoY(nanvals)=0;
                thisdFoY=conv(thisdFoY,filter,'same');
                thisdFoY(nanvals)=NaN;
            end
            
            if length(thisfield)~=length(thisdFoY)
                thisfield = thisfield(1:2:end); % Usually, PF are determined on 80 bin DS
            end
            
            % Center of mass of avg dF/F in the placefield
            data.cells{n}.Placefield_COM(c) =...
                sum(find(thisfield)'.*(thisdFoY(thisfield)/sum(thisdFoY(thisfield))));
        else
            data.cells{n}.Placefield_COM(c) = NaN;
        end
    end
end
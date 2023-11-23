function data = recalculate_dFoY(data,moving,transients,bins,zscored)
% Repeat the calculation of dFoY maps for an input dataset.
if nargin < 5
    zscored = false;
end
if nargin < 4
    bins = 0.1:0.025:2.1;
end
if nargin < 3
    % Use only significant transients.
    transients = true;
end
if nargin < 2
    % Use only moving periods
    moving = true;
end

for n = 1:length(data.cells)
    for c = 1:length(data.cells{n}.categories)
        for tr = 1:length(data.cells{n}.categories{c}.dFoT)
            % To even out inconsitencies in vector orientation
            vshape = size(data.cells{n}.categories{c}.dFoT{tr});
            thisy = data.metadata.categories{c}.y{tr};
            
            % CAVE: Check if this works with the current y-traces - a
            % factor of 2 may be missing
            
            if iscell(data.metadata.categories{c}.acquisition_rate)
                fr = data.metadata.categories{c}.acquisition_rate{tr};
            else
                fr = data.metadata.categories{c}.acquisition_rate(tr);
            end
            
            speedbins = (0:.01:.3)/fr; % 2V = 4m, 1cm/s increments max 30cm/s
            
            if zscored
                thistrace = data.cells{n}.categories{c}.zscored{tr};
            else
                thistrace = data.cells{n}.categories{c}.dFoT{tr};
            end
            
            if moving && transients
                movmask = reshape(data.metadata.categories{c}.moving{tr},vshape);
                transmask = reshape(data.cells{n}.categories{c}.transientmask{tr},vshape);
                thisdFoY = SBdiscretize(thistrace.*movmask.*transmask,thisy,bins);
                %thisdFodY = dFoSpeed(thistrace.*transmask,thisy,speedbins);
            elseif moving
                movmask = reshape(data.metadata.categories{c}.moving{tr},vshape);
                thisdFoY = SBdiscretize(thistrace.*movmask,thisy,bins);
                %thisdFodY = dFoSpeed(thistrace.*movmask,thisy,speedbins);
            elseif transients
                transmask = reshape(data.cells{n}.categories{c}.transientmask{tr},vshape);
                thisdFoY = SBdiscretize(thistrace.*transmask,thisy,bins);
                %thisdFodY = dFoSpeed(thistrace.*transmask,thisy,speedbins);
            else
                thisdFoY = SBdiscretize(thistrace,thisy,bins);
                %thisdFodY = dFoSpeed(thistrace,thisy,speedbins);
            end
            
            if zscored
                data.cells{n}.categories{c}.dZoY{tr} = thisdFoY;
                %data.cells{n}.categories{c}.dZodY{tr}=thisdFodY;
                % ---- TODO: Make compatible with new speed filtering!
            else
                data.cells{n}.categories{c}.dFoY{tr}=thisdFoY;
                %data.cells{n}.categories{c}.dFodY{tr}=thisdFodY;
                % ---- TODO: Make compatible with new speed filtering!
            end
        end
    end
end
            

function ds = recalculate_dFodY(ds,varargin)
% Up-to date version for interneuron data
p = inputParser;
addParameter(p,'transients',false,@islogical);
addParameter(p,'Yunit_to_cm',200,@isnumeric);
addParameter(p,'speedbins',0:31,@isnumeric);
parse(p,varargin{:})
transients = p.Results.transients;
Yunit_to_cm = p.Results.Yunit_to_cm;
speedbins = p.Results.speedbins;

for c = 1:length(ds.metadata.categories)
    for r = 1:length(ds.metadata.categories{c}.y)
        fr = ds.metadata.categories{c}.acquisition_rate(r);
        y = ds.metadata.categories{c}.y{r} * Yunit_to_cm;
        v = movmedian([0 abs(diff(movmedian(y,round(fr/3))))]*fr,fr);
        ds.metadata.categories{c}.v{r} = v;
        
        for n = 1:length(ds.cells)
            tr = ds.cells{n}.categories{c}.dFoT{r};
            Ztr = ds.cells{n}.categories{c}.zscored{r};
            
            if transients
                mask = ds.cells{n}.categories{c}.transientmask{r};
                tr = tr*mask;
                Ztr = Ztr*mask;
            end
            
            % CAVE: This refers to the CUSTOM DISCRETIZE FUNCTION (not
            % the Matlab function of the same name)!
            ds.cells{n}.categories{c}.dFodY{r} = SBdiscretize(tr,v,speedbins);
            ds.cells{n}.categories{c}.dZodY{r} = SBdiscretize(Ztr,v,speedbins);
        end
    end
end


end
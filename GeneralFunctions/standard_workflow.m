function data = standard_workflow(data,varargin)
% Take a dataset obtained with get_lineartracks and corrected category
% settings and apply the standard working procedures.
args=struct('bins',0.1:0.025:2.1,'minspeed',1.67E-3,'PVcells',false,... % For data coversion
    'sigmarange',2:0.5:5,'timerange',1:5,'cutoff',0.05,... % For event detection
    'nshuffle',1000,'ywindow',50,... % For spatial information
    'method','Dombeck','dFabsThresh',0.0,'dFrelThresh',7,'Ftransient',0.2,... % For placefields
    'sigma',1,'active',1/60,'pSI',0.05); % For place field correlations
    % will appear in the order specified in input argument.

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    args.(pair{1}) = pair{2};
end

fprintf('Converting data...\n');
data=convert_data(data,args.bins,args.minspeed,args.PVcells);
fprintf('Detecting events...\n');
data=significant_events(data,args.sigmarange,args.timerange,args.cutoff); 
    % Recalculate dFoY included in significant events
data=transientrates(data,true);
fprintf('Calculating spatial information...\n');
data=spatial_info_bootstrap(data,args.nshuffle,args.ywindow,args.bins);
fprintf('Detecting placefields...\n');
data=placefields(data,args.method,'bins',args.bins,'nshuffle',args.nshuffle,...
    'dFabsThresh',args.dFabsThresh,'dFrelThresh',args.dFrelThresh,...
    'Ftransient',args.Ftransient);
fprintf('Corrlating spatial information across sessions...\n');
data=spatial_corr(data,args.sigma,args.active,args.pSI);

end
function [dataset] = revise_eye(dataset, varargin)
% Rerun the lick detection. Should be run from the folder in which the
% '*.ephys' files are stored.

args=struct('folder',[],'display',false,'sensitivity',.92,'radii',[12 48],...
    'repeat_prep',false, 'separation',.5, 'frequency', 30);
% repeat preprocessing: set true, if previous values are insufficient
% separation: threshold for dividing pupil from iris pixels in % of total
% range
% frequency: recording freq, 30Hz for most recordings, for filtering

% Overwrite default parameters if required
for pair = reshape(varargin,2,[])
    if isfield(args,pair{1})
        args.(pair{1})=pair{2};
    else
        error('Input argument does not exist\n')
    end
end

if isfield(dataset.metadata,'eyemask') && isfield(dataset.metadata,'eyerange') && ~args.repeat_prep
    % Load preprocessing if exsitant
    mask = dataset.metadata.eyemask;
    range = dataset.metadata.eyerange;
else
    % Manual preprocessing
    for c=1:length(dataset.metadata.categories)
        fn = dataset.metadata.categories{c}.filename{3};
        vid = load([fn,'_eye.mat'],'data');
        vid = vid.data;
        
        mask(:,:,c) = roipoly(mean(vid,4));
        range(:,c) = input('Type min and max intensity value [min max]: ');
        close all;
    end
    
    dataset.metadata.eyemask = mask;
    dataset.metadata.eyerange = range;
end

% 1 Hz Lowpass to remove eyeblink artifacts
iir = designfilt('lowpassiir','FilterOrder',7, ...
    'HalfPowerFrequency',1,'SampleRate',30);

for c = 1:length(dataset.metadata.categories)
    for run = 1:length(dataset.metadata.categories{c}.ft)
        
        fn = dataset.metadata.categories{c}.filename{run};
        
        try
            vid = load([fn,'_eye.mat'],'data');
            vid = vid.data;
            
            close all;
            figure('Position',[700 500, 500, 350]);
            
%             centres = NaN(size(vid,4),2);
%             [radii, area] = deal(NaN(size(vid,4),1));
            area = deal(NaN(size(vid,4),1));
            
            for n=size(vid,4):-1:1
                % Area detection
                aim = vid(:,:,1,n).*uint8(mask(:,:,c));
                lb = prctile(aim(~aim==0),1);
                ub = prctile(aim(~aim==0),99);
                pupil = aim > (lb + args.separation*(ub-lb));
                area(n) = length(find(pupil));
                
%                 thisim = ((squeeze(vid(:,:,1,n)).*uint8(mask(:,:,c)))-range(1,c))...
%                     *(255/(range(2,c)-range(1,c)));
%                 [thiscent, thisrad] = imfindcircles(thisim,...
%                     args.radii,'ObjectPolarity','bright','Sensitivity',...
%                     args.sensitivity);
%                 if thiscent
%                     centres(n,:) = thiscent(1,:); % Accept only first detected circle
%                     radii(n) = thisrad(1);
%                 else
%                     centres(n,:) = [NaN, NaN];
%                     radii(n) = NaN;
%                 end
                
                % Display every 50th frame for control
                if args.display && ~mod(n,50)
%                     figure('Position',[100 500, 500, 350]); imshow(thisim); caxis([0 255]);
%                     h = viscircles(centres(n,:),radii(n)); 
                    imagesc(pupil); colormap('gray');
                    drawnow
                end
            end
            
            % Interpolate missing values
%             centres(:,1) = interp1(find(~isnan(centres(:,1))),...
%                 centres(~isnan(centres(:,1)),1), 1:length(centres(:,1)),...
%                 'nearest','extrap');
%             centres(:,2) = interp1(find(~isnan(centres(:,2))),...
%                 centres(~isnan(centres(:,2)),2), 1:length(centres(:,2)),...
%                 'nearest','extrap');
%             radii = interp1(find(~isnan(radii)), radii(~isnan(radii),1),...
%                 1:length(radii), 'nearest', 'extrap');
            area = interp1(find(~isnan(area)), area(~isnan(area),1),...
                1:length(area), 'nearest', 'extrap');
            area = filtfilt(iir,area);
            
            if args.display
                figure; hold on; 
                %plot(centres(:,1)); plot(centres(:,2)); plot(radii); plot(area/100);
                plot(area/100);
                drawnow
            end
            
            %dataset.metadata.categories{c}.pupil_radii{run} = radii;
            %dataset.metadata.categories{c}.pupil_centres{run} = centres;
            dataset.metadata.categories{c}.pupil_area{run} = area;
            
        catch
            warning('%s eye tracking failed',fn)
        end
    end
%    setup = false; % Re-adjust for every category
end
end
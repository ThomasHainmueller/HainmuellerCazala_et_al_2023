function plot_activity(neuron)
% Works for struc arrays (cells) returned by 'analyze'. Plots space as well
% as delay-period associated acitivity of an individual neuron over all
% runs.
SampLcorr = [];
SampLerr = [];
SampRcorr = [];
SampRerr = [];

DelLcorr = {};
DelLerr = {};
DelRcorr = {};
DelRerr = {};

MatLcorr = [];
MatLerr = [];
MatRcorr = [];
MatRerr = [];

% Append all traces to their appropriate fields
for tr = 1:length(neuron.trace)
    if strcmp(neuron.trace(tr).type,'sample')
        thistrace(1,:)=neuron.trace(tr).data;
        thistrace(2,:)=neuron.trace(tr).frametimes;
        thistrace(3,:)=neuron.trace(tr).xcord;
        thistrace(4,:)=neuron.trace(tr).ycord;
        try
            if strcmp(neuron.trace(tr).sample_dir,'left')
                if neuron.trace(tr).same_dir==0
                    SampLcorr = cat(2,SampLcorr,thistrace);
                elseif neuron.trace(tr).same_dir==1
                    SampLerr = cat(2,SampLerr,thistrace);
                end
            elseif strcmp(neuron.trace(tr).sample_dir,'right')
                if neuron.trace(tr).same_dir==0
                    SampRcorr = cat(2,SampRcorr,thistrace);
                elseif neuron.trace(tr).same_dir==1
                    SampRerr = cat(2,SampRerr,thistrace);
                end
            end
        clear thistrace;
        catch
            warning([num2str(tr) 'th trace couldnt be added']);
        end
    end    
    if strcmp(neuron.trace(tr).type,'delay')
        thistrace(1,:)=neuron.trace(tr).data;
        thistrace(2,:)=neuron.trace(tr).frametimes;
        try
            if strcmp(neuron.trace(tr).match_dir,'left')
                if neuron.trace(tr).same_dir==0
                    DelLcorr{end+1} = thistrace;
                elseif neuron.trace(tr).same_dir==1
                    DelLerr{end+1} = thistrace;
                end
            elseif strcmp(neuron.trace(tr).match_dir,'right')
                if neuron.trace(tr).same_dir==0
                    DelRcorr{end+1} = thistrace;
                elseif neuron.trace(tr).same_dir==1
                    DelRerr{end+1} = thistrace;
                end
            end
        clear thistrace;
        catch
            warning([num2str(tr) 'th trace couldnt be added']);
        end
    end    
    if strcmp(neuron.trace(tr).type,'match')
        thistrace(1,:)=neuron.trace(tr).data;
        thistrace(2,:)=neuron.trace(tr).frametimes;
        thistrace(3,:)=neuron.trace(tr).xcord;
        thistrace(4,:)=neuron.trace(tr).ycord;
        try
            if strcmp(neuron.trace(tr).match_dir,'left')
                if neuron.trace(tr).same_dir==0
                    MatLcorr = cat(2,SampLcorr,thistrace);
                elseif neuron.trace(tr).same_dir==1
                    MatLerr = cat(2,SampLerr,thistrace);
                end
            elseif strcmp(neuron.trace(tr).sample_dir,'right')
                if neuron.trace(tr).same_dir==0
                    MatRcorr = cat(2,SampRcorr,thistrace);
                elseif neuron.trace(tr).same_dir==1
                    MatRerr = cat(2,SampRerr,thistrace);
                end
            end
        clear thistrace;
        catch
            warning([num2str(tr) 'th trace couldnt be added']);
        end
    end
end

% Arrange Data for plots.
AllSamp = cat(2,SampLcorr,SampLerr,SampRcorr,SampRerr);
AllSampMap=map_data(AllSamp(1,:),AllSamp(3,:),AllSamp(4,:));
AllSampMap=AllSampMap(10:end,:);

AllMat = cat(2,MatLcorr,MatLerr,MatRcorr,MatRerr);
AllMatMap=map_data(AllMat(1,:),AllMat(3,:),AllMat(4,:));
AllMatMap=AllMatMap(10:end,:);

AllDel = [DelLcorr DelLerr DelRcorr DelRerr];

% Plot the data
subplot(4,6,[1 2 7 8])
imagesc(AllSampMap')
title('All Sample runs')

subplot(4,6,[3 4 9 10])
imagesc(AllMatMap')
title('All Match runs')

subplot(4,6,[5 6 11 12])
imagesc(trial_heatmap(AllDel))
title('All Delay periods')

subplot(4,6,13)
try
    SampLcorrMap=map_data(SampLcorr(1,:),SampLcorr(3,:),SampLcorr(4,:));
    imagesc(SampLcorrMap(12:end,:)')
end
title('Left sample before correct')

subplot(4,6,14)
try
    SampRcorrMap=map_data(SampRcorr(1,:),SampRcorr(3,:),SampRcorr(4,:));
    imagesc(SampRcorrMap(9:end,:)')
end
title('Right sample before correct')

subplot(4,6,15)
try
    MatLcorrMap=map_data(MatLcorr(1,:),MatLcorr(3,:),MatLcorr(4,:));
    imagesc(MatLcorrMap(12:end,:)')
end
title('Correct left matches')

subplot(4,6,16)
try
    MatRcorrMap=map_data(MatRcorr(1,:),MatRcorr(3,:),MatRcorr(4,:));
    imagesc(MatRcorrMap(9:end,:)')
end
title('Correct right matches')

subplot(4,6,17)
try
    imagesc(trial_heatmap(DelLcorr))
end
title('Delay before correct left')

subplot(4,6,18)
try
    imagesc(trial_heatmap(DelRcorr))
end
title('Delay before correct right')

subplot(4,6,19)
try
    SampLerrMap=map_data(SampLerr(1,:),SampLerr(3,:),SampLerr(4,:));
    imagesc(SampLerrMap(12:end,:)')
end
title('Left sample before errors')

subplot(4,6,20)
try
    SampRerrMap=map_data(SampRerr(1,:),SampRerr(3,:),SampRerr(4,:));
    imagesc(SampRerrMap(9:end,:)')
end
title('Right sample before errors')

subplot(4,6,21)
try
    MatLerrMap=map_data(MatLerr(1,:),MatLerr(3,:),MatLerr(4,:));
    imagesc(MatLerrMap(12:end,:)')
end
title('Erroneous left matches')

subplot(4,6,22)
try
    MatRerrMap=map_data(MatRerr(1,:),MatRerr(3,:),MatRerr(4,:));
    imagesc(MatRerrMap(9:end,:)')
end
title('Erroneous right matches')

subplot(4,6,23)
try
    imagesc(trial_heatmap(DelLerr))
end
title('Delay before left error')

subplot(4,6,24)
try
    imagesc(trial_heatmap(DelRerr))
end
title('Delay before right error')
end
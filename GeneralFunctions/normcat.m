function all = normcat(traces)
all = [];
for iter = [1:length(traces)]
    tr = traces{iter};
    tnorm = tr/mean(tr)-mean(tr/mean(tr));
    all = cat(1,all,tnorm);
end
return
end
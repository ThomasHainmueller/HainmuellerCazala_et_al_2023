function vectors = shuffle(vector,nsamples,window)
% Shuffle an input vector randomly for 'nsamples' times and return an
% output matrix with the shuffled vectors. Usefull for bootstrapping.
% 'window' defines the minimum size of the fragments to be shuffled.
cutpoints=1:window:length(vector);
cutpoints(end)=length(vector);

vectors=double.empty(length(vector),nsamples,0);

for n = 1:nsamples
    randorder = randperm(length(cutpoints)-1);
    index=1;
    for seg = 1:length(randorder)
        fragment=vector(cutpoints(randorder(seg)):cutpoints(randorder(seg)+1)-1);
        vectors(index:index+length(fragment)-1,n,1)=fragment';
        index=index+length(fragment);
    end
end
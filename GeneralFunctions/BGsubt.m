function stack = BGsubt(stack, Fbackground)
% Subtract the background value (the 'Fbackground' % of darkest pixels in
% the mean image from a 3D stack (xyt).

if nargin<2
    Fbackground=0.01;
end

hist = mean(stack,3);
hist = reshape(hist,(size(hist,1)*size(hist,2)),1);
hist = sort(hist);

BGval = mean(hist(1:round(Fbackground*length(hist))));

stack = stack-BGval;
stack(stack<0)=0;

end
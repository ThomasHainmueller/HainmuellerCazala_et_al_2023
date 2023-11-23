function a = gorconvert(gorobj)

for i = 1:length(gorobj),
[a(:,1,i),a(:,2,i)] = get(gorobj(i),'extract');
end
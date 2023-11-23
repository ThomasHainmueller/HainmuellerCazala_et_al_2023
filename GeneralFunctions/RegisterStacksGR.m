function registered = RegisterStacksGR(moving,fixed)
% Register two image stacks according to a reference image.
% moving, fixed = w x h x 2 matrices, 1=green and 2=red channel.

MovingRedM = mean(moving(:,:,2),3);
FixedRedM = mean(fixed(:,:,2),3);

[optimizer, metric]=imregconfig('monomodal');
optimizer.MaximumIterations = 5e3;
optimizer.MaximumStepLength = 1.0;

MovingRedReg = imregister(MovingRedM,FixedRedM,'translation',optimizer,metric);
tform = InferTranslation(MovingRedReg);
test=imtransform(MovingRedM,tform,'XData',[1 size(MovingRedM,2)],'YData',[1 size(MovingRedM,1)]);

green = moving(:,:,:,1);
red = moving(:,:,:,2);

registered(:,:,:,1)=imtransform(green,tform,'XData',[1 size(fixed,2)],'YData',[1 size(fixed,1)]);

registered(:,:,:,2)=imtransform(red,tform,'XData',[1 size(fixed,2)],'YData',[1 size(fixed,1)]);
end
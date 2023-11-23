function [greenreg redreg] = RegisterStacksRef(moving,movingRef,fixed,fixedRef)
% Register two image stacks according to a reference image.
% moving, fixed = w x h x 2 matrices, 1=green and 2=red channel.

[optimizer, metric]=imregconfig('monomodal');
optimizer.MaximumIterations = 1e5;
optimizer.MaximumStepLength = 0.8;

MovingRedReg = imregister(movingRef,fixedRef,'translation',optimizer,metric);
tform = InferTranslation(MovingRedReg);
%test=imtransform(MovingRedM,tform,'XData',[1 size(MovingRedM,2)],'YData',[1 size(MovingRedM,1)]);

greenreg=imtransform(moving,tform,'XData',[1 size(fixed,2)],'YData',[1 size(fixed,1)]);
clear moving;

redreg=imtransform(movingRef,tform,'XData',[1 size(fixed,2)],'YData',[1 size(fixed,1)]);
clear moving movingRef fixed fixedRef;
end
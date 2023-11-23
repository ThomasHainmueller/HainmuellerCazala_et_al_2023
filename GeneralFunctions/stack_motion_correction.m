function aligned = stack_motion_correction(stack, mode, winsize)
% Brute force variant, takes ages
if nargin < 3
    winsize = 50;
end

if nargin < 2
    mode = 'translation';
end

[optimizer, metric] = imregconfig('monomodal');

if strcmp(mode,'translation')
    %optimizer.MaximumIterations = 400;
    %optimizer.MaximumStepLength = .02;
    optimizer.MaximumIterations = 200;
    optimizer.MaximumStepLength = .02;
    optimizer.GradientMagnitudeTolerance = 3e-4;
    optimizer.MinimumStepLength = 3e-3;
    
    aligned(:,:,size(stack,3)) = stack(:,:,end);
    
    % First 50 frames
    for n=size(stack,3)-1:-1:size(stack,3)-winsize
        aligned(:,:,n) = imregister(stack(:,:,n),mean(stack(:,:,(n+1):end),3),...
            'translation',optimizer,metric); %replaced 'affine' 170513
    end
    
    % from thereon
    for n=size(stack,3)-winsize-1:-1:1
        aligned(:,:,n) = imregister(stack(:,:,n),mean(aligned(:,:,n+1:n+winsize+1),3),...
            'translation',optimizer,metric);
    end

elseif strcmp(mode,'affine')
    optimizer.MaximumIterations = 400;
    optimizer.MaximumStepLength = .02;
    
    aligned(:,:,size(stack,3)) = stack(:,:,end);
    
    % First 15 frames
    for n=size(stack,3)-1:-1:size(stack,3)-winsize
        aligned(:,:,n) = imregister(stack(:,:,n),mean(stack(:,:,(n+1):end),3),...
            'affine',optimizer,metric);
    end
    
    % from thereon
    for n=size(stack,3)-winsize-1:-1:1
        aligned(:,:,n) = imregister(stack(:,:,n),mean(stack(:,:,n+1:n+winsize+1),3),...
            'affine',optimizer,metric);
    end
end

end
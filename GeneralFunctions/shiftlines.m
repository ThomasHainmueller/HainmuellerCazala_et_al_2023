function corr = shiftlines(array,shift)
% for correction of forward-backward scan shift in framescans, shifts every
% second line in the scan by 'shift'; for movies or images.
if shift == 0
    corr = array;
    return
elseif shift >0
    for frame = size(array,3):-1:1
        for w=size(array,2)/2:-1:1
            % even numbered lines (2,4,6,...)
            corr(:,(w*2),frame)=array(shift:end,(w*2),frame);
            % odd numbered lines (1,3,5,...)
            corr(:,(w*2)-1,frame)=array(1:(size(array,1)+1-shift),(w*2)-1,frame);
        end
    end
    return
else
    for frame = size(array,3):-1:1
        for w=size(array,2)/2:-1:1
            % even numbered lines (2,4,6,...)
            corr(:,(w*2),frame)=array(1:size(array,1)+1+shift,(w*2),frame);
            % odd numbered lines (1,3,5,...)
            corr(:,(w*2)-1,frame)=array(-shift:end,(w*2)-1,frame);
        end
    end
    return    
end
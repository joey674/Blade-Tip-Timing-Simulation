% MDOF_SetBoundary get the data info and find best freq range, which fit
% the model best and ignore the rest of the data.
% the precess is like Weight.


function [freq,magn,weights_idx,peaks_idx] = SetBoundary(freq,magn,peaks_idx,weights_idx,boundary_idx)
    if ~isempty(peaks_idx)
        valid_idx = boundary_idx(1):boundary_idx(2);

        freq = freq(valid_idx);
        magn = magn(valid_idx);

        if weights_idx(1,1) < boundary_idx(1)
            weights_idx(1,1) = boundary_idx(1);
        end
        if weights_idx(end,2) > boundary_idx(2)
            weights_idx(end,2) = boundary_idx(2);
        end
        weights_idx = weights_idx -valid_idx(1) + 1; 
        peaks_idx = peaks_idx - valid_idx(1) + 1;
    end      
end


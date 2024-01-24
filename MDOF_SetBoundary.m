% MDOF_SetBoundary get the data info and find best freq range, which fit
% the model best and ignore the rest of the data.

% Input variables:
%   - freq: 
%   - magn: 
%   - err: 
%   - peaks_idx: 

% Notice:
%   - all vectors(freq,magn) will need to be long enough to have
%   boundary_idx in both direction

%   - all vector_idx(peaks_idx,weights_idx) can over the boundary(if is
%   overflow, will be cut to boundary)


function [freq,magn,weights_idx,peaks_idx] = MDOF_SetBoundary(freq,magn,peaks_idx,weights_idx,boundary_idx)
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


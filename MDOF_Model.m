% %% Model of MDOF(to calculate magn)
% function H_kl = MDOF_Model(P, w)
%     % cut the input back to n_modes vectors
%     n_modes = length(P) / 4;
%     ParVecArray = cell(1, n_modes);
%     for i = 1:n_modes
%         startIdx = (i-1)*4 + 1;  
%         endIdx = i*4;            
%         ParVecArray{i} = P(startIdx:endIdx);  
%     end
%     % assemble the model
%     H_kl = 0;
%     for i = 1:n_modes
%         P = ParVecArray{i};
%         W_m = P(1); D_m = P(2); r_re = P(3); r_im = P(4);
%         H_kl = H_kl + ...
%             (2 * (W_m * (r_re*D_m - r_im*sqrt(1-D_m*D_m)) + 2i*r_re*w))./...
%             (W_m*W_m - w.*w + 1i*2*D_m*W_m*w);
%     end
% end

%% Model of MDOF(to calculate magn)
function H_kl = MDOF_Model(P, w)
    % cut the input back to n_modes vectors
    n_modes = length(P) / 3;
    ParVecArray = cell(1, n_modes);
    for i = 1:n_modes
        startIdx = (i-1)*3 + 1;  
        endIdx = i*3;            
        ParVecArray{i} = P(startIdx:endIdx);  
    end
    % assemble the model
    H_kl = 0;
    for i = 1:n_modes
        P = ParVecArray{i};
        W_m = P(1); D_m = P(2); k_j = P(3);
        H_kl = H_kl + ...
            (k_j)./...
            (W_m*W_m - w.*w + 1i*2*D_m*W_m*w);
    end
end
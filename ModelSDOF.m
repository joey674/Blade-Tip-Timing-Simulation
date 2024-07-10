% Model of MDOF(to calculate magn) Model 2
function H_kl = ModelSDOF(P, w)
    W_m = P(1); D_m = P(2); k_j = P(3);
    H_kl = (k_j)./...
        (W_m*W_m - w.*w + 1i*2*D_m*W_m*w);
end


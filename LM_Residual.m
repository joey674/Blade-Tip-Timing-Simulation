
function residual = LM_Residual(params,freq,magn)
    residual = MDOF_Model(params,freq) - magn;
end
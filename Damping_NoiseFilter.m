function [output] = Damping_NoiseFilter(input)
    output = smoothdata(input,"rlowess",200);
end


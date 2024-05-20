function [output] = Damping_NoiseFilter(input)
    output = smoothdata(input,"sgolay",100);
end


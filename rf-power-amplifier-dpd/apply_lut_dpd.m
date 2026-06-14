function xPre = apply_lut_dpd(lut, x)

    amp = abs(x);

    xPre = zeros(size(x));

    for n = 1:length(x)

        idx = find(amp(n) >= lut.edges(1:end-1) & ...
                   amp(n) <  lut.edges(2:end), 1);

        if isempty(idx)
            gain = 1;
        else
            gain = lut.gain(idx);
        end

        xPre(n) = gain * x(n);
    end
end
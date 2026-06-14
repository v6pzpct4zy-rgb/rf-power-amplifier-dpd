function lut = train_lut_dpd(xTrain, ibo, cfg)

    % Number of LUT bins
    Nb = 32;

    % Input amplitude
    amp = abs(xTrain);

    % Reference linear signal
    xRef = linear_reference(xTrain, ibo);

    % PA output without DPD
    yPA = pa_model(xTrain, ibo, cfg);

    % Desired correction
    gainCorr = abs(xRef) ./ (abs(yPA) + eps);

    % LUT edges
    edges = linspace(min(amp), max(amp), Nb+1);

    lut.edges = edges;
    lut.gain  = zeros(Nb,1);

    for k = 1:Nb

        idx = amp >= edges(k) & amp < edges(k+1);

        if any(idx)
            lut.gain(k) = mean(gainCorr(idx));
        else
            lut.gain(k) = 1;
        end
    end
end
function dpd = train_dpd_model(xTrain, yTrain, ibo, cfg)
% Simple polynomial DPD:
% xPre = x + c3 |x|^2 x

    % Student-ID based parameter
    d = 18;   % sum of digits of student ID: 21070005021

    % IBO-dependent search range to reduce overcompensation
    if ibo >= 6
        c3_list = linspace(0, 0.1 + 0.002*d, 61);
    else
        c3_list = linspace(0, 0.3 + 0.01*d, 61);
    end

    bestNMSE = inf;
    best_c3 = 0;

    xRef = linear_reference(xTrain, ibo);

    for c3 = c3_list
        xPre = xTrain + c3*(abs(xTrain).^2).*xTrain;
        y = pa_model(xPre, ibo, cfg);

        nmseLin = sum(abs(y - xRef).^2) / sum(abs(xRef).^2);

        if nmseLin < bestNMSE
            bestNMSE = nmseLin;
            best_c3 = c3;
        end
    end

    dpd.c3 = best_c3;
end
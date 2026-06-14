function x_out = generate_modulated_waveform(Nsym, cfg)
% Generate pulse-shaped 16-QAM waveform

    data = randi([0 15], Nsym, 1);
    x = qammod(data, 16, 'UnitAveragePower', true);

    h = rcosdesign(cfg.rolloff, cfg.span, cfg.sps, 'sqrt');
    x_up = upfirdn(x, h, cfg.sps);

    delay = cfg.span * cfg.sps / 2;
    x_out = x_up(delay+1:end-delay);
end
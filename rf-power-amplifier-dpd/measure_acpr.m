function acpr_dB = measure_acpr(y, cfg)
% Simple ACPR estimate using PSD integration
% Main band: |f| <= 15 MHz
% Adjacent band: 20 MHz to 35 MHz

    nfft = 4096;
    [Pxx, f] = pwelch(y, [], [], nfft, cfg.fs, 'centered');

    mainMask = abs(f) <= 15e6;
    adjMask  = (abs(f) >= 20e6) & (abs(f) <= 35e6);

    Pmain = sum(Pxx(mainMask)) + eps;
    Padj  = sum(Pxx(adjMask)) + eps;

    acpr_dB = 10*log10(Padj / Pmain);
end
% P17: PA Nonlinearity, Polynomial DPD, LUT-Based DPD
% Final main.m with generalization test

rng(3787);
clear; clc; close all;

%% White figure background
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesColor','w');
set(groot,'defaultAxesXColor','k');
set(groot,'defaultAxesYColor','k');
set(groot,'defaultTextColor','k');

%% Configuration
cfg.inputBackoffdBList = [0 3 6 9];

cfg.nTrain = 5800;
cfg.nTest  = 2000;

cfg.rolloff = 0.25;
cfg.span = 10;
cfg.sps = 4;
cfg.fs = 100e6;

cfg.a1 = 1.0;
cfg.a3 = -0.15;
cfg.a5 = -0.05;

if ~exist('results','dir')
    mkdir('results');
end

Nibo = numel(cfg.inputBackoffdBList);

%% Pre-allocation
nmseNoDPD = zeros(Nibo,1);
nmsePoly  = zeros(Nibo,1);
nmseLUT   = zeros(Nibo,1);

acprNoDPD = zeros(Nibo,1);
acprPoly  = zeros(Nibo,1);
acprLUT   = zeros(Nibo,1);

evmNoDPD = zeros(Nibo,1);
evmPoly  = zeros(Nibo,1);
evmLUT   = zeros(Nibo,1);

c3Learned = zeros(Nibo,1);

%% Main loop
for iB = 1:Nibo

    ibo = cfg.inputBackoffdBList(iB);

    %% Training
    xTrain = generate_modulated_waveform(cfg.nTrain, cfg);

    dpdPoly = train_dpd_model(xTrain, [], ibo, cfg);
    c3Learned(iB) = dpdPoly.c3;

    lutDPD = train_lut_dpd(xTrain, ibo, cfg);

    %% Test
    xTest = generate_modulated_waveform(cfg.nTest, cfg);
    xRef  = linear_reference(xTest, ibo);

    %% No DPD
    yNo = pa_model(xTest, ibo, cfg);

    %% Polynomial DPD
    xPrePoly = apply_dpd(dpdPoly, xTest, cfg);
    yPoly = pa_model(xPrePoly, ibo, cfg);

    %% LUT DPD
    xPreLUT = apply_lut_dpd(lutDPD, xTest);
    yLUT = pa_model(xPreLUT, ibo, cfg);

    %% Metrics
    nmseNoDPD(iB) = compute_nmse(yNo, xRef);
    nmsePoly(iB)  = compute_nmse(yPoly, xRef);
    nmseLUT(iB)   = compute_nmse(yLUT, xRef);

    acprNoDPD(iB) = measure_acpr(yNo, cfg);
    acprPoly(iB)  = measure_acpr(yPoly, cfg);
    acprLUT(iB)   = measure_acpr(yLUT, cfg);

    %% EVM
    xSym     = xRef(1:cfg.sps:end);
    yNoSym   = yNo(1:cfg.sps:end);
    yPolySym = yPoly(1:cfg.sps:end);
    yLUTSym  = yLUT(1:cfg.sps:end);

    L = min([length(xSym), length(yNoSym), length(yPolySym), length(yLUTSym)]);

    xSym     = xSym(1:L);
    yNoSym   = yNoSym(1:L);
    yPolySym = yPolySym(1:L);
    yLUTSym  = yLUTSym(1:L);

    gNo   = (xSym' * yNoSym)   / (xSym' * xSym);
    gPoly = (xSym' * yPolySym) / (xSym' * xSym);
    gLUT  = (xSym' * yLUTSym)  / (xSym' * xSym);

    yNoAligned   = yNoSym   / gNo;
    yPolyAligned = yPolySym / gPoly;
    yLUTAligned  = yLUTSym  / gLUT;

    evmNoDPD(iB) = sqrt(mean(abs(yNoAligned - xSym).^2) / mean(abs(xSym).^2)) * 100;
    evmPoly(iB)  = sqrt(mean(abs(yPolyAligned - xSym).^2) / mean(abs(xSym).^2)) * 100;
    evmLUT(iB)   = sqrt(mean(abs(yLUTAligned - xSym).^2) / mean(abs(xSym).^2)) * 100;

    %% Detailed plots at IBO = 3 dB
    if ibo == 3

        fprintf('IBO = 3 dB | EVM No DPD = %.4f %% | Poly DPD = %.4f %% | LUT DPD = %.4f %%\n', ...
            evmNoDPD(iB), evmPoly(iB), evmLUT(iB));
        fprintf('IBO = 3 dB | Learned c3 = %.4f\n', dpdPoly.c3);

        %% AM/AM
        amp_in   = abs(xRef);
        amp_no   = abs(yNo);
        amp_poly = abs(yPoly);
        amp_lut  = abs(yLUT);

        [amp_in_sorted, idx] = sort(amp_in);

        figure;
        plot(amp_in_sorted, amp_in_sorted, '--k', 'LineWidth', 1.5); hold on;
        plot(amp_in_sorted, amp_no(idx), 'b', 'LineWidth', 1.3);
        plot(amp_in_sorted, amp_poly(idx), 'r', 'LineWidth', 1.3);
        plot(amp_in_sorted, amp_lut(idx), 'g', 'LineWidth', 1.3);
        grid on;
        xlabel('Input amplitude |x|');
        ylabel('Output amplitude |y|');
        title('AM/AM Before and After DPD (IBO = 3 dB)');
        legend('Ideal Linear','No DPD','Polynomial DPD','LUT DPD','Location','best');
        saveas(gcf, fullfile('results','fig_AMAM_IBO3.png'));

        %% PSD
        nfft = 2048;
        [pxx_ref, f]  = pwelch(xRef, [], [], nfft, cfg.fs, 'centered');
        [pxx_no, ~]   = pwelch(yNo, [], [], nfft, cfg.fs, 'centered');
        [pxx_poly, ~] = pwelch(yPoly, [], [], nfft, cfg.fs, 'centered');
        [pxx_lut, ~]  = pwelch(yLUT, [], [], nfft, cfg.fs, 'centered');

        figure;
        plot(f/1e6, 10*log10(pxx_ref + eps), 'k--', 'LineWidth', 1.2); hold on;
        plot(f/1e6, 10*log10(pxx_no + eps), 'b', 'LineWidth', 1.2);
        plot(f/1e6, 10*log10(pxx_poly + eps), 'r', 'LineWidth', 1.2);
        plot(f/1e6, 10*log10(pxx_lut + eps), 'g', 'LineWidth', 1.2);
        grid on;
        xlabel('Frequency (MHz)');
        ylabel('PSD (dB/Hz)');
        title('PSD Before and After DPD (IBO = 3 dB)');
        legend('Reference','No DPD','Polynomial DPD','LUT DPD','Location','best');
        saveas(gcf, fullfile('results','fig_PSD_IBO3.png'));

        %% PSD zoom
        figure;
        plot(f/1e6, 10*log10(pxx_no + eps), 'b', 'LineWidth', 1.2); hold on;
        plot(f/1e6, 10*log10(pxx_poly + eps), 'r', 'LineWidth', 1.2);
        plot(f/1e6, 10*log10(pxx_lut + eps), 'g', 'LineWidth', 1.2);
        grid on;
        xlim([18 35]);
        xlabel('Frequency (MHz)');
        ylabel('PSD (dB/Hz)');
        title('PSD Zoom-In Around Adjacent Band (IBO = 3 dB)');
        legend('No DPD','Polynomial DPD','LUT DPD','Location','best');
        saveas(gcf, fullfile('results','fig_PSD_zoom_IBO3.png'));

        %% Constellation
        Ns = min(800, length(xSym));

        figure;

        subplot(1,3,1);
        plot(real(yNoAligned(1:Ns)), imag(yNoAligned(1:Ns)), '.', 'MarkerSize', 5); hold on;
        plot(real(xSym(1:Ns)), imag(xSym(1:Ns)), 'ko', 'MarkerSize', 4);
        grid on; axis equal;
        xlabel('In-Phase');
        ylabel('Quadrature');
        title(sprintf('No DPD\nEVM = %.4f %%', evmNoDPD(iB)));
        legend('Received','Ideal','Location','best');

        subplot(1,3,2);
        plot(real(yPolyAligned(1:Ns)), imag(yPolyAligned(1:Ns)), '.', 'MarkerSize', 5); hold on;
        plot(real(xSym(1:Ns)), imag(xSym(1:Ns)), 'ko', 'MarkerSize', 4);
        grid on; axis equal;
        xlabel('In-Phase');
        ylabel('Quadrature');
        title(sprintf('Polynomial DPD\nEVM = %.4f %%', evmPoly(iB)));
        legend('Received','Ideal','Location','best');

        subplot(1,3,3);
        plot(real(yLUTAligned(1:Ns)), imag(yLUTAligned(1:Ns)), '.', 'MarkerSize', 5); hold on;
        plot(real(xSym(1:Ns)), imag(xSym(1:Ns)), 'ko', 'MarkerSize', 4);
        grid on; axis equal;
        xlabel('In-Phase');
        ylabel('Quadrature');
        title(sprintf('LUT DPD\nEVM = %.4f %%', evmLUT(iB)));
        legend('Received','Ideal','Location','best');

        sgtitle('Constellation Comparison at IBO = 3 dB');
        saveas(gcf, fullfile('results','fig_constellation_IBO3.png'));
    end
end

%% Summary NMSE and ACPR
figure;

subplot(2,1,1);
plot(cfg.inputBackoffdBList, nmseNoDPD, '-o', 'LineWidth', 1.2); hold on;
plot(cfg.inputBackoffdBList, nmsePoly, '-s', 'LineWidth', 1.2);
plot(cfg.inputBackoffdBList, nmseLUT, '-d', 'LineWidth', 1.2);
grid on;
xlabel('Input backoff [dB]');
ylabel('NMSE [dB]');
legend('No DPD','Polynomial DPD','LUT DPD','Location','best');
title('NMSE versus IBO');

subplot(2,1,2);
plot(cfg.inputBackoffdBList, acprNoDPD, '-o', 'LineWidth', 1.2); hold on;
plot(cfg.inputBackoffdBList, acprPoly, '-s', 'LineWidth', 1.2);
plot(cfg.inputBackoffdBList, acprLUT, '-d', 'LineWidth', 1.2);
grid on;
xlabel('Input backoff [dB]');
ylabel('ACPR [dB]');
legend('No DPD','Polynomial DPD','LUT DPD','Location','best');
title('ACPR versus IBO');

saveas(gcf, fullfile('results','fig_NMSE_ACPR_vs_IBO.png'));

%% EVM
figure;
plot(cfg.inputBackoffdBList, evmNoDPD, '-o', 'LineWidth', 1.2); hold on;
plot(cfg.inputBackoffdBList, evmPoly, '-s', 'LineWidth', 1.2);
plot(cfg.inputBackoffdBList, evmLUT, '-d', 'LineWidth', 1.2);
grid on;
xlabel('Input backoff [dB]');
ylabel('EVM [%]');
legend('No DPD','Polynomial DPD','LUT DPD','Location','best');
title('EVM versus IBO');
saveas(gcf, fullfile('results','fig_EVM_vs_IBO.png'));

%% Learned c3
figure;
plot(cfg.inputBackoffdBList, c3Learned, '-o', 'LineWidth', 1.2);
grid on;
xlabel('Input backoff [dB]');
ylabel('Learned c_3');
title('Learned Polynomial DPD Coefficient versus IBO');
saveas(gcf, fullfile('results','fig_c3_vs_IBO.png'));

%% Generalization Test
% Train both DPD models only at IBO = 3 dB and test them at all IBO values.
% This evaluates whether a DPD model trained at one operating point can
% still compensate the PA when the input backoff changes.
trainIBO = 3;

xTrainGen = generate_modulated_waveform(cfg.nTrain, cfg);

dpdPolyGen = train_dpd_model(xTrainGen, [], trainIBO, cfg);
lutDPDGen  = train_lut_dpd(xTrainGen, trainIBO, cfg);

genNMSEPoly = zeros(Nibo,1);
genNMSELUT  = zeros(Nibo,1);

genEVMPoly = zeros(Nibo,1);
genEVMLUT  = zeros(Nibo,1);

genACPRPoly = zeros(Nibo,1);
genACPRLUT  = zeros(Nibo,1);

for iB = 1:Nibo

    testIBO = cfg.inputBackoffdBList(iB);

    %% Test waveform for generalization
    xTestGen = generate_modulated_waveform(cfg.nTest, cfg);
    xRefGen  = linear_reference(xTestGen, testIBO);

    %% Polynomial DPD trained at IBO = 3 dB
    xPrePolyGen = apply_dpd(dpdPolyGen, xTestGen, cfg);
    yPolyGen = pa_model(xPrePolyGen, testIBO, cfg);

    %% LUT DPD trained at IBO = 3 dB
    xPreLUTGen = apply_lut_dpd(lutDPDGen, xTestGen);
    yLUTGen = pa_model(xPreLUTGen, testIBO, cfg);

    %% NMSE
    genNMSEPoly(iB) = compute_nmse(yPolyGen, xRefGen);
    genNMSELUT(iB)  = compute_nmse(yLUTGen, xRefGen);

    %% ACPR
    genACPRPoly(iB) = measure_acpr(yPolyGen, cfg);
    genACPRLUT(iB)  = measure_acpr(yLUTGen, cfg);

    %% EVM
    xSymGen     = xRefGen(1:cfg.sps:end);
    yPolySymGen = yPolyGen(1:cfg.sps:end);
    yLUTSymGen  = yLUTGen(1:cfg.sps:end);

    Lgen = min([length(xSymGen), length(yPolySymGen), length(yLUTSymGen)]);

    xSymGen     = xSymGen(1:Lgen);
    yPolySymGen = yPolySymGen(1:Lgen);
    yLUTSymGen  = yLUTSymGen(1:Lgen);

    gPolyGen = (xSymGen' * yPolySymGen) / (xSymGen' * xSymGen);
    gLUTGen  = (xSymGen' * yLUTSymGen)  / (xSymGen' * xSymGen);

    yPolyAlignedGen = yPolySymGen / gPolyGen;
    yLUTAlignedGen  = yLUTSymGen  / gLUTGen;

    genEVMPoly(iB) = sqrt(mean(abs(yPolyAlignedGen - xSymGen).^2) / mean(abs(xSymGen).^2)) * 100;
    genEVMLUT(iB)  = sqrt(mean(abs(yLUTAlignedGen  - xSymGen).^2) / mean(abs(xSymGen).^2)) * 100;
end

%% Generalization EVM plot
figure;
plot(cfg.inputBackoffdBList, genEVMPoly, '-s', 'LineWidth', 1.2); hold on;
plot(cfg.inputBackoffdBList, genEVMLUT, '-d', 'LineWidth', 1.2);
grid on;
xlabel('Test IBO [dB]');
ylabel('EVM [%]');
legend('Polynomial DPD trained at 3 dB', 'LUT DPD trained at 3 dB', 'Location','best');
title('Generalization Test: EVM when DPD is Trained at IBO = 3 dB');
saveas(gcf, fullfile('results','fig_generalization_EVM.png'));

%% Generalization NMSE plot
figure;
plot(cfg.inputBackoffdBList, genNMSEPoly, '-s', 'LineWidth', 1.2); hold on;
plot(cfg.inputBackoffdBList, genNMSELUT, '-d', 'LineWidth', 1.2);
grid on;
xlabel('Test IBO [dB]');
ylabel('NMSE [dB]');
legend('Polynomial DPD trained at 3 dB', 'LUT DPD trained at 3 dB', 'Location','best');
title('Generalization Test: NMSE when DPD is Trained at IBO = 3 dB');
saveas(gcf, fullfile('results','fig_generalization_NMSE.png'));

%% Generalization results table
generalizationTable = table(cfg.inputBackoffdBList(:), ...
                            genNMSEPoly, genNMSELUT, ...
                            genEVMPoly, genEVMLUT, ...
                            genACPRPoly, genACPRLUT, ...
    'VariableNames', {'Test_IBO_dB', ...
                      'NMSE_PolyDPD_TrainedAt3dB_dB', ...
                      'NMSE_LUTDPD_TrainedAt3dB_dB', ...
                      'EVM_PolyDPD_TrainedAt3dB_percent', ...
                      'EVM_LUTDPD_TrainedAt3dB_percent', ...
                      'ACPR_PolyDPD_TrainedAt3dB_dB', ...
                      'ACPR_LUTDPD_TrainedAt3dB_dB'});

disp(generalizationTable);
writetable(generalizationTable, fullfile('results','generalization_test_trained_IBO3.csv'));


%% Results table
resultsTable = table(cfg.inputBackoffdBList(:), ...
                     nmseNoDPD, nmsePoly, nmseLUT, ...
                     evmNoDPD, evmPoly, evmLUT, ...
                     acprNoDPD, acprPoly, acprLUT, ...
                     c3Learned, ...
    'VariableNames', {'IBO_dB', ...
                      'NMSE_NoDPD_dB', 'NMSE_PolyDPD_dB', 'NMSE_LUTDPD_dB', ...
                      'EVM_NoDPD_percent', 'EVM_PolyDPD_percent', 'EVM_LUTDPD_percent', ...
                      'ACPR_NoDPD_dB', 'ACPR_PolyDPD_dB', 'ACPR_LUTDPD_dB', ...
                      'Learned_c3'});

disp(resultsTable);

writetable(resultsTable, fullfile('results','results_table_P17_final.csv'));

%% Save all
save(fullfile('results','results_P17_final.mat'), ...
     'cfg', ...
     'nmseNoDPD', 'nmsePoly', 'nmseLUT', ...
     'acprNoDPD', 'acprPoly', 'acprLUT', ...
     'evmNoDPD', 'evmPoly', 'evmLUT', ...
     'c3Learned', 'resultsTable', ...
     'genNMSEPoly', 'genNMSELUT', ...
     'genEVMPoly', 'genEVMLUT', ...
     'genACPRPoly', 'genACPRLUT', ...
     'generalizationTable');

fprintf('\nFinal P17 simulation with generalization test completed. Results saved in the results folder.\n');
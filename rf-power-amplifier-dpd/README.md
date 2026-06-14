# RF Power Amplifier Digital Predistortion

This project presents a MATLAB-based simulation framework for modeling nonlinear RF power amplifier behavior and improving linearity using digital predistortion techniques.
The system uses 16-QAM modulation with root-raised cosine pulse shaping and evaluates the performance of polynomial-based and LUT-based DPD methods under different input back-off conditions.

## Key Features
- 16-QAM modulated waveform generation
- Root-raised cosine pulse shaping
- Memoryless nonlinear RF power amplifier modeling
- Polynomial-based digital predistortion
- LUT-based digital predistortion
- NMSE, EVM, and ACPR performance evaluation
- Constellation and PSD-based visualization

## Tools
- Communications Toolbox
- Signal Processing Toolbox

## Project Structure

main.m                         Main simulation script
generate_modulated_waveform.m  Generates the 16-QAM waveform
pa_model.m                     Nonlinear RF power amplifier model
train_dpd_model.m              Trains polynomial DPD coefficient
apply_dpd.m                    Applies polynomial DPD
train_lut_dpd.m                Trains LUT-based DPD
apply_lut_dpd.m                Applies LUT-based DPD
compute_nmse.m                 Computes NMSE performance
measure_acpr.m                 Measures ACPR
results/                       Simulation outputs and plots
report.pdf                     Project report
The simulation automatically generates figures and saves all outputs into the results/ folder.

Approximate runtime: 1–2 minutes.

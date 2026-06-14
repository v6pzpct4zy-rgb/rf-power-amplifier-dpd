function nmse_dB = compute_nmse(y, xRef)

    nmseLin = sum(abs(y - xRef).^2) / sum(abs(xRef).^2);
    nmse_dB = 10*log10(nmseLin);
end
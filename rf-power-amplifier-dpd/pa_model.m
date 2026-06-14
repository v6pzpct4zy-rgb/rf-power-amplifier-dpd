function y = pa_model(x, ibo_dB, cfg)
% Memoryless nonlinear PA

    x_in = x * 10^(-ibo_dB/20);

    y = cfg.a1*x_in + ...
        cfg.a3*(abs(x_in).^2).*x_in + ...
        cfg.a5*(abs(x_in).^4).*x_in;
end
function xPre = apply_dpd(dpd, x, cfg)

    xPre = x + dpd.c3*(abs(x).^2).*x;
end
function xRef = linear_reference(x, ibo_dB)

    xRef = x * 10^(-ibo_dB/20);
end
function T=getTFromDelta(delta)
    % Computes the electron temperature (in eV) from the
    % delta (=v_th/c) parameter used in CODE.
    
    e = 1.602176e-19;
    c = 2.997925e8; %m/s
    m0 = 9.10938e-31; % kg
    
    T = (delta * c) ^ 2 * m0 / (2 * e);
end
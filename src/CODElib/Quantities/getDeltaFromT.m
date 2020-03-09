function delta = getDeltaFromT(T)
    % Computes the delta (=v_th/c) parameter used in CODE from the
    % electron temperature (in eV).

    e = 1.602176e-19;
    c = 2.997925e8; %m/s
    m0 = 9.10938e-31; % kg

    v_th = sqrt(2*e*T/m0);
    delta = v_th/c;
end
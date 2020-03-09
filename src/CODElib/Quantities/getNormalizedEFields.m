function [EOverEC,EOverED,EHat] = getNormalizedEFields(E,T,density)
    %   Normalizes an electric field in V/m to the critical field EOverEc,
    %   and to the dreicer field EOverED
    % E is the field in V/m, T is the electron temperature in eV and
    % density is the electron density in m^{-3}.

    e = 1.602176e-19;
    eCube = e^3;
    eps0Sq = (8.854188e-12)^2;
    cSq = (2.997925e8)^2; %m/s %%A: really in m²/s², but c is expressed in m/s
    m0 = 9.10938e-31; %kg

    couLog = 14.9-0.5*log(density/1e20)+log(T/1e3);
    delta = getDeltaFromT(T);
    eCritical = density.*couLog*eCube / (4*pi*eps0Sq*m0*cSq);
    EOverEC = E./eCritical;
    EOverED = e*T/(m0*cSq) .* EOverEC;
    EHat = EOverEC .* delta .* delta * 3 * sqrt(pi) / 4;
end
function [delta,lnLambda,collFreq,BHat] ...
            = getDerivedParameters(T,density,B)
    % Calculates the parameters delta, coloumb logarithm lnLambda and the
    % collision frequency (and time) for given E (in V/m), T (in eV) and density
    % (in m^{-3}). Collission frequency uses Matt's definition
    % Collision time = 1/Collision Frequency

    e = 4.80320425e-10; %Gaussian units!
    m0 = 9.10938e-28; %Gaussian
    c = 2.997925e+10; %Gaussian

    delta = getDeltaFromT(T);
    lnLambda = 14.9-0.5*log(density/1e20)+log(T/1e3); %T in eV
    TG = 1.602176e-12 * T; %Gaussian convertion
    densityG = density * 1e-6; %Gaussian convertion
    collFreq = 4/3*sqrt(2*pi) * e^4 *densityG .* lnLambda / sqrt(m0) .* TG.^(-1.5); %From Matt's definition of collision time
    %         collisionTime = 1/collisionFreq;
    if isempty(B)
        BHat = 0;
    else
        BHat = sqrt( 2*e^4./(3*m0^3*c^5.*collFreq) ) .* 1e4.*B;
    end
end
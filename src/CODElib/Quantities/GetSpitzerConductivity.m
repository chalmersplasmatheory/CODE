function sigma = GetSpitzerConductivity(T,n,Z)
    %GETSPITZERCONDUCTIVITY calculates the plasma Spitzer conductivity in 
    %units of (Ohm m)^-1, from the formula on p. 72 of Helander & Sigmar. 
    %Usage:
    %       sigma = GetSpitzerConductivity(T,n,Z)
    %
    % n must be in units of m^(-3), T in eV.
    
    %There is a prefactor that depends on Z which is only known numerically
    %(table on p.74 in Helander & Sigmar and Table III in Spitzer & Härm).
    %Let's interpolate to find the value for any Z. Since there is a data
    %point at infinity, lets use 1/Z for the interpolation.
    invZs = [1,0.5,0.25,1/16,0];
    vals = 32/(3*pi)*[0.5816,0.6833,0.7849,0.9225,1];    
    preFac = interp1(invZs,vals,1./Z,'pchip');

    lnLambda = 14.9-0.5*log(n/1e20)+log(T/1e3); 
    constants = 9.695919807026154e3; % = 12/sqrt(2)*pi^1.5*eps0^2/(sqrt(e)*sqrt(m0))
    sigma = preFac*constants*T.^(3/2)./(Z.*lnLambda); %Helander and Sigmar

    %Simplified expression in Chang (Eq. [5-76]):
%     sigma = 19.2307e3*T.^(1.5)./Z./lnLambda;

%Written by Adam Stahl
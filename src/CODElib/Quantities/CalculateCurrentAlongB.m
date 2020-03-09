function i = CalculateCurrentAlongB(f,y,yWeights, Ny, nRef ,deltaRef)
%CALCULATECURRENTALONGB Calculates current density in SI along magnetic
%axis
%specifically calculates integral^3 f hat n e v d^3v
% which reduces to
% 4 nRef e VRef / (3 sqrtPi m^3) integral_0^inf y^3 / gamma dy
%where e is electron charge, and v is the electron speed

% %%%%%%%%%%%%%%
% Input:
% %%%%%%%%%%%%%%%
% f - distribution function with normalization as in CODE
% and Nxi legandre modes projections. first Ny points are for first
% legandre mode.
% y - is vector of momentum points as y = gamma v / v_Ref where gamma is
%   relativistic gamma function, v electron speed, and v_Ref reference
%   thermal speed of electrons used in y
% yWeights - contains dy weights for integrating over y
% nRef - referece density of electrons in SI
% deltaRef - vRef/c where c is speed of light

% Written by Albert Johansson 2019

e = -1.60217662e-19; %SI
c = 299792458; %SI

f1 = f(Ny+1:2*Ny);
f1 = f1(:);
y = y(:);
yWeights = yWeights(:)';

gamma = sqrt(1+deltaRef^2*y.^2);
preFac = 4 * e * nRef * deltaRef * c / (3* sqrt(pi));

i = preFac * yWeights * (f1 .* y.^3 ./ gamma);
end


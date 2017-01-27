function [L, vo, r1, r2, V, T, Et, S] = lambda(Ealpha, A, P)
% Calculates the decay constant from the formula lambda = P*fo*exp(-2S),
% where P is the preformation factor (probability that an alpha particle
% exists at a certain time), fo is the frequency with which the alpha
% particle hits the potential barrier and S is the WKB- factor. fo is
% evaluated as fo = sqrt(2*e/m)/(2*r1), where e = Ealpha-(Vo+V_em(r=0)) is
% an approximation of the kinetic energy of the alpha particle in the
% nucleus and m is the alpha particles' mass. 
%
% Input arguments:
% Ealpha: Energy of the alpha particle
% A: Mass number of the DAUGHTER nucleus
% P: Preformation factor
%
% Output arguments:
% L: Decay constant
% vo: Vo term of the Woods-Saxon potential
% r1: Beginning of the tunnel
% r2: End of the tunnel
% V: Array of values for V(r) from 0.01->100 fm
% T: Array of the energy levels, bound and first quasibound
% Et: First quasibound energylevel
% S: WKB- factor

tic 
Z = 92; % Proton number of the parent nucleus (Uranium)

% Calculate the WKB- factor by calling sigma()
[S, vo, r1, r2, V, T, Et] = sigma(Ealpha, A);

V_em0 = 3*(Z-2)*1.4399764*(1.3*(A^(1/3))^(-1)); % Coulomb potential at r=0, units MeV
e = Ealpha-(vo+V_em0);                          % Kinetic energy of the alpha particle, units MeV
c = 299792458*10^(15);                          % Speed of light in units fm/s 
m = (3.72737924*10^3)/(c^2);                    % Mass of the alpha particle, units MeVs^2/(m^2)
fo = sqrt(2*e/m)/(2*r1);                        % Frequency by which the alpha particle hits the potential barrier

L = P*fo*exp(-2*S); % Decay constant
toc
end



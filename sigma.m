function [S, vo, r1, r2, X, T, Et]  = sigma(Ealpha, A)
% Calculates the WKB- factor, i.e the tunneling integral. This function
% calls the function Vo() in order to determine the Vo term for the
% Woods-Saxon potential needed in the function radii(), which is also
% called in order to determine the values of the potential energy in the
% interval of the tunnel, and the values r1, and r2 for the ends of the
% tunnel.
%
% Input arguments
% Ealpha: The alpha particle energy
% A: Mass number of the DAUGHTER nucleus
%
% Output arguments:
% S: The WKB- factor
% vo: The Vo term of the Woods-Saxon potential
% r1: Beginning of the tunnel
% r2: End of the tunnel
% X: Array of the values of the potential energy has between 0.01->100 fm
% T: Array of the energy levels, bound and the first quasibound
% Et: Energy of the first positive (quasibound) energy level

H = 1/(6.582119514*10^(-22)); % 1 over the reduced Planck constant, units (MeV*s)^-1
c = 299792458*10^(15);        % Speed of light in the units fm/s 
M = 931.494095*(4*A/((A+4)*(c^2)));  % Reduced mass of the alpha particle in units MeVs^2/(fm^2)

% Determine the Vo term of Woods-Saxon potetnial by calling Vo() with the
% beginning value of -100 MeV for Vo
V0 = -100;           
[vo, T, Et] = Vo(Ealpha, A, V0);

% Determine the values of V(r) in the interval between the tunnel, and the
% ends of the tunnel r1 and r2
[r1, r2, V, X] = radii(Ealpha, A, vo);

% Calculate the WKB- factor by numerical integration using trapetzoidal
% method
x = sqrt(2*M*abs(V-Ealpha)); % The integrand
z = (r1:0.01:r2);            % The integration interval
S = H*trapz(z,x);          % WKB- factor

end

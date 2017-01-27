function [V, T, Et] = Vo(Ealpha, A, V0)
% Calculates the Vo term for the Woods-Saxon potential. Also gives out an
% array of all the energy levels including the quasi bound energy level,
% and the alpha particle energy (this so one can check that the calculation
% is correct). This function calls two other functions, written below as is
% done in MATLAB, called Elevel() and schr(); details of these below. This
% code is written for even isotopes of Uranium from U(238) to U(222), so
% incase one wants to use this for other isotopes, modification is
% necessary.
%
% Input arguments:
% Ealpha: Measured energy of the alpha particle
% A: Mass number of the DAUGHTER nucleus
% V0: Some beginnign value for Vo, need to start the search somewhere
%
% Output arguments:
% V: The Vo term calculated below
% T: Array of calculated energy levels
% Et: Energy of the quasi bound state, should be as close as possible to
%     the alpha particle energy

% Defining variables
l = 0;   % Variable for the while-loop
dV = -1; % Starting value for the change in Vo term
Et = 0;  % Quasi bound energy level from the previous iteration, starting value set to zero
n1 = 10; % Desired quantum number based on previous code runs, this is to cut down on computing time

while l<1
    
    % Calculate the energy levels up to the first positive (quasibound)
    % level, the quantum number n and array of all energy levels
    [Ea, n, T] = Elevel(A, V0);
    
    % If n is too large:
    if n>n1
        V0 = V0+2; % Lift the bottom of the potential well in order to decrease n 
    end
    
    % If n is too small:
    if n<n1
        Et = Ea;
        V0 = V0-2; % Push down the bottom of the potential well in order to increase n
    end
    
    % If n is exactly what we want:
    if n==n1
        
        % If  the new and old calculated quasi bound energy levels on the 
        % same side of alpha particle energy, change V0 in order to get Ea
        % colser to Ealpha:
        if (Ealpha-Ea)*(Ealpha-Et)>0
            Et = Ea;    % Update Et, the previous calculated (quasi)bound energy level
            V0 = V0+dV; % Update V0
            
        % Or if the are on different sides, we crossed the desired value
        % for V0. In that case:
        elseif (Ealpha-Ea)*(Ealpha-Et)<0
            
            % In case the energy Et is smaller than Ealpha, meaning Ea, the
            % energylevel now calculated is slightly above Ealpha:
            if (Et<Ealpha)  
                    
                % If we have not reached the desired accuracy:
                if abs(dV)>10^(-6)
                    Et = Ea;          % Update Et
                    dV = (dV/2)*(-1); % Update dV
                    V0 = V0+dV;       % And update V0
                    
                % If we have reached the desired accuracy
                elseif ((n1==n) && (abs(dV)<10^(-6)))
                    Et = Ea; 
                    V = V0;  % V0 has the value that reproduces the measured alpha particle energy
                    l = 2;   % Break the while loop
                    
                end
            
            % Or if Et is above Ealpha, meaning Ea is below it then:
            elseif Et>Ealpha
                Et = Ea;
                dV = (dV/2)*(-1);
                V0 = V0+dV;
                
            end
        end
            
    end
end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% The functions that are calle by the main function Vo().
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function r1 = schr(E, A, V0)
% This function iterates the discretized Schrödingers equation
% u[i+1]=2u[i]+2m(d/h)^2(V[i]-E)u[i]-u[i-1], where h is the reduced plancks
% constant, d is length of the iteration step and m is mass of the alpha 
% particle, to the distance of 20 fm. 
%
% Input arguments:
% E: Energy level of the state
% A: Mass number of the DAUGHTER nucleus
% V0: Vo term of the Woods-Saxon potential
%
% Output arguments:
% r1: Value of the radial wave function at 20 fm

q = 1.4399764;        % e^2 in units of MeV*fm
a = 0.8;              % Surface "thickness" of the parent nucleus, units fm
w = 1.3*(A^(1/3));    % "Radius" of the daughter nucleus, units fm
X = 1.9145*(10^(-5)); % The constant term "2m(d/h)^2" in units of 1/MeV
r1 = 0;               % Value of u[i] at the distance of 20 fm, starting value arbitrary
u0 = 0;               % u[i-1]- term
u1 = 1;               % u[i]- term
r = 0.01:0.01:20;     % Vector of iteration steps
Z = 92;               % Proton number of the parent nucleus (Uranium)

for i = 1:length(r)
    
    % First, calculate the value of the potential energy terms
    V_em1 = 2*(Z-2)*q*(3*(w^2)-(r(i)^2))/(2*(w^3)); % r<=w
    V_em2 = 2*(Z-2)*q/r(i);                         % r>w
    x = (r(i)-w)/a;
    V_ws = V0/(exp(x)+1);
    
    % Calculate u2=u[i+1] 
    if r(i)<=w
        u2 = 2*u1+X*(V_em1+V_ws-E)*u1-u0; 
        
    elseif r(i)>w
        u2 = 2*u1+X*(V_em2+V_ws-E)*u1-u0;
            
    end
    % Update variables for next iteration step
    u0 = u1;     % u[i-1] --> u[i]
    u1 = u2;     % u[i] ---> u[i+1]
  
    % If we are done, u[i+1]=u2==r1
    if r(i)==r(end)
        r1 = u2; % u[i+1] --> u(r1)
        
    end
    
end
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function [Ea, n, T] = Elevel(A, V0)
% Calculates the first positive (quasibound) energy level. All the energy
% levels are calculated by means of the bisection methotd described in the
% thesis paper. Iteration of the Schrödinger equation is done by calling 
% the function schr() defined above.
%
% Input arguments:
% A: Mass number of the DAUGHTER nucleus
% V0: A value for the Vo term in Woods-Saxon potential
%
% Output arguments:
% Ea: The first positive (quasibound) energy level
% n: Radial quantum number
% T: Vector of all energy levels

z = 0;        % Parameter for the while- loop
r0 = 0;       % Result of the previous run of schr(), starting value arbitrary 
dE = 0.5;     % Staring value for the energy increment
E = V0+1;     % Energy level in the beginning of the iteration
E1 = V0+1;    % Last bound energy level calculated, starting value arbitrary 
n = 0;        % Starting value for the radial quantum number 
T = zeros(1); % Will be a vector of energy levels

while z<1
    
    % Calculate r1 from schr()
    r1 = schr(E, A, V0);
    
    % If the last energy level calculated is negative:
    if E1<0
        % And if r0 is 0 or r0 and r1 have the same sign
        if ((r0==0) || (r0*r1>0))
            r0 = r1;  % Update r0
            E = E+dE; % Update E
        
             
        % or if ro and r1 have different signs, but desired accuracy has 
        % not been acheived    
        elseif ((r0*r1<0) && (abs(dE)>10^(-6)))
                r0 = r1;          % Update r0
                dE = (dE/2)*(-1); % Update dE
                E = E+dE;         % Update E
                
        % or if r0 and r1 have different signs and the desired accuracy has 
        % been reached
        elseif ((r0*r1<0) && (abs(dE)<10^(-6)))
                n = n+1;    % Update the radial quantum number
                T(1,n) = E; % Update the array for energy levels
                E1 = E;     % Update the value for most recent found bound state
                % Reset variables for next round
                r0 = 0;
                E = E+1;
                dE = 0.5;
                
        end
        
    end
    
    % If the most recent found energy level was positive, then we are
    % finished
    if E1>0
        Ea = T(end); % Output Ea gets the value of last element in the array T 
        z=2;         % update  the while-loop parameter in order to break the loop
        
    end
    
end
end
%--------------------------------------------------------------------------







                
        
    
                
            
















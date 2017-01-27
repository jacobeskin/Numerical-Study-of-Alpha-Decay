function [r1, r2, V, X] = radii(Ealpha, A, Vo)
% Determines the ends of the tunnel through which the alpha particle
% penetrates the potential barrier. Also gives the values of the potential
% energy from the 0.01 fm up to 100 fm and the tunneling interval. The ends
% of the tunnel satisfy the relation Ealpha=V(r1)=V(r2).
%
% Input arguments:
% Ealpha: Alpha particle energy
% A: Mass number of the DAUGHTER nucleus
% Vo: Vo term of the Woods-Saxon potential
%
% Output arguments:
% r1: Beginning of the tunnel
% r2: End of the tunnel
% V: Array of values of the potential in the tunneling interval
% X: Array of the values the potential energy gets in the interval
%    0.01->100 fm

q = 1.4399764;          % e^2 in units MeV*fm
a = 0.8;                % Surface "thickness" of the nucleus 
w = 1.3*(A^(1/3));      % Radius of the daughter nucleus
r = 0.01:0.01:100;      % Array of the discrete steps to which 100 fm are divided 
Z = 92;                 % Proton number of the parent nucleus (Uranium)
X = zeros(1,length(r)); % Array for the values of the potential energy between the interval 0.01->100 fm

% Calculate the value of the potenrial between 0.01->100 fm
for i = 1:length(r)
    
    V_em1 = 2*(Z-2)*q*(3*(w^2)-(r(i)^2))/(2*(w^3)); % r<=w
    V_em2 = 2*(Z-2)*q/r(i);                         % r>w
    x = (r(i)-w)/a;          
    V_ws = Vo/(exp(x)+1);
    
    % Calculating the value depending if r<=w or r>w    
    if r(i)<=w     
        v = V_em1+V_ws;
    elseif r(i)>w
        v = V_em2+V_ws;
    end
    
    X(1,i) = v; % Assigning the value of the potential to the array X
    
end

% Declairing these variables now, avoiding a weird error....
r1 = 0;
r2 = 0;
I1 = 0;
I2 = 0;

% Determining the values of r for which the relation Ealpha=V(r1)=V(r2) 
% holds as closely as possible
for j = 1:length(X)-1
    
    % Follow 2 consecutive values of the potenital energy
    x1 = X(j);
    x2 = X(j+1);
    
    % If Ealpha is either of these values of falls in between for V(r1)
    if ((x1<=Ealpha) && (x2>=Ealpha))
        % We save the position and value of x2, since there may be many 
        % values of x1 and x2 this might apply; we take the largest value
        r1 = r(j+1); 
        I1 = j+1;
        
    % If Ealpha is eithr x1 or x2 or falls in between for V(r2)
    elseif ((x1>=Ealpha) && (x2<=Ealpha))
        % We save the position and value of x1; reason above 
        r2 = r(j);
        I2 = j;
    end
end

% Create the array for values of the potential in the tunneling interval
V = X(1,I1:I2);

end
        
        
        
        
        
        
        
        
        
        
        
    
   
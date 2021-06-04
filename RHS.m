function dx = RHS(x,y)
% Evaluate RHS of compressible B.L. similarity solution
% Self-similar formulation
%    (Cf'')' + f f'' = 0
%    (Cg')' + Pr f g' + Pr C (gamma-1) M^2 f''^2 = 0
% where f'=u/ue, g=T/Te and C=(rho*mu/(rhoe*mue)=g^(n-1)
% inputs:
% t - variable required by ode45
% x - state vector = (f,f',f'',g,g')
% params - structure with 4 parameters described below
%-----------------------------------------------------

% n = params.n; % power law coefficient for viscosity μ/μ0=(T/T0)^n
% gamma = params.gamma; % specific heat ratio
% Pr = params.Pr; % Pr - Prandtl number
% M = params.M; % Mach

dx = zeros(5,1);
dx(1) = y(2); 
dx(2) = y(3); 
dx(3) = y(2)^2-2*y(1)*y(3)-(y(4)+1)^2;
dx(4) = y(5); 
dx(5) = 2*(y(4)+1)*y(2)-2*y(1)*y(5);%V''
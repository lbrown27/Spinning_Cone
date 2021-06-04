function [eta,f]=solve_comp(M,flag_ad,Tw,Pr,gamma,n,etamax)
% Solve compressible B.L. similarity equations
%    (Cf'')' + f f'' = 0
%    (Cg')' + Pr f g' + Pr C (gamma-1) M^2 f''^2 = 0
% where f'=u/ue, g=T/Te and C=(rho*mu)/(rhoe*mue)=g^(n-1)
% the iterative solution is using Broydens method
% inputs:
% flag_ad - 1=adiabatic wall (Tw not used); 0=const. temperature (Tw is used)
% M - Mach
% Pr  Prandtl number
% gamma - adiabatic coefficient
% n - power law coefficient for viscosity mu/mue=(T/Te)^n
% Tw - wall temperature (not used in adiabatic case)
% etamax - max. value for eta (10 should work fine)
% outputs:
% eta - similarity variable
% f - state vector = (f,f',f'',g,g')
%-----------------------------------------------------
% initiate params:
params.M = M;
params.Pr = Pr;
params.gamma = gamma;
params.n = n;

ic = @(q) [0 0 q(1) 0 q(2)];

% initial guess for f''(0) and g(0)
q0 = [6.16;-2.65];
% solve for current guess:
[eta,f]=ode15s(@(t,y) RHS(t,y,params),[0 etamax],ic(q0));
% compare values to required freestream values:
F0 = [f(end,1); f(end,4)+1];
% Guess for jacobian:



%%%% NEW AND IN PROGRESS
dx = .01;
q1 = [q0(1)+dx,q0(2)];
[eta,f1]=ode15s(@(t,y) RHS(t,y,params),[0 etamax],ic(q1));
F1 = [f1(end,1); f1(end,4)+1];

q2 = [q0(1),q0(2)+ dx];
[eta,f2]=ode15s(@(t,y) RHS(t,y,params),[0 etamax],ic(q2));
F2 = [f2(end,1); f2(end,4)+1];

J_val(1) = (F1(1)-F0(1))/dx;
J_val(2) = (F2(1)-F0(1))/dx;
J_val(3) = (F1(2)-F0(2))/dx;
J_val(4) = (F2(2)-F0(2))/dx;
J = [J_val(1),J_val(2);J_val(3),J_val(4)];



%%%




J = eye(2) ;

k = 1;
stop = 50;
while sum(abs(F0))>1e-6 && k<stop
    % calc next guess:
    d = -J\F0;
    q = real(q0+d);
    
%     q(1) = min(q(1),1e8);
%     q(1) = max(q(1),-1e8);
    % solve for current guess:
    fprintf('\nnew guess: \n');
    q
    [eta,f]=ode15s(@(t,y) RHS(t,y,params),[0 etamax],ic(q));
    % compare values to required freestream values:
    F = [f(end,1); f(end,4)+1];
    % update Jacobian
    J = J + ((F-F0)-J*d)*d'/(d'*d);
    q0 = q;
    F0 = F;
    k = k + 1;
    k
end

if k==stop
    error('Failed to converge, try changing initial guess')
else
    fprintf('converged in %d iterations\n',k);
end
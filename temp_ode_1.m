function dx = temp_ode_1(x,y,psi_t, psi_double_prime,V_prime,psi,Pr)
global streamfn Pr
if ~exist('streamfn', 'var')

[streamfn,~,~] = streamfunction();
end


psi_t = linspace(0,100,length(streamfn.x));
psi_double_prime = streamfn.y(3,:);
V_prime = streamfn.y(5,:);
psi = streamfn.y(1,:);

psi_double_prime = interp1(psi_t,psi_double_prime,x); % Interpolate the data set (ft,f) at time t
V_prime = interp1(psi_t,V_prime,x); % Interpolate the data set (gt,g) at time t
psi = interp1(psi_t,psi,x); % Interpolate the data set (gt,g) at time t

dx = zeros(2,1);
dx(1) = y(2);
dx(2) = 2*Pr*(psi_double_prime^2+V_prime^2-psi*y(2)+y(1));
end
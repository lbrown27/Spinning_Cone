function dx = temp_ode_2(x,y)
global streamfn Pr
if ~exist('streamfn', 'var')

[streamfn,~,~] = streamfunction();
end


psi_t = linspace(0,100,length(streamfn.x));
psi = streamfn.y(1,:);

psi = interp1(psi_t,psi,x);

dx = zeros(2,1);
dx(1) = y(2);
dx(2) = -2*Pr*psi*y(2);
end
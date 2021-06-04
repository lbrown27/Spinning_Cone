function [T,xi,n] = temp_field()

global streamfn Pr Ma gamma Tw
Pr = 1;
if isempty(streamfn)
    
    [streamfn,~,~] = streamfunction();
end
MAX = 10;
xmesh = linspace(0,MAX,400);
solinit = bvpinit(xmesh,[0 -20]);


sol1 = bvp4c(@temp_ode_1, @bcfcn1, solinit);
%  figure();
%  plot(sol1.x,sol1.y(1,:));
%  title('f function');



solinit = bvpinit(xmesh,[1 0]);


sol2 = bvp4c(@temp_ode_2, @bcfcn2, solinit);
% figure();
% plot(sol2.x,sol2.y(1,:));
% title('q function')



%% find temp distribution


n = min(length(sol1.y(1,:)),length(sol2.y(1,:)));
xi = linspace(0,MAX, n);
f = interp1(sol1.x,sol1.y(1,:),xi);
q = interp1(sol2.x,sol2.y(1,:),xi);
for i = 1:n
    T(i) = 1- (gamma - 1)/2 * Ma^2*f(i)+(Tw - 1)* q(i);
end



    function res = bcfcn1(ya,yb)
        res = [ya(1)
            yb(1)];
    end

    function res = bcfcn2(ya,yb)
        res = [ya(1)-1
            yb(1)];
    end
end
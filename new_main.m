close all;
clc;

global streamfn Pr Ma gamma Tw

Pr = .7;
Ma = 8;
gamma = 1.4;
Tw = 1.5;

%[T,xi,n] = temp_field();
figure();
plot(xi,T);
hold on;
ylim([0 2]);
title('Temperature Distribution, Pr = .7, Ma = 8, Tw = 1.5');
xlabel('y');
[~,ind] = min(abs(streamfn.x(:) - 10));

figure();
plot(streamfn.x(1:ind),streamfn.y(1,1:ind));
hold on;
plot(streamfn.x(1:ind),streamfn.y(2,1:ind));

plot(streamfn.x(1:ind),-1 * streamfn.y(4,1:ind));

legend('psi','U','-V');
xlabel('y');

for Ma = [1,2,4,8]
    figure();
   for Tw =  0:.4:2
        [T,xi,n] = temp_field();
        plot(xi,T);
        hold on;
       
   end
    title("Mach Number: " + Ma);
    legend('Tw = 0', 'Tw = .4','Tw = .8','Tw = 1.2', 'Tw = 1.6', 'Tw = 2');
        ylabel('Temp');
        xlabel('y');

end

for Tw = .5:1.5
    figure();
   for Ma = 0:2:10
        [T,xi,n] = temp_field();
        plot(xi,T);
        hold on;
        xlabel('y')

       
   end
    title("Wall Temp: " + Tw);
    legend('Ma = 0', 'Ma = 2','Ma = 4','Ma = 6', 'Ma = 8', 'Ma = 10')
    ylabel('Temp');

end
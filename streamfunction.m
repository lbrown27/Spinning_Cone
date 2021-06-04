function [smallest_sol, guess1_best, guess2_best] = streamfunction()

xmesh = linspace(0,450,500);
smallest_val = 10000000;
bests1 = 0;
bests2 = 0;
for guess1 = -.3:.03:0
    for guess2 = -.5:.03:.0
        
        q0 = [guess1;guess2]
        bc = @(q) [0 0 q(1) 0 q(2)];
        solinit = bvpinit(xmesh, bc(q0));
        options = bvpset('RelTol',1e-8);
        try
            sol = bvp5c(@RHS, @bcfcn, solinit,options);
            my_sum = sum(abs(sol.y(1,:)));
            k = abs(sol.y(1,end)^2+(sol.y(4,end)+1)^2+sol.y(5,end)^2 + sol.y(2,end)^2);
            if k < smallest_val
                guess1_best = guess1;
                guess2_best = guess2;
                smallest_val = k;
                smallest_sol = sol;
                fprintf('\nsmallest val is now %.4f at %.3f , %.3f\n', smallest_val, guess1,guess2);
                plot(sol.x, sol.y(1,:));
                smallest_sum = my_sum;
                title("guess 1: " + guess1 + " guess 2: " + guess2);
                drawnow;
                bests1 = [bests1,guess1];
                bests2 = [bests2,guess2];
            end
            %              plot(sol.x, sol.y(1,:));
            %                 title("guess 1: " + guess1 + " guess 2: " + guess2);
            %                 drawnow;
        catch
        end
    end
end

F0 = [smallest_sol.y(1,end); smallest_sol.y(4,end)+1];

J = eye(2);
k = 1;
stop = 50;
while sum(abs(F0))>1e-8 && k<stop
    % calc next guess:
    d = -J\F0;
    q = real(q0+d);
    
    %     q(1) = min(q(1),1e8);
    %     q(1) = max(q(1),-1e8);
    % solve for current guess:
    fprintf('\nnew guess: \n');
    q
    solinit = bvpinit(xmesh, bc(q));
    sol = bvp5c(@RHS, @bcfcn, solinit,options);
    % compare values to required freestream values:
    F = [sol.y(1,end); sol.y(4,end)];
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


    function res = bcfcn(ya,yb)
        res = [ya(1) yb(1) ya(2) ya(4) yb(4)+1];
    end

    function g = guess()
        g = [6 -2.5];
    end
end

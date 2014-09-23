function [x, u] = exc4f(dt, t_start, t_end, dx, x_start, x_end, initial, f, s, method, method2, limiter)

x = x_start:dx:x_end;
t = t_start:dt:t_end;

u(1,:) = initial(x); % Add the initial condition as the first result

for i=1:numel(t) %Loop through time and get their results
    u(i+1,:) = method(dt, dx, u(i,:), f, s , method2, limiter);
    g = figure(3);
    plot(x,u(i+1,:));
    title(strcat('t = ', num2str(t(i))));
     
    for ix=1:size(u,2)
        if(abs(dt/dx * f(u(i+1,ix))) > 1) % Apply the CFL-condition to check if the used dt and dx give us stable results
            error('Solution is invalid, timestep restriction is not OK!'); % Stop the script and give an error
        end
    end
    
end

end
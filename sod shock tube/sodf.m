function [x, u, rho, v, p] = hoofdopdrachtf(dt, t, dx, x, u, gamma)

[rho(:,1), v(:,1), p(:,1), H(:,1)] = q(u(:,:,1), gamma); % Get the values for rho, v, p and H

for i=2:numel(t) %Loop through time and get their results
    
    k1 = dt * muscl(u(:,:,i-1), dx, gamma);

    ustar = u(:,:,i-1) - 0.5 * k1;
    k2 = dt * muscl(ustar, dx, gamma);

    ustar = u(:,:,i-1) - 0.5 * k2;
    k3 = dt * muscl(ustar, dx, gamma);

    ustar = u(:,:,i-1) - k3;
    k4 = dt * muscl(ustar, dx, gamma);

    u(:,:,i) = u(:,:,i-1) - 1/6 * (k1 + 2*k2 + 2*k3 + k4); %Runge-Kutta
    
    [rho(:,i), v(:,i), p(:,i), H(:,i)] = q(u(:,:,i), gamma);  % Get the values for rho, v, p and H
    
    h = figure(1);
    
    subplot(1,3,1); plot(x,rho(:,i)); xlabel('x'); ylabel('\rho'); axis([x(1) x(end) 0 2]);
    subplot(1,3,2); plot(x,v(:,i)); xlabel('x'); ylabel('u'); axis([x(1) x(end) 0 2]);
    subplot(1,3,3); plot(x,p(:,i)); xlabel('x'); ylabel('p'); axis([x(1) x(end) 0 2]);
    title(strcat('t = ', num2str(t(i))));
    drawnow
    
end

close(h);

end
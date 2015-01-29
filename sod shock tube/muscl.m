function [ux] = muscl(u, dx, gamma)

nx = size(u,2);

ux = u;

for i = 2:nx-1
    if(i == 2 || i == (nx-1))
        ulm05 = u(:,i-1);
        urm05 = u(:,i);
        
        ulp05 = u(:,i);
        urp05 = u(:,i+1);
    else
        ulm05 = ustarL(u, i-1); % uL_{i+1/2}
        urm05 = ustarR(u, i); % uR_{i-1/2}

        ulp05 = ustarL(u, i); % uL_{i+1/2}
        urp05 = ustarR(u, i+1); % uR_{i-1/2}
    end
    
    ux(:,i) = (froe(ulp05, urp05, gamma)-froe(ulm05, urm05, gamma))/dx;
end

ux(:,1) = ux(:,2); ux(:,nx) = ux(:,nx-1); % Neumann boundary conditions

end

function f = ustarL(u, i)

r = ratio(u, i);

%f = u(i) + 1/2*r*limiter(1/r)*(u(i+1) - u(i)); % Equation (6.16) from reader
f = u(:,i) + 1/2 * vanalbada(r) .* (u(:,i+1) - u(:,i)); % Equation (6.15) from reader

end

function f = ustarR(u, i)

r = ratio(u, i);

%f = u(i) - 1/2*r*limiter(1/r) * (u(i+1) - u(i)); % Equation (6.16) from reader
f = u(:,i) - 1/2 * vanalbada(r) .* (u(:,i+1) - u(:,i)); % Equation (6.15) from reader

end

function r = ratio(u, i)
r = (u(:,i)-u(:,i-1))./(u(:,i+1)-u(:,i));
end
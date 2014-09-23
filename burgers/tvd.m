function ux = tvd(u, i, dx, f, s, limiter)

% uiL1 = ustarL(u, i, limiter); % uL_{i+1/2}
% uiR1 = ustarR(u, i+1, limiter); % uR_{i-1/2}
% 
% fip12 = f(  riemann(uiL1, uiR1, s) ); % f_{i+1/2}
% 
% uiL2 = ustarL(u, i-1, limiter); % uL_{i+1/2}
% uiR2 = ustarR(u, i, limiter); % uR_{i-1/2}
% 
% fim12 = f( riemann(uiL2, uiR2, s) ); % f_{i-1/2}
% 
% ux = ( fip12 - fim12 ) / dx; % (f_{i+1/2} - f_{i-1/2})/dx

uiLm32 = ustarL(u, i-1, limiter); % uL_{i+1/2}
uiRm12 = ustarR(u, i, limiter); % uR_{i-1/2}

fim1 = f(  riemann(uiLm32, uiRm12, s) ); % f_{i+1/2}

uiLm12 = ustarL(u, i, limiter); % uL_{i+1/2}
uiRp12 = ustarR(u, i+1, limiter); % uR_{i-1/2}

fi = f( riemann(uiLm12, uiRp12, s) ); % f_{i-1/2}

ux = ( fi - fim1 ) / dx; % (f_{i+1/2} - f_{i-1/2})/dx

end

function f = ustarL(u, i, limiter)

r = ratio(u, i);

%f = u(i) + 1/2*r*limiter(1/r)*(u(i+1) - u(i)); % Equation (6.16) from reader
f = u(i) + 1/2 * limiter(r) * (u(i+1) - u(i)); % Equation (6.15) from reader

end

function f = ustarR(u, i, limiter)

r = ratio(u, i);

%f = u(i) - 1/2*r*limiter(1/r) * (u(i+1) - u(i)); % Equation (6.16) from reader
f = u(i) - 1/2 * limiter(r) * (u(i+1) - u(i)); % Equation (6.15) from reader

end

function r = ratio(u, i)
r = (u(i)-u(i-1))/(u(i+1)-u(i));
end
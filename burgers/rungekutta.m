function [q] = rungekutta(dt, dx, q, f, s, method, limiter)

k1 = dt * method(q, dx, f, s, limiter);

ustar = q - 0.5 * k1;
k2 = dt * method(ustar, dx, f, s, limiter);

ustar = q - 0.5 * k2;
k3 = dt * method(ustar, dx, f, s, limiter);

ustar = q - k3;
k4 = dt * method(ustar, dx, f, s, limiter);

q = q - 1/6 * (k1 + 2*k2 + 2*k3 + k4);

end
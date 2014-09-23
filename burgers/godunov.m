function ux = godunov(u, i, dx, f, s)

ux = (f(riemann(u(i), u(i+1), s)) - f(riemann(u(i-1), u(i), s)))/dx;

end
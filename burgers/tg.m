function ux = tg(u, dx, f, s, limiter)

ux = u;

for i=[2, numel(u)-1] % Use Godunov next to the boundaries
    ux(i) = godunov(u, i, dx, f, s);
end

for i=3:numel(u)-2 % Use TVD for all the other points in the grid
    ux(i) = tvd(u, i, dx, f, s, limiter);
end

ux(1) = 0; % We don't have to calculate the new values on the boundaries
ux(numel(u)) = 0; % We don't have to calculate the new values on the boundaries

end
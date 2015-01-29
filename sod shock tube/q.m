function [rho, u, p, H] = q(U, gamma)

rho = U(1,:,:);
u = U(2,:,:)./U(1,:,:);
E = U(3,:,:)./U(1,:,:);
p = (gamma-1)*(rho.*E-1/2*rho.*u.^2);
H = E + p./rho;

end
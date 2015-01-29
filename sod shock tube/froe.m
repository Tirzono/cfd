function flux = froe(ul, ur, gamma)

global F;
global R;
global lambda;
global Z_bar;

[rhoL, vL, pL, HL] = q(ul, gamma);
[rhoR, vR, pR, HR] = q(ur, gamma);

Z = Z_bar(rhoL,vL,HL,rhoR,vR,HR,gamma)';
R_bar = R(Z);
lambda_bar = lambda(Z);

flux = 1/2*( F(rhoR, vR, pR, HR) + F(rhoL, vL, pL, HL) - (R_bar*abs(diag(lambda_bar))*inv(R_bar))*(ur-ul) );

%F_{i+1/2} = 1/2(F_{i+1}+F_{i}) - 1/2(R*lambda*R^(-1))

end
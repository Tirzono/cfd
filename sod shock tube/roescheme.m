clear all; clc; close all;

%% Constants
dt = 0.001;
dx = 0.01;

t_start = 0;
t_end = 0.15;

x_start = -1/2;
x_end = 3/2;

snapshots = 0:0.025:t_end; % At what time steps do we have to plot a result

gamma = 1.4;

split = 0; %The place where the split is between u_L and u_R

rho_L = 1;
p_L = 1;
u_L = 0;

rho_R = 0.125;
p_R = 0.1;
u_R = 0;

%% Domain
x = linspace(-1/2,3/2,dx);          % cells nodes
nx = numel(x);

%% Time 
t  = t_start:dt:t_end;
nt = numel(t);

%% Initial condition

for i = 1:numel(x)
    if(x(i) < split)
        u(1,i,1) = rho_L; 
        u(2,i,1) = rho_L*u_L;
        u(3,i,1) = p_L/(gamma - 1)+rho_L*(1/2*u_L^2);
    else
        u(1,i,1) = rho_R;
        u(2,i,1) = rho_R*u_R;
        u(3,i,1) = p_R/(gamma - 1)+rho_R*(1/2*u_R^2);
    end
end

L = 1:nx-1; R = 2:nx; % inner nodes / middle points

for n = 2:nt; 
    % Compute Roe Averages
    % Velovity 'u'
    u_bar = (sqrt (u(1,L,n-1)).*u(L) + sqrt (r(R)).*u(R)) ...
                ./ (sqrt(r(L))+sqrt(r(R)));
    % Total Entalpy 'H'
    H_bar = (sqrt (r(L)).*H(L) + sqrt (r(R)).*H(R)) ...
                ./ (sqrt(r(L))+sqrt(r(R)));
    % Sound Speed 'a'
    a_bar = sqrt((gamma-1)*(H_bar-0.5*u_bar.^2));

    % Compute Delta U's
    % dU = U(R)-U(L) at the cell boundaries { x_{i},x_{i+1}, ... }
    dU = U(:,R)-U(:,L);

    % Compute Fluxes at the cell centers
    % F = [F1 F2 F3]
    F = [r.*u; r.*u.^2 + p; r.*u.*H]; 
    % Fluxes to the left and right of 'F_ {i+1/2}'
    FL = F(:,L); FR = F(:,R);
    
    % Scaled Right Eigenvectors are given by:
    k1_bar = [1*ones(1,nx-1); u_bar-a_bar; H_bar-u_bar.*a_bar];
    k2_bar = [1*ones(1,nx-1); u_bar      ; 1/2*u_bar.^2      ];
    k3_bar = [1*ones(1,nx-1); u_bar+a_bar; H_bar+u_bar.*a_bar];
    
    % compute Roe waves strength alpha_bar
    alpha2_bar = (gamma-1)./(a_bar.^2).*(dU (1,:).*(H_bar-u_bar.^2) ...
        + u_bar.*dU(2,:) - dU(3,:));
    alpha1_bar = 1./(2*a_bar).*(dU (1,:).*(u_bar+a_bar) ...
        - dU(2,:)-a_bar.*alpha2_bar);
    alpha3_bar = dU(1,:)-(alpha1_bar + alpha2_bar);
    
    % Eigenvalues of A_bar are (same as the original A mat)
    lambda_bar(1,:) = abs(u_bar - a_bar);
    lambda_bar(2,:) = abs(u_bar);
    lambda_bar(3,:) = abs(u_bar + a_bar);
    
%     % Entropy Fix
%     boolv1 = lambda_bar(1,:) < etpfix;
%         lambda_bar(1,:) = 0.5*(etpfix + lambda_bar (1,:).^2/etpfix).*boolv1 ...
%             + lambda_bar(1,:).*(1-boolv1);
%     boolv2 = lambda_bar(3,:) < etpfix;
%         lambda_bar(3,:) = 0.5*(etpfix + lambda_bar (3,:).^2/etpfix).*boolv2 ...
%             + lambda_bar(3,:).*(1-boolv2);
        
    % Conditioning data
    alpha1_bar = repmat(alpha1_bar,3,1);
    alpha2_bar = repmat(alpha2_bar,3,1);
    alpha3_bar = repmat(alpha3_bar,3,1);
    lambda1_bar = repmat(lambda_bar(1,:),3,1);
    lambda2_bar = repmat(lambda_bar(2,:),3,1);
    lambda3_bar = repmat(lambda_bar(3,:),3,1);
    
    % Roe Fluxes
    Flux = 0.5*(FL+FR)-0.5*(alpha1_bar.*lambda1_bar.*k1_bar + ...
                            alpha2_bar.*lambda2_bar.*k2_bar + ...
                            alpha3_bar.*lambda3_bar.*k3_bar );
                        
    % Compute next time step
    
    for i = 2:nx-1
% if (x(1,i) < 0.335)
%      BGK-CODE;
%    elseif( x(1,i) > 0.675)
%        Roe-Euler
     U_next(:,i) = U(:,i) - dtdx*(1-h(1,i))*(Flux(:,i) - Flux(:,i-1));
%    else
%    Roe-Euler &  BGK Coupling      
%    h(1,i)= (x(1,i)-0.675)/(0.335-0.675);
   end
%end
    
    
   % BCs
    U_next(:,1)  = U(:,2); % Neumann condition to the left
    U_next(:,nx) = U(:,nx-1); % Neumann condition to the right
    
    % Compute variables of the new time step
    r_next = U_next(1,:);               % Density
    u_next = U_next (2,:)./U_next(1,:);  % Velocity
    E_next = U_next(3,:);               % Total Energy
    p_next = (gamma-1).*(E_next-r_next.*u_next.^2/2);  % Pressure
    e_next = 1/(gamma-1)*(p_next./r_next);      % Internal Energy
    a_bar_next = sqrt(gamma*p_next./r_next);    % sound speed
    m_next = u_next./a_bar_next;        % Mach 
    s_next = log(p_next./r_next.^gamma);% Entropy
    H_next = (E_next + p_next)./r_next; % Enthalpy
    
    % Update info
    U = U_next;
    r = r_next;
    u = u_next;
    e = e_next;
    p = p_next;
    m = m_next;
    s = s_next;
    H = H_next;
        
    
    % Plot figure
    if plot_fig == 1;
        subplot(2,3,1); plot(x,r); title('Velocity');
        subplot(2,3,2); plot(x,u); title('ernal Energy');
        subplot(2,3,3); plot(x,p); title('Pressure');
        subplot(2,3,4); plot(x,m); title('Mach number');
        subplot(2,3,5); plot(x,s); title('Entropy');
        subplot(2,3,6); plot(x,e); title('IntDensity');
    end
    drawnow
end

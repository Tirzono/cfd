%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Numerical solution to Sod Shock Tube problem with use of Roe's
% Approximate Riemann Solver
%
% Boundary conditions:
% pL = 1, rhoL = 1, uL = 0
% pR = 0.125, rhoR = 0.1, uR = 0
%
% Time-discretization: fourth order Runge-Kutta scheme
% Space-discretization: second order TVD-scheme with MC-limiter (Godunov
% for points near boundaries)
%
% Limiters possible: 
% @mc - MC-limiter
% @vanalbada - Van Albada limiter
% @minmod - MinMod limiter
% @superbee - SuperBee limiter
% @vanleer - Van Leer limiter
%
% Script by Sebastiaan ten Pas and Martin Goossens
% 29-1-2015 | University of Twente
%
% For more information: info@sebastiaantenpas.nl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Start script

clear all;
close all;
clc;
format long e;

%% Constants

dx_start = 0.02;
dxn = 5; 

dx = dx_start./2.^((1:dxn)-1); % Spatial step
dt = dx/2; % Time step

t_start = 0;
t_end = 0.15;

x_start = -1/2;
x_end = 3/2;

gamma = 1.4;

x0 = 0.5; %The place where the split is between u_L and u_R

%% Initial conditions

rho_L = 1;
p_L = 1;
u_L = 0;

rho_R = 0.125;
p_R = 0.1;
u_R = 0;

%% Determine the functions

global F;
global R;
global lambda;
global Z_bar;

F = @(rho, u, p, H) [rho*u; rho*u^2 + p; rho*u*H];
U = @(rho, u, p, gamma) [rho; rho*u; p/(gamma - 1)+rho*(1/2*u^2)];
R = @(Z) [1, 1, 1; % 1, 1, 1
    Z(1)-Z(3), Z(1), Z(1)+Z(3); %u - a, u, u+a
    Z(2)-Z(1)*Z(3), 1/2*Z(1)^2, Z(2)+Z(1)*Z(3)]; %H-u*a, 1/2*u^2, H+u*a
lambda = @(Z) [Z(1) - Z(3), Z(1), Z(1) + Z(3)]; %u - a, u, u + a
Z_bar = @(rhoL,vL,HL,rhoR,vR,HR,gamma) [(sqrt (rhoL) * vL + sqrt (rhoR) * vR) / (sqrt(rhoL)+sqrt(rhoR));
(sqrt (rhoL) * HL + sqrt (rhoR) * HR) / (sqrt(rhoL)+sqrt(rhoR));
sqrt((gamma-1)*(((sqrt (rhoL) * HL + sqrt (rhoR) * HR) / (sqrt(rhoL)+sqrt(rhoR)))- 1/2*((sqrt (rhoL) * vL + sqrt (rhoR) * vR) / (sqrt(rhoL)+sqrt(rhoR)))^2))];

%% Make a new folder to save the data to

tp = mfilename('fullpath'); %Look up the path
dir = {tp(1:(end-3))};      %Strip the name of the script from the path
tclock = fix(clock);        %Get the current date and time
nmap = strcat(date, '-', num2str(tclock(4)), num2str(tclock(5)),num2str(tclock(6))); %Name of the new map
path = strcat(dir{1}, 'data-files/'); %Path of the new map
mkdir(path, nmap);          %Make a new directory to place the results in

disp(sprintf('dx\t\t\t\tdt\t\t\t\tL1 rho\t\t\tr rho\t\t\tL1 u\t\t\tr u\t\t\t\tL1 p\t\t\tr p')); % Print the table header

for j = 1:numel(dx)
    x = x_start:dx(j):x_end;
    t = t_start:dt(j):t_end;

    for i = 1:numel(x)
        if(x(i) < x0)
            u(:,i,1) = U(rho_L, u_L, p_L, gamma); % Calculate the intial solution before the interface
        else
            u(:,i,1) = U(rho_R, u_R, p_R, gamma); % Calculate the intial solution after the interface
        end
    end

    %% Call the function for the approximated solution

    [x, u, rho, v, p] = sodf(dt(j), t, dx(j), x, u, gamma);
    
    E = p./(rho * (gamma - 1)) + 1/2 * v.^2;

    %% Exact solution

    [x1, x2, x3, x4, rhoe, ve, pe] = exactSOD(rho_L, u_L, p_L, rho_R, u_R, p_R, gamma, x, x0, t_end);

    Ee = pe./(rhoe * (gamma - 1)) + 1/2 * ve.^2;

    %% Plot the function
    
    c = hsv(2);

    fi(j, 1) = figure(4*j-2);
    plot(x, rho(:,end), 'o', 'LineWidth', 2, 'Color', c(1,:)); xlabel('x'); ylabel('\rho');  axis([x(1) x(end) 0 2]); title(strcat('Solution for \rho for dx = ', num2str(dx(j)),' and dt = ', num2str(dt(j))));
    saveas(fi(j, 1), strcat(path, nmap, '/excH-rho-dx-', num2str(dx(j)) ,'.eps'), 'epsc');
    hold on; plot(x, rhoe, 'LineWidth', 2, 'Color', c(2,:)); legend(strcat('t = ', num2str(t_end), ' s'), strcat('Exact at t = ', num2str(t_end), ' s'));
    saveas(fi(j, 1), strcat(path, nmap, '/excH-rho-exact-dx-', num2str(dx(j)) ,'.eps'), 'epsc'); 
    
    fi(j, 2) = figure(4*j-1);
    plot(x, v(:,end), 'o', 'LineWidth', 2, 'Color', c(1,:)); xlabel('x'); ylabel('u');  axis([x(1) x(end) 0 2]); title(strcat('Solution for u for dx = ', num2str(dx(j)),' and dt = ', num2str(dt(j))));
    saveas(fi(j, 2), strcat(path, nmap, '/excH-v-dx-', num2str(dx(j)) ,'.eps'), 'epsc'); 
    hold on; plot(x, ve, 'LineWidth', 2, 'Color', c(2,:)); legend(strcat('t = ', num2str(t_end), ' s'), strcat('Exact at t = ', num2str(t_end), ' s'));
    saveas(fi(j, 2), strcat(path, nmap, '/excH-v-exact-dx-', num2str(dx(j)) ,'.eps'), 'epsc'); 
    
    fi(j, 3) = figure(4*j);
    plot(x, p(:,end), 'o', 'LineWidth', 2, 'Color', c(1,:)); xlabel('x'); ylabel('p');  axis([x(1) x(end) 0 2]); title(strcat('Solution for p for dx = ', num2str(dx(j)),' and dt = ', num2str(dt(j))));
    saveas(fi(j, 3), strcat(path, nmap, '/excH-p-dx-', num2str(dx(j)) ,'.eps'), 'epsc'); 
    hold on; plot(x, pe, 'LineWidth', 2, 'Color', c(2,:)); legend(strcat('t = ', num2str(t_end), ' s'), strcat('Exact at t = ', num2str(t_end), ' s'));
    saveas(fi(j, 3), strcat(path, nmap, '/excH-p-exact-dx-', num2str(dx(j)) ,'.eps'), 'epsc'); 
    
    fi(j, 4) = figure(4*j+1);
    plot(x, E(:,end), 'o', 'LineWidth', 2, 'Color', c(1,:)); xlabel('x'); ylabel('E');  axis([x(1) x(end) 1.5 3.5]); title(strcat('Solution for E for dx = ', num2str(dx(j)),' and dt = ', num2str(dt(j))));
    saveas(fi(j, 4), strcat(path, nmap, '/excH-E-dx-', num2str(dx(j)) ,'.eps'), 'epsc'); 
    hold on; plot(x, Ee, 'LineWidth', 2, 'Color', c(2,:)); legend(strcat('t = ', num2str(t_end), ' s'), strcat('Exact at t = ', num2str(t_end), ' s'));
    saveas(fi(j, 4), strcat(path, nmap, '/excH-E-exact-dx-', num2str(dx(j)) ,'.eps'), 'epsc'); 
    
    d = abs([rho(:,end)'; v(:,end)'; p(:,end)']-[rhoe; ve; pe]); % Calculate the difference between the exact and approximate solution
    
    for l=1:size(d,1)
        vL1(l,j) = 1/size(d,2)*sum(d(l,:)); % Calculate the L1-norm
        vL2(l,j) = sqrt(1/size(d,2)*sum(d(l,:))^2); % Calculate the L2-norm
        
        if(j > 1)
            L1(l,j-1) = vL1(l,j-1)/vL1(l,j); % Calculate the error reduction for the L1-norm
            L2(l,j-1) = vL2(l,j-1)/vL2(l,j); % Calculate the error reduction for the L2-norm
        end
    end
    
    %% Print the results (for the L1-norm) to the screen
    
    if(j > 1)
        disp(sprintf('%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f', dx(j), dt(j), vL1(1,j), L1(1,j-1), vL1(2,j), L1(2,j-1), vL1(3,j), L1(3,j-1)));
    else
        disp(sprintf('%f\t\t%f\t\t%f\t\t-\t\t\t\t%f\t\t-\t\t\t\t%f\t\t-', dx(j), dt(j), vL1(1,j), vL1(2,j), vL1(3,j)));
    end    
    
end

%% Print the results (for the L2-norm) to the screen

disp(sprintf('dx\t\t\t\tdt\t\t\t\tL2 rho\t\t\tr rho\t\t\tL2 u\t\t\tr u\t\t\t\tL2 p\t\t\tr p')); % Print the table header

for j = 1:numel(dx)
    if(j > 1)
        disp(sprintf('%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f', dx(j), dt(j), vL2(1,j), L2(1,j-1), vL2(2,j), L2(2,j-1), vL2(3,j), L2(3,j-1)));
    else
        disp(sprintf('%f\t\t%f\t\t%f\t\t-\t\t\t\t%f\t\t-\t\t\t\t%f\t\t-', dx(j), dt(j), vL2(1,j), vL2(2,j), vL2(3,j)));
    end  
end
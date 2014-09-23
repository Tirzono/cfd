%% Start script

clear all;
close all;
clc;
format long e;

%% Constants
dt = 0.0025;
dx = 0.01;

t_start = 0;
t_end = 0.8;

x_start = 0;
x_end = 1;

snapshots = 0:0.1:t_end; % At what time steps do we have to plot a result

%% Determine the function

initial = @(x) 3/2 - 2*x; % Initial condition
f = @(u) 0.5*u^2; %du/dt + df/dx = 0
s = @(ul, ur) 0.5*(ul+ur); % Shock speed

%% Call the function

[x, u] = burgersf(dt, t_start, t_end, dx, x_start, x_end, initial, f, s, @rungekutta, @tg, @mc);

%% Plot the function

fi = figure;
ns = numel(snapshots);
c = hsv(ns);

for s=1:ns
    plot(x, u((round(snapshots(s)/dt)+1),:), 'LineWidth', 2, 'Color', c(s,:)); hold on;
    ulegend{s} = strcat('t = ', num2str(snapshots(s)), ' s');
end

legend(ulegend{:}, 'Location','NorthEastOutside');
xlabel('x');
ylabel('u');
title(strcat('Solutions for u at different t for dt = ', num2str(dt), ' and dx = ', num2str(dx)));
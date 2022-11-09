% This code is the solution to Assignment 4 from course AESM1511. 
%
% Authors (student numbers): 
%  Sjaak van Meulebrouck (5845424), 
%  Carleyn den Ouden (5168457),
%  Kilian Glatz (),
%  Malte Schade ()
% 
%

%% Initial settings and parameters 
% First, start with clean system:
clear variables; % free up workspace, clear system memory
close all; % close any open figures
clc; % clear the command window

%% Defining parameters

% free space values
c0 = 299792458 ;        % [m/s] velocity
u0 = 4*pi*10^(-7);      % [H/m] magnetic permeability
eps_0 = 1/(u0*c0^2) ;   % [F/m] electric permitivity

t_w = 12*10^(-9);                                           %time window in [s]
f_c = [0.5*10^(12), 1*10^12, 2*10^12];                      % center frequencies in [Hz] 
d_r = [0.24, 0.1, 0.003, 0.05, 0.15, 0];                    % layer thickness in [m]
eps_r = [1, 6, 2, 16, 6, 9];                                % relative permitivity
u_r = [1, 1, 1, 1, 1, 1];                                   % relative magnetic permeability in [H/m]
sigma = [0, 10^(-3), 10^(-4), 10^(-2), 10^(-3), 10^(-3)];   % electric conductivity in [S/m]

%% Basic calculations
Z0 = sqrt(u0/eps_0);        
tau = sqrt(2)./f_c ;
fmax = 4*f_c ;              %maxiumum freq to consider

eps_i = eps_0.*eps_r ;      % electric permitivit per layer
u_i = u0.*u_r ;             % magnetic permeability per layer
d2 = sum(d_r(2:end-1));     % thickness of D2

%% Task 1

% compute the maximum velocity in the grid
cmax_list = sqrt(1./(eps_i.*u_i));  % calculate the velocity per layer
cmax = max(cmax_list(2:end-1));     % take the maxium velocity of D2
disp(cmax)

%compute the maximum step size, z
z = min(d_r(2:end-1))/2; 
disp(z)
delta_t = z/cmax;
disp(delta_t)


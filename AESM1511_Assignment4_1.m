% This code is the solution to Assignment 4.1 from course AESM1511. 
%
% Authors (student numbers): 
%  Sjaak van Meulebrouck (5845424), 
%  Carleyn den Ouden (5168457),
%  Kilian Glatz (),
%  Malte Schade (5850282)

%% Initial settings and parameters 
% First, start with clean system:
clear variables; % free up workspace, clear system memory
close all; % close any open figures
clc; % clear the command window

%% Defining parameters
T0 = 4 ;            %periodic with period T0
T1 = T0/4 ; 
w0 = (2*pi)/T0;     %fundamental angular frequency

points = 201;
lower = 51;         %Lower boundary of the box function
upper = 150;        %Upper boundary of the box function

nmax = 80 ; % limit for the computation of the fourier coefficients

% Given numbers of N for which an approximate box function is calculated:
% the approximate B(t) results for N = 1, 3, 7, 21, 47, 79
N = [1, 3, 7, 21, 47, 79];

linewidth = 2;      % general linewidth used in the plots
%% Plot box function
t =linspace(-T0/2, T0/2, points);   % x-axis for the box function

%Define the box function
y = zeros(points,1);
y(lower: upper) = 1;

%Plot the box function
plot(t, y, 'LineWidth', linewidth)
title('Box function')
xlabel('Time [s]')
ylabel('Amplitude')
grid on
%set the x- & y-axis limits to get better visualization
xlim([-T0/2, T0/2])
ylim([-0.25, 1.25])

%% Dirichlet boundary conditions
disp('Dirichlet his first condition is that the function must be absolutely ') 
disp('integrable over T0. This will be satisfied as the value of the function ') 
disp('is finite for all values of t and the integral over T0 can be computed,') 
disp(' which equals to 2*T1.')
disp('The second condition is that for any time interval, the function is of') 
disp('bounded variation, which means that the total number of maximums and')
disp('minimums is finite during one period of the signal. A box function will')
disp('satisfy this condition.')
disp('The third condition is that in any finite time inteval, the function')
disp('can only have a finite number of finite-jump discontinuities. The box')
disp('function has two discontnuities per period with a jump of value 1 (which') 
disp('is finite). Therefore the box function also satisfies this condition.')
%% Compute Fourier Coefficients
% for n = 1, 2, 3, ..., 80

% a0 = (2*T1)/T0                
% an = sin(n*w0*T1)/(n*pi)

a0 = (2*T1)/T0 ;                                    %eq. 15
an = zeros(nmax, 1); 
for n = 1:nmax
    an(n) = sin(n*w0*T1)/(n*pi);                    %eq. 16
end

%disp(a0)
%disp(an)

%% Approximate B(t) by finite number N
D = zeros(1,length(N));                             %Matrix to store the energy in the differences
figure()
for m = 1:numel(N)                                  % considers all values of N
       B = zeros (points,n);
       for n = 1:N(m)
           B(1:numel(t),n)= an(n)*exp(1i*n*w0*t);   % eq. 17
       end
   B = a0 + 2*sum(B,2);                             % eq. 11
   
   % Calculate the value at ±T1:
   B_T1 = real(B(round(T1*points/T0)+1));
   disp('The value at t=±T1 does not change over N')

   % Calculate the peak value of B:
   B_peak = max(real(B));
   disp('The peak value of B does change over N, the higher N the more the peak')
   disp('value will converge to the peak value of the box function')

   % Compute the energy in the differences
   % the sum of the squares of the local differences between the numerical
   % approximate fourier series representation B(t) and the original box
   % function as a function of N
   diff_B_sq = (real(B) - y).^2;    % squared difference B(t) and the box
   D(m) = sum(diff_B_sq);           % sum all differences and store the value in the vector D

   % Plot B(t) and the box function
   subplot(2,3,m)
   % Original box function
   plot(t, y, 'LineWidth', linewidth)           
   hold on
   % Approximated box function
   plot(t, real(B), 'LineWidth', linewidth)
   title(['N = ', num2str(N(m)), ', B(+/-T_1) = ' num2str(B_T1) ', B_{max}(t) = ' num2str(B_peak)]) 
   xlabel('Time [s]')
   ylabel('Amplitude')
   legend ('Original box', 'Approximated','Location','south')
end

%% Display the energy in the differences
figure ()
plot(N,D,'LineWidth',linewidth)
title('Energy in the difference for different N for T_0 = 4*T_1')
xlabel('N')
ylabel('Energy difference')

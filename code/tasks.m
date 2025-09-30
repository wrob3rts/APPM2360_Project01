% tasks.m
% "Project 1 Task Sets"
% Authors: Ryan Jenkins, Evan Miller, William Roberts
% Date Last Modified: 09/30/2025

%% =================================================================
% --- Task Set B (RK4 Method vs. Exact Solution) ---
% =================================================================
clc;
clear;
close all;

disp('--- Starting Task Set B ---');

%Defining rk4 variables
ti=0;
tf=24;
npts= 240;
y0=50;
f= @(t,w) 0.25*(75-w); 

%Calling rk4 function
[t1,w1] = rk4(ti,tf,npts,y0,f);

%Defining exact solution eqn
T_exact = 75 - 25 * exp(-0.25 * t1);


figure;
%Plotting rk4 approximation
plot(t1, w1, 'ro');
%Graphing exact solution on same plot
hold on;
plot(t1, T_exact, 'b-');
%Labeling plots
xlabel('t');
ylabel('T(t)');
legend('RK4 Approximation', 'Exact Solution');
title('RK4 vs Exact Solution');
grid on;
hold off;

%Displays max calculated temperature
disp(['Max Temperature (RK4): ', num2str(max(w1))]);
disp(' ');

%Plotting error (|exact - rk4 approximation|)
figure;
plot(t1, abs(T_exact-w1));
xlabel('t');
ylabel('| Approximation - Exact |');
legend('Error');
title('Error of RK4 Numerical Approximation')
%% =================================================================
% --- Task Set C (Analytical Solution) ---
% =================================================================
disp('--- Starting Task Set C ---');

%Defining rk4 variables 
ti=0;
tf=24;
npts= 240;
y0=65;
% other equation variables -- M0=75
M0=75
k=0.25
% Equation 3 with specefic conditions
f= @(t,w) k*((M0-12*cos((pi*(t-5))/12))-w);
%Creating vector of t values for graphing M(t)
tC = 0:0.1:24;
%Defining M(t)
M= M0-12*cos((pi*(tC-5))/12); 

%Calling rk4 function
[t1,w1] = rk4(ti,tf,npts,y0,f);

%Plotting rk4 approximation
figure;
plot(t1,w1 , 'ro', 'DisplayName','T(t)');
%Graphing W(t) on same plot
hold on;
plot(tC, M, 'b-', 'DisplayName','Mt(t)');
%labeling plots
xlabel('t'); ylabel('Temperature');
legend();
title('Analytical Temperature Model when H(t)=Q(t)=0, M_0 = 75');
grid on;
hold off;

%Calculating and displaying minimum and maximum of each
disp(['Min M(t): ', num2str(min(M))]);
disp(['Max M(t): ', num2str(max(M))]);
disp(['Min T: ', num2str(min(w1))]);
disp(['Max T: ', num2str(max(w1))]);
disp(' ');



%Repeating for M0=35


% other equation variables -- M0=35
M0=35
k=0.25
% Equation 3 with specefic conditions
f= @(t,w) k*((M0-12*cos((pi*(t-5))/12))-w);
%Creating vector of t values for graphing M(t)
tC2 = 0:0.1:24;
%Defining M(t)
M2= M0-12*cos((pi*(tC-5))/12); 

%Calling rk4 function
[t2,w2] = rk4(ti,tf,npts,y0,f);

%Plotting rk4 approximation
figure;
plot(t2,w2 , 'ro', 'DisplayName','T(t)');
%Graphing W(t) on same plot
hold on;
plot(tC2, M2, 'b-', 'DisplayName','Mt(t)');
%labeling plots
xlabel('t'); ylabel('Temperature');
legend();
title('Analytical Temperature Model when H(t)=Q(t)=0, M_0 = 35');
grid on;
hold off;

%Calculating and displaying minimum and maximum of each
disp(['Min M(t): ', num2str(min(M2))]);
disp(['Max M(t): ', num2str(max(M2))]);
disp(['Min T: ', num2str(min(w2))]);
disp(['Max T: ', num2str(max(w2))]);
disp(' ');

%% =================================================================
% --- Task Set D (Analytical Function) ---
% =================================================================
disp('--- Starting Task Set D ---');

%Defining rk4 variables 
ti=0;
tf=24;
npts= 240;
y0=65;
% Equation 3 with specefic conditions
f= @(t,w) 7*sech((3/4)*(t-10));

%Calling rk4 function
[t1,w1] = rk4(ti,tf,npts,y0,f);

%Plotting rk4 approximation
figure;
plot(t1,w1 , 'ro', 'DisplayName','T(t)');

%labeling plots
xlabel('t'); ylabel('Temperature');
legend();
title('Analytical Temperature Model when A(t)=Q(t)=0');
grid on;
hold off;

%Creating vector of t values for graphing H(t)
tH = 0:0.1:24;
%Defining H(t)
Ht = 7*sech((3/4)*(tH-10));
%Plotting H(t)
figure;
plot(tH, Ht, 'b');
xlabel('t');
ylabel('H(t)');
title('H(t)');
grid on;

%Calculating and displaying maximum and minimum temps
disp(['Max T: ', num2str(max(w1))]);
disp(['Min T: ', num2str(min(w1))]);
disp(' ');

%% =================================================================
% --- Task Set E (Four Analytical Functions) ---
% =================================================================
disp('--- Starting Task Set E ---');

%Defining consistant rk4 variables 
ti=0;
tf=24;
npts= 240;
T_d=77

% Equation 3 with specefic conditions
%T(0)=65, k_d = 0.2:
y0 = 65;
k=0.2
f= @(t,w) k*(T_d - w);

    %Calling rk4 function
[t1,w1] = rk4(ti,tf,npts,y0,f);

%T(0)=65, k_d = 2.0:
y0 = 65;
k=2.0
f= @(t,w) k*(T_d - w);

    %Calling rk4 function
[t2,w2] = rk4(ti,tf,npts,y0,f);


%T(0)=95, k_d = 0.2:
y0 = 95;
k=0.2
f= @(t,w) k*(T_d - w);

    %Calling rk4 function
[t3,w3] = rk4(ti,tf,npts,y0,f);

%T(0)=95, k_d = 2.0:
y0 = 95;
k=2.0
f= @(t,w) k*(T_d - w);

    %Calling rk4 function
[t4,w4] = rk4(ti,tf,npts,y0,f);

%Plotting rk4 approximations
figure;
plot(t1, w1, 'r', 'DisplayName','T(0)=65, k_d=0.2');
hold on;
plot(t2, w2, 'g', 'DisplayName','T(0)=65, k_d=2.0');
plot(t3, w3, 'b', 'DisplayName','T(0)=95, k_d=0.2');
plot(t4, w4, 'y', 'DisplayName','T(0)=95, k_d=2.0');
xlabel('t'); ylabel('Temperature');
legend();
title('Four Analytical Functions when A(t)=H(t)=0');
grid on;
hold off;
disp('Task Set E plots four distinct analytical solutions.');
disp(' ');

%% =================================================================
% --- Task Set F - Q1 (RK4) ---
% =================================================================
disp('--- Starting Task Set F - Q1 ---');

%Defining rk4 variables 
ti=0;
tf=24;
npts= 240;
y0=75;
k=2;
f= @(t,w) 7*sech((3/4)*(t-10))+k*(77-w); 

%Calling rk4 function
[t1,w1] = rk4(ti,tf,npts,y0,f);


%Plotting rk4 approximation
figure;
plot (t1,w1, 'mo')
xlabel('t'); ylabel('T(t)');
legend("RK4 Approximation");
title('Temperature when A(t)=0');
grid on;
disp(['Max T: ', num2str(max(w1))]);
disp(' ');

%% =================================================================
% --- Task Set F - Q2 (RK4) ---
% =================================================================
disp('--- Starting Task Set F - Q2 ---');

%Defining rk4 variables 
ti=0;
tf=24;
npts= 240;
y0=75;
f= @(t,w) 0.25*(85-10*cos((pi*(t-5))/12)-w); 

%Calling rk4 function
[t1,w1] = rk4(ti,tf,npts,y0,f);


%Plotting rk4 approximation
figure;
plot (t1,w1, 'co')
xlabel('t'); ylabel('T(t)');
legend('RK4 Approximation');
title('Temperature when H(t)=Q(t)=0');
grid on;
disp(['Max T: ', num2str(max(w1))]);
disp(' ');

%% =================================================================
% --- Task Set F - Q3 (RK4) ---
% =================================================================
disp('--- Starting Task Set F - Q3 ---');

%Defining consistant rk4 variables 
ti=0;
tf=24;
npts= 240;
y0=75;

%Defining function when k=0.5
k=0.5;
f= @(t,w) 0.25*(85-10*cos((pi*(t-5))/12)-w) + k*(77-w);
%Calling rk4 function 
[t1,w1] = rk4(ti,tf,npts,y0,f);

% rkf approximation when k=0.5
figure;

plot (t1,w1, 'ro');

%Defining function when k=2
k=2;
f2 = @(t,w) 0.25*(85-10*cos((pi*(t-5))/12)-w) + k*(77-w);
%Calling rk4 function
[t2,w2] = rk4(ti,tf,npts,y0,f2);


hold on
plot (t2,w2, 'b');
hold off

xlabel('t'); ylabel('T(t)');
legend('RK4 Approximation k = 0.5','RK4 Approximation k = 2');
title('Temperature when H(t)=0');
grid on;
disp(['Max T(k=0.5): ', num2str(max(w1))]);
disp(['Max T(k=2): ', num2str(max(w2))]);
disp(' ');

%% =================================================================
% --- Task Set F - Q4 (RK4) ---
% =================================================================
disp('--- Starting Task Set F - Q4 ---');

%Defining rk4 variables 
ti=0;
tf=72;
npts= 720;
y0=75;
f= @(t,w) 0.25*(85-10*cos((pi*(t-5))/12)-w) + 7*sech((3/4)*(t-10))+2*(77-w);

%Calling rk4 function
[t1,w1] = rk4(ti,tf,npts,y0,f);


%Plottingg rk4 approximation
figure;
plot (t1,w1,'o')

%Graphing M(t) on same plot
hold on
plot (t1,85-10*cos((pi*(t1-5))/12));
hold off

%Labeling plots
xlabel('t'); 
ylabel('T(t)');
legend('RK4 Approximation', 'M(t)');
title('Long-term Temperature');
grid on;
disp(['Max T: ', num2str(max(w1))]);
disp(' ');

%% =================================================================
% --- RK4 Method (Local Function) ---
% =================================================================
function [t, y] = rk4(t0, tf, n, y0, f)
    h = (tf - t0) / n;
    t = linspace(t0, tf, n+1);
    y = zeros(1, n+1);
    y(1) = y0;

    for i = 1:n
        k1 = f(t(i), y(i));
        k2 = f(t(i) + h/2, y(i) + h*k1/2);
        k3 = f(t(i) + h/2, y(i) + h*k2/2);
        k4 = f(t(i) + h, y(i) + h*k3);
        y(i+1) = y(i) + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end

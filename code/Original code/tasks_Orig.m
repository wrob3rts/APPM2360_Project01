% tasks.m
% "Project 1 Task Sets"
% Authors: Ryan Jenkins, Evan Miller, William Roberts
% Date Last Modified: 09/29/2025

%% =================================================================
% --- Task Set B (RK4 Method vs. Exact Solution) ---
% =================================================================
clc;
clear;
close all;

disp('--- Starting Task Set B ---');
ti=0;
tf=24;
npts= 240;
y0=50;
f= @(t,w) 0.25*(75-w); 

[t1,w1] = rk4(ti,tf,npts,y0,f);


T_exact = 75 - 25 * exp(-0.25 * t1);


figure;
plot(t1, w1, 'ro');
hold on;
plot(t1, T_exact, 'b-');
xlabel('t');
ylabel('T(t)');
legend('RK4 Approximation', 'Exact Solution');
title('RK4 vs Exact Solution');
grid on;
hold off;

disp(['Max Temperature (RK4): ', num2str(max(w1))]);
disp(' ');

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
tC = 0:0.1:24;
M = 75;
T0 = 65;
Mt = M - 12 * cos((pi * (tC - 5)) / 12);
k = 0.25;
T_C = Mt - (Mt - T0) .* exp(-k * tC);

figure;
plot(tC, T_C, 'b-', 'DisplayName','T(t)');
hold on;
plot(tC, Mt, 'ro', 'DisplayName','Mt(t)');
xlabel('t'); ylabel('Temperature');
legend();
title('Analytical Temperature Model');
grid on;
hold off;

disp(['Min Mt: ', num2str(min(Mt))]);
disp(['Max Mt: ', num2str(max(Mt))]);
disp(['Min T: ', num2str(min(T_C))]);
disp(['Max T: ', num2str(max(T_C))]);
disp(' ');

%% =================================================================
% --- Task Set D (Analytical Function) ---
% =================================================================
disp('--- Starting Task Set D ---');
tD = 0:0.1:24;
TD = (56/3) * atan(tanh((3/8) * (tD - 10))) + 79.6504;

figure;
plot(tD, TD, 'b');
xlabel('t');
ylabel('T(t)');
title('Analytical Temperature Function');
grid on;

disp(['Max T: ', num2str(max(TD))]);
disp(['Min T: ', num2str(min(TD))]);
disp(' ');

%% =================================================================
% --- Task Set E (Four Analytical Functions) ---
% =================================================================
disp('--- Starting Task Set E ---');
tE = 0:0.1:24;
Ta = 77 - 12 * exp(-0.2 * tE);
Tb = 77 - 12 * exp(-2 * tE);
Tc = 77 + 18 * exp(-0.2 * tE);
Td = 77 + 18 * exp(-2 * tE);

figure;
plot(tE, Ta, 'r', 'DisplayName','Ta');
hold on;
plot(tE, Tb, 'g', 'DisplayName','Tb');
plot(tE, Tc, 'b', 'DisplayName','Tc');
plot(tE, Td, 'y', 'DisplayName','Td');
xlabel('t'); ylabel('Temperature');
legend();
title('Four Analytical Functions');
grid on;
hold off;
disp('Task Set E plots four distinct analytical solutions.');
disp(' ');

%% =================================================================
% --- Task Set F - Q1 (RK4) ---
% =================================================================
disp('--- Starting Task Set F - Q1 ---');

ti=0;
tf=24;
npts= 240;
y0=75;
f= @(t,w) 7*sech((3/4)*(t-10))+2*(77-w); 

[t1,w1] = rk4(ti,tf,npts,y0,f);


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
ti=0;
tf=24;
npts= 240;
y0=75;
f= @(t,w) 0.25*(85-10*cos((pi*(t-5))/12)-w); 

[t1,w1] = rk4(ti,tf,npts,y0,f);


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

ti=0;
tf=24;
npts= 240;
y0=75;
f= @(t,w) 0.25*(85-10*cos((pi*(t-5))/12)-w) + 0.5*(77-w);

[t1,w1] = rk4(ti,tf,npts,y0,f);


figure;

plot (t1,w1, 'ro');


f2 = @(t,w) 0.25*(85-10*cos((pi*(t-5))/12)-w) + 2*(77-w);
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
ti=0;
tf=72;
npts= 720;
y0=75;
f= @(t,w) 0.25*(85-10*cos((pi*(t-5))/12)-w) + 7*sech((3/4)*(t-10))+2*(77-w);


[t1,w1] = rk4(ti,tf,npts,y0,f);


figure;
plot (t1,w1,'o')
hold on
plot (t1,85-10*cos((pi*(t1-5))/12));
hold off

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

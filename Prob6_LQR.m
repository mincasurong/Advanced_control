%% Advanced control I MATLAB
% Tutorial code for the lecture in department of mechanical engineering,
% Sogang university, Korea, Republic of
% Made by Minsu Chang, Ph.D candidate
clear all; close all; clc;

%% System model

% Transfer function in continuous time
G = tf(1.295,[1,6,34]) 
[b, a] = tfdata(G,'v');

% State space representation in observable canonical form
A = [-a(3) 1; -a(2) 0];
B = [b(3); b(2)]; 
C = [1 0]; 
D = 0;
sys = ss(A, B, C, D); % === Empty

% Sampling period: 1 ms
T = 0.001; 

% Transfer function in discrete time 
sysz = c2d(sys,T,'zoh') 

Ad = sysz.a;
Bd = sysz.b;
Cd = sysz.c;
Dd = sysz.d;

% Controllability matrix
rank(ctrb(Ad, Bd))

%% Define controller and estimator
% Design controller based on digital LQR
r = 0.01;  % ============================================= Most important 
Q = Cd' * Cd;  % === Empty
Klq = dlqr(Ad, Bd, Q, r)   % === Empty

% estimation gain of state observer
L = place(Ad', Cd', [0.01, 0.011])' % ==== estimator gain is also important

%% Define variables
t = 0:T:0.5; % Time
N = length(t);

x1 = zeros(N,1);   x2 = zeros(N,1);   % state variables
xh1 = zeros(N,1);  xh2 = zeros(N,1);  % estimated state
u1 = zeros(N,1);   u2 = zeros(N,1);  u = zeros(N,1);  % Control input
y = zeros(N,1);    yh = zeros(N,1); % Output

x1(1) = 1; % Initial value

%% Simulation

for k=1:N-1
    u1(k) = Klq(1)*xh1(k); 
    u2(k) = Klq(2)*xh2(k);
    u(k) = u1(k) + u2(k);
    
    % Actual system
    x1(k+1) = Ad(1,1)*x1(k) + Ad(1,2)*x2(k) - Bd(1)*u1(k) - Bd(1)*u2(k) + Bd(1)*0*(1+sin(2*pi*10*t(k)));
    x2(k+1) = Ad(2,1)*x1(k) + Ad(2,2)*x2(k) - Bd(2)*u1(k) - Bd(2)*u2(k);
    y(k) = Cd(1)*x1(k) + Cd(2)*x2(k) + 0*randn(1,1);
    
    % Estimated system based on state observer
    xh1(k+1) = (Ad(1,1)-Bd(1)*Klq(1))*xh1(k) + (Ad(1,2)-Bd(1)*Klq(2))*xh2(k) + ...
               L(1)*( y(k) - ( Cd(1)*xh1(k) + Cd(2)*xh2(k) ) );
    xh2(k+1) = (Ad(2,1)-Bd(2)*Klq(1))*xh1(k) + (Ad(2,2)-Bd(2)*Klq(2))*xh2(k) + ...
               L(2)*( y(k) - ( Cd(1)*xh1(k) + Cd(2)*xh2(k) ) );
    yh(k) = Cd(1)*xh1(k) + Cd(2)*xh2(k);
end

%% Graph

figure('color','w')
subplot(211);
plot(t,y,'b','linewidth',2); hold on;
plot(t,yh,'r:','linewidth',2);
xlabel('time (sec.)'); 
ylabel('Output (deg)');
legend('Actual output','Estimated output')
axis([0 t(end) -2 2]); grid on;

subplot(212);
plot(t,u ,'b','linewidth',2); hold on;
xlabel('time (sec.)'); 
ylabel('Input (V)');
xlim([0 t(end)]); grid on;



clear all;
close all;
clc;

T = 0.001;                      % Sampling period
w = logspace(-1,1,100)*2*pi;    % Frequency samples for obtaining Bode plots

%% System identification
t = 0:0.001:5;
N = length(t);
y = zeros(1,N);
u = zeros(1,N);
r= zeros(1,N);
e=zeros(1,N);
dhat = zeros(1,N);  % Estimated disturbance
ua = zeros(1,N);    % Actual input incluing a disturnace and an estimated disturbance
d = zeros(1,N);     % Disturbance
yd = zeros(1,N);    % Desired output

num = [719.3  1.208e04];
den = [1 101.2 1956 4805];
Gs = tf(num,den)
Gz = c2d(Gs,T)
F=[-101.2 1 0;
   -1956 0 1;
   -4805 0 0];
G=[0; 719.3; 12080];
H=[1 0 0];
J=0;

sys=ss(F,G,H,J);
sysz=c2d(sys,T);
phi=sysz.a
Gamma=sysz.b
H=sysz.c
J=sysz.d

[a b]=tfdata(Gz,'v')

%% Trajectory design
% Design of a desired output
for k=1:N;
yd(k) = 1*sin(0.5*2*pi*k*T)-0.25*cos(2*pi*k*T);    % A point-to-point trajectory
end
% figure('position',[410 420 400 300]);
% plot(t,yd,'r');

% Smoothening of the desired output
% Gsmooth = tf([50],[1 50]);  % A filter for smoothning (cut-off frequency is 100 rad/s)
% Gsmoothz = c2d(Gsmooth,T,'tustin');  % Tustin = bilinear method
% [A,B] = tfdata(Gsmoothz,'v');
% yd = filtfilt(A,B,yd);  % Smoothning
% hold on;
% plot(t,yd,'k');
% title('Desired outputs');
% axis([0 10 -0.001 0.021]);

%% Controller design
% Setting feedback control (PD) gains
Kp=50; Kd=0.1
C=tf([Kd Kp],[1]);
Cz=c2d(C,T,'matched')
Gcz=feedback(Gz*Cz,1)

% Feedforward filter design

Gcz=minreal(Gcz)
inv_Gcz=minreal(Gcz^-1)
[numz denz]=tfdata(Gcz,'v')
[inumz idenz]=tfdata(inv_Gcz,'v')

% Design of a Q filter (Notice that the relative order of Gdnew(z) is 2)
CutoffFreq = 35*2*pi;    % Cutoff frequency of Q filter
Q = zpk([],[exp(-CutoffFreq*T) exp(-CutoffFreq*T)],1,T);    % Pole-zero matching between Q(s) and Q(z) 
Q = Q / dcgain(Q);      % To make the Q(exp(j0T)) = 1
[aQ,bQ] = tfdata(Q,'v');    % Model parameters of Q filter

% Design of DOB
GDOB = minreal(Q*Gz^-1);   
[aD,bD] = tfdata(GDOB,'v');    % Model parameters of Q*G^-1
% Notice that bD and bQ are the same

%% Simulation
for k = 5:N-2,
    % Simulation by G(z), i.e., actual dynamics
    y(k) = -b(2)*y(k-1)-b(3)*y(k-2)-b(4)*y(k-3)+a(1)*ua(k)+a(2)*ua(k-1)+a(3)*ua(k-2)+a(4)*ua(k-3);

    r(k)=-(idenz(2)*r(k-1)+idenz(3)*r(k-2)+idenz(4)*r(k-3))+inumz(1)*yd(k)+inumz(2)*yd(k-1)+inumz(3)*yd(k-2)+inumz(4)*yd(k-3);

    e(k)=yd(k)-y(k);
    % Disturbance observer
%     dhat(k) = - bQ(2)*dhat(k-1) - bQ(3)*dhat(k-2) + (aD(1)*y(k) + aD(2)*y(k-1) + aD(3)*y(k-2)+aD(4)*y(k-3)+aD(5)*y(k-4)) - (aQ(1)*u(k) + aQ(2)*u(k-1) + aQ(3)*u(k-2));
    dhat(k)=+1.6162681719029*dhat(k-1)+0.294116295467*dhat(k-2)+-1.5275146845572*dhat(k-3)+0.61585019396264*dhat(k-4)+111.43245178213*y(k-1)+-323.36471327886*y(k-2)+312.64010126596*y(k-3)+-100.70733062093*y(k-4)+-0.038970773065752*u(k-2)+0.00043212745584801*u(k-3)+0.037258622385322*u(k-4);
    % Feedback control
    uc(k) = Kp*e(k) + Kd*(e(k)-e(k-1))/T;   % PD control
    u(k) = uc(k) - dhat(k);     % Disturbance rejection
    
    % Saturation of control input
    if u(k) > 10,
        u(k) = 10;
    end
    if u(k) < -10,
        u(k) = -10;
    end
    
    % Disturbance (it is unknown in practice)
    d(k) = sin(5*pi*k*T)+cos(0.2*pi*k*T)+exp(0.1*T*cos(k*2*pi));    % If you do not want to inject a disturbance, set d(k) = 0;
    ua(k) = u(k) + d(k);     % Actual input with a disturbance. The control algorithm has no information on d(k)
end

% figure(2);
% plot(t,r,'b');
% legend('Unsmoothened','Smoothened','Feedforward filtered');

figure('position',[10 50 400 300]);
subplot(313)
plot(t,yd,'r*',t,y,'b','linewidth',2);
xlabel('Time (sec.)');
ylabel('Position(m)');
% title('Position control result');
legend('Desired output','Model output')
grid on
subplot(312)
plot(t,u,'linewidth',2);
xlabel('Time (sec.)');
ylabel('Voltage(V)');
% title('Input signals');
legend('Control input');
grid on
subplot(311)
plot(t,dhat,'r',t,d,'b','linewidth',2);
xlabel('Time (sec.)');
ylabel('Magnitude');
% title('Input signals');
legend('Estimated disturbance','Disturbance');
grid on
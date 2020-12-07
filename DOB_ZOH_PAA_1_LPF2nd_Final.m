% This system is designed by MS. C, 20.9.24
% Disturbance observer of SUBAR's legf leg
% Comparison of Experimental data version

clear all;  close all;  clc;

global num den a b aQ bQ Kp Kd spd str kh yd bD aD Joint Simul TFmodel Tfinal PAA DOB_Rejection var M C K an bn w
Tfinal = 20;
Simul = 1;         % 0: Experimental data, 1: Simulation data
TFmodel = 0;       % 0: Experimenta model, 1: Adopted model
PAA = 0;           % 0: OFF, 1: Four parameter
DOB_Rejection = 1; % 0: PD control only, 1: Disturbance rejection ON
%% Trajectory design
switch(Simul)
    case 0,
        
        P = data(:,1);
        t_new = data(:,2);
        Trj = data(:,3);
        Ang = data(:,4);
        dq = data(:,5);
        ud = data(:,9);
        yd = Trj(:,1);
        
        T = 0.001;
        N = length(data);
        t = linspace(0,(N-1)*T,N);
        
    case 1,
        T = 0.001;
        t = 0:T:Tfinal;
        yd = 20*sin(2*pi*1*t);
end

%% Initialization
N = length(t); y = zeros(1,N); ud = zeros(1,N); e=zeros(1,N);
ua = zeros(1,N);    % Actual input incluing a disturnace and an estimated disturbance
uc = zeros(1,N);    % controller output
d1 = zeros(1,N);     d2 = zeros(1,N);     % Disturbance
dy = zeros(1,N);    % ADOB
dhat = zeros(1,N);  % Estimated disturbance
d = zeros(1,N);     % Disturbance
duq = zeros(1,N);     % Disturbance
dyq = zeros(1,N);     % Disturbance
wf = logspace(-2,2,100);    % Frequency samples for obtaining Bode plots


%% Transfer function
num = [0,0,1.29497547158315]; den = [1,6,34];
Gain = num;
M = den(1); C = den(2); K = den(3);

% Zero-order-holder equivalent
G = tf(num,den)
Gz = c2d(G,T,'zoh')

[a b] = tfdata(Gz,'v')

% After Parameter adaptation algorithm (Gzz)
switch(TFmodel)
    case 0, an = a; bn = b;  % Original nominal model
    case 1,
        numz = [0,2.14030136941629e-06,7.79441524828860e-08]; % Adopted nominal model
        denz = [1,-1.99417566310256,0.994511602123872];
        Gzz = tf(numz,denz,T);
        [an bn] = tfdata(Gzz,'v');
end
switch(PAA)
    case 1, theta = zeros(4,1); phi = zeros(4,N); var = zeros(4,N);
end

%% Controller & DoB design
% Setting feedback control (PD) gains
Kp=500; Kd=4;
C=tf([Kd Kp],[1]);
Cz=c2d(C,T,'matched');
Gcz=feedback(Gz*Cz,1);

% Design of a Q filter (Notice that the relative order of Gdnew(z) is 2)
w = 50*2*pi;  Q1 = tf([w],[1 w]);

syms z; z =tf('z',T); alpha = exp(-w*T);  QQ = ((1-alpha)/2*(z+1)/(z-alpha))^2;

Q = Q1*Q1;
Qz = c2d(Q,T,'matched')
Qz = Qz / dcgain(Qz);      % To make the Q(exp(j0T)) = 1
[aQ,bQ] = tfdata(Qz,'v');    % Model parameters of Q filter
% Design of DOB
GDOB = minreal(Qz*Gz^-1);
[aD,bD] = tfdata(GDOB,'v');    % Model parameters of Q*G^-1


%% Disturbance
v = zeros(1,N);
for k=2:N, v(k) = (yd(k)-yd(k-1))/T; end
for k=1:N
    g = 9.81;
    gear_ratio = 1; l = 5;    m = 20;
    d1(k) = -m*g*l*gear_ratio*sin(yd(k)*10*pi/180); % Gravity
        d2(k) = 800*(sin(2*pi*0.3*t(k))) + 500*(sin(2*pi*1.3*t(k))) + 200*(sin(2*pi*2*t(k)));
    d(k) = d1(k) + d2(k);
end

%% Simulation
y(1:3) = 0; iter = 0;
for k=1:4,
    switch(PAA),
        case 1, var(:,k) = [b(2) b(3) a(2) a(3)];   F = eye(4)*10^10; ramda = 0.8;
    end;
end
for k = 4:N,
    % Simulation by G(z), i.e., actual dynamics
    y(k) = -b(2)*y(k-1) -b(3)*y(k-2) +a(1)*ua(k) +a(2)*ua(k-1) +a(3)*ua(k-2);
    e(k) = yd(k) - y(k);
    
    dhat(k) = 1/(an(2)*bQ(1))*(...
        (aQ(2)*bn(1)*y(k) + (aQ(3)*bn(1) + aQ(2)*bn(2))*y(k-1) + (aQ(2)*bn(3)+aQ(3)*bn(2))*y(k-2) + aQ(3)*bn(3)*y(k-3) ) +...
        - ( an(2)*aQ(2)*ud(k-1) +(an(3)*aQ(2) + an(2)*aQ(3))*ud(k-2) + an(3)*aQ(3)*ud(k-3)   ) +...
        - ( an(3)*bQ(1) + an(2)*bQ(2) )*dhat(k-1) - (an(2)*bQ(3)+an(3)*bQ(2))*dhat(k-2) - an(3)*bQ(3)*dhat(k-3) ...
        );
    
    % Control
    uc(k) = Kp*e(k) + Kd*(e(k)-e(k-1))/T;   % PD control
    switch(DOB_Rejection)
        case 0, ud(k) = uc(k);     % Disturbance rejection
        case 1, ud(k) = uc(k) - dhat(k);     % Disturbance rejection
    end
    % Saturation of control input
    if ud(k) > 5000, ud(k) = 5000; end;  if ud(k) < -5000, ud(k) = -5000; end
    
    switch(PAA)
        case 1,
            phi(:,k) = [-y(k-1) -y(k-2) ud(k-1) ud(k-2)]'; % adaptive algorithm
            F = F - (F*phi(:,k)*phi(:,k)'*F)/(1+phi(:,k)'*F*phi(:,k));
            output = y(k);
            theta = theta + F*phi(:,k)*(output - theta'*phi(:,k));
            var(:,k) = theta;
            iter = 1;
    end
    
    % Disturbance (it is unknown in practice)
    ua(k) = ud(k) + d(k);     % Actual input with a disturbance. The control algorithm has no information on d(k)
end

%% Figure
figure('color','w')
subplot(511)
plot(t,yd,'r','linewidth',2); hold on;
plot(t,y,'b:','linewidth',2); hold on;
ylabel('Position(m)');
legend('Desired output','Model output','Exp output')
grid on
axis([0 t(end) min(yd)-1 max(yd)+1])

subplot(512)
plot(t,e,'b','linewidth',2); hold on;
% plot(t,err(:,1),'r:','linewidth',2);
ylabel('Position(m)');
legend('Feedback error','Ex error');
grid on
axis([0 t(end) -2 2])

subplot(513)
plot(t,ud*0.001,'b','linewidth',2); hold on;
plot(t,uc*0.001,'g:','linewidth',2); hold on;
% plot(t,udata(:,1),'r:','linewidth',2); hold on;
ylabel('Voltage(V)');
legend('Control input','PD controller output','Exp output');
grid on
axis([0 t(end) -5 5])

subplot(514)
plot(t,dhat*0.001,'r','linewidth',3); hold on;
plot(t,d*0.001,'b','linewidth',2); hold on;
plot(t,d1*0.001,'c:','linewidth',2); hold on;
plot(t,d2*0.001,'m:','linewidth',2); hold on;
ylabel('Magnitude');
legend('Estimated disturbance','Disturbance','gravity','external');
grid on
axis([0 t(end) -5 5])

subplot(515)
de = (d-dhat)*0.001;
plot(t,de,'b','linewidth',2);
ylabel('Magnitude');
legend('d-dhat');
grid on
axis([0 t(end) -5 5])
drawnow;
%% Transfer function after PAA

switch(PAA)
    case 1, Gzz = minreal(tf([0 var(3:4,end)'],[1 var(1:2,end)'],T));
end

if PAA > 0
    w = logspace(-2,2,1000);
    t1 = 0:0.001:10;
    figure('color','w');
    subplot(221); step(Gz,t1); hold on; step(Gzz,t1)
    subplot(222); impulse(Gz,t1); hold on; impulse(Gzz,t1)
    subplot(223); bode(Gz,w); hold on; bode(Gzz,w)
    subplot(224); pzmap(Gz); hold on; pzmap(Gzz)
end

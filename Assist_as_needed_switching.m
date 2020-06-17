clear all; close all; clc;

% System of Position Feedback
global animation Assist
animation = 100;
%% Initial Setting
T=0.001;
t=0:T:20;
N=length(t);

y=zeros(1,N); yd= zeros(1,N); e=zeros(1,N); r=zeros(1,N); yimp=zeros(N,1);
u=zeros(1,N); uc=zeros(1,N); u1=zeros(1,N); u2=zeros(1,N); ua=zeros(1,N); ur=zeros(1,N);
d=zeros(1,N); h1=zeros(1,N); h2=zeros(1,N); h=zeros(1,N); v=zeros(1,N); eh = zeros(1,N);
ua_resist=zeros(1,N);  ua_assist=zeros(1,N); n = zeros(1,N);
Beta=zeros(1,N); rv = zeros(1,N);
dhat1 = zeros(1,N);  dhat2 = zeros(1,N); dhat = zeros(1,N);

%% System Identification
m = 2; k = 1;
num=[k/m]; den=[1 2*sqrt(k/m) k/m];
Gs=tf(num,den);
Gz=c2d(Gs,T);
[a b]=tfdata(Gz,'v');

% PD Controller design
Kp=150; Kd=10;


%% Desired Position
A = 0.3;
for k=1:N;
    r(k)=A*(cos(0.5*pi*t(k))-1);
end


%% Disturbance
B = 1;
for k=1:N;
    d(k)=B*(cos(1*pi*t(k))-1);
end

%% DoB
% Design of a Q filter (Notice that the relative order of Gdnew(z) is 2)
CutoffFreq = 35*2*pi;    % Cutoff frequency of Q filter
Q = zpk([],[exp(-CutoffFreq*T) exp(-CutoffFreq*T)],1,T);    % Pole-zero matching between Q(s) and Q(z)
Q = Q / dcgain(Q);      % To make the Q(exp(j0T)) = 1
[aQ,bQ] = tfdata(Q,'v');    % Model parameters of Q filter

% Design of DOB
GDOB = minreal(Q*Gz^-1);
[aD,bD] = tfdata(GDOB,'v');    % Model parameters of Q*G^-1

%% Hybrid controller
Beta_t1 = 10;
Beta_t2 = 20;
Beta_td = Beta_t2 - Beta_t1;

for k=1:N, 
    if t(k)>=Beta_t1 && t(k)<Beta_t2
        Beta(k) = 0.5*(1-cos(pi*(t(k)-Beta_t1)/Beta_td));
    elseif t(k)>=Beta_t2
        Beta(k) = 1;
    end
    
    if Beta(k) <0, Beta(k) = 0; elseif Beta(k) > 1, Beta(k) = 1; end
end

%% Human controller
Hp = 30;
Hd = 1;

j = 1;
%% Simulation
for k=4:N-1,
    % measurement noise
    n(k) = 0.000001*sin(80*t(k));
    
    % Transfer function (Gz)
    y(k) = -b(2)*y(k-1)-b(3)*y(k-2)+a(1)*u(k)+a(2)*u(k-1)+a(3)*u(k-2) + n(k);
    
    % Velocity
    v(k) = (y(k)-y(k-1))/T;
%     if Beta(k) == 0, 
    
    % Virtual trajectory
    rv(k) = (1-Beta(k))*y(k) + Beta(k)*r(k);
    
    % Feedback error & PD controller
    e(k)=rv(k)-y(k);
    u1(k)=Kp*Beta(k)*e(k);             if u1(k) > 5, u1(k) = 5; end; if u1(k) < -5, u1(k) = -5; end
    u2(k)=Kd*Beta(k)*(e(k)-e(k-1))/T;  if u2(k) > 5, u2(k) = 5; end; if u2(k) < -5, u2(k) = -5; end   
    uc(k)=u1(k)+u2(k);                 
    if uc(k) > 5, uc(k) = 5; end; if uc(k) < -5, uc(k) = -5; end
    
    % PD controller
    u1(k)=Kp*Beta(k)*e(k);
    u2(k)=Kd*Beta(k)*(e(k)-e(k-1))/T;
    uc(k)=u1(k)+u2(k);
    if uc(k) > 5, uc(k) = 5; end; if uc(k) < -5, uc(k) = -5; end
    if u1(k) > 5, u1(k) = 5; end; if u1(k) < -5, u1(k) = -5; end
    if u2(k) > 5, u2(k) = 5; end; if u2(k) < -5, u2(k) = -5; end
    
    % DoB
    dhat1(k) = -bD(2)*dhat1(k-1) -bD(3)*dhat1(k-2) - bD(4)*dhat1(k-3) + aD(1)*y(k) + aD(2)*y(k-1) + aD(3)*y(k-2) + aD(4)*y(k-3);
    dhat2(k) = -bQ(2)*dhat2(k-1) -bQ(3)*dhat2(k-2) + aQ(1)*ua(k) + aQ(2)*ua(k-1) + aQ(3)*ua(k-2);
%     dhat(k) = dhat1(k)-dhat2(k);
    
    % Actuator force (Tracking)
    ua_assist(k) = uc(k)*30/5*13.7*0.001*100 - dhat(k);
    
    % Actuator force (Virtual damping)
    ua_resist(k) = -100*v(k);
    
    % Hybrid controller
    Assist = 2; % 0: Resist, 1: Assist, 2: Hybrid
    switch(Assist)
        case 0, ua(k) = ua_resist(k);
        case 1, ua(k) = ua_assist(k);
        case 2, ua(k) = (1-Beta(k))*ua_resist(k) + Beta(k)*ua_assist(k);
    end
    
    % Human force
    eh(k) = r(k)-y(k);
    h1(k)=Hp*eh(k);
    h2(k)=Hd*(eh(k)-eh(k-1))/T;
    h(k)=h1(k)+h2(k);
    
    u(k) = ua(k) + d(k) + h(k);
end


%% Data Plotting

figure('color','w');
subplot(511); plot(t,r,'b','linewidth',2); hold on; hold on; plot(t,y,'r:','linewidth',2);
legend('yd','y');
ylabel('Position (m)');
axis([0 t(end) -1 1])

subplot(512); plot(t,v,'b','linewidth',2);
ylabel('Velocity (m/s)');
axis([0 t(end) -3 3])

subplot(513); plot(t,e,'b','linewidth',2); legend('Position Error')
ylabel('Position Error (m)');
axis([0 t(end) -1 1])

subplot(514); 
plot(t,u1,'c','linewidth',2); hold on; 
plot(t,u2,'m','linewidth',2); hold on; 
plot(t,uc,'b','linewidth',2);
ylabel('Input(V)'); legend('u_p','u_d','u');
axis([0 t(end) -6 6])

subplot(515);   
plot(t,ua,'b','linewidth',2);  hold on; 
plot(t,d,'r','linewidth',2);  hold on; 
plot(t,dhat,'k:','linewidth',2);  hold on; 
plot(t,h,'m:','linewidth',2);  hold on;
legend('Actuation force','Disturbance','DOB','Human force')
ylabel('Force (N)');
xlabel('Time(sec)')
axis([0 t(end) -50 50])

for k=1:5
    subplot(5,1,k); set(gca,'fontsize',10);
end
%% Animation
N=length(t);
figure('color','w');
for k=1:animation:N;
    % Moving Block
    subplot(211);
    plot([r(k)-0.1 r(k)+0.1],[0 0],'r','linewidth',2); hold on;
    plot([r(k)-0.1 r(k)+0.1],[0.5 0.5],'r','linewidth',2); hold on;
    plot([r(k)+0.1 r(k)+0.1],[0 0.5],'r','linewidth',2); hold on;
    plot([r(k)-0.1 r(k)-0.1],[0 0.5],'r','linewidth',2); hold on;
    plot([y(k)-0.1 y(k)+0.1],[0 0],'b:','linewidth',2); hold on;
    plot([y(k)-0.1 y(k)+0.1],[0.5 0.5],'b:','linewidth',2); hold on;
    plot([y(k)+0.1 y(k)+0.1],[0 0.5],'b:','linewidth',2); hold on;
    plot([y(k)-0.1 y(k)-0.1],[0 0.5],'b:','linewidth',2); hold on;
    plot([-1 1],[0 0],'k');
    title(sprintf('Desired Position: %.3f (m). Actual Position: %.3f (m). Error: %.3f (m)',r(k),y(k),r(k)-y(k)))
    axis([-1 1 -1 2]);
    ylabel('y-axis'); xlabel('x-axis');     drawnow;
    hold off;
    
    % Graph
    subplot(212);
    plot(t(1:k),r(1:k),'r','linewidth',2); hold on;
    plot(t(1:k),y(1:k),'b:','linewidth',2);
    ylabel('Position (m)'); xlabel('time (sec)')
    axis([0 t(end) -10 10]);     drawnow;
    hold off;
end
legend('Referece (r)','Output (y)')
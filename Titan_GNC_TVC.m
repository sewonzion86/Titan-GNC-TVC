clear; clc; close all;

gTitan = 1.352;
rho0 = 5.43;
Hrho = 22000;
p0 = 146700;

m0 = 1200;
mDry = 700;

Tmax = 10000;
Isp = 310;
g0 = 9.80665;

Cd = 0.8;
Aref = 0.9;

LTVC = 1.0;
Iy = 1200;

deltaMax = deg2rad(12);
tauTVC = 0.08;

wn = 1.5;
zeta = 0.85;

Kx = 0.010;
Kvx = 0.20;
Kv = 0.80;

dt = 0.05;
tMax = 1000;

t = 0:dt:tMax;
N = numel(t);

state = [15000; 30000; -35; -180; deg2rad(8); 0; 0; m0];

X = zeros(8,N);
U = zeros(5,N);

for k = 1:N
    x = state(1);
    h = state(2);
    Vx = state(3);
    Vz = state(4);
    theta = state(5);
    q = state(6);
    delta = state(7);
    mass = state(8);

    V = hypot(Vx,Vz);
    rho = rho0*exp(-max(h,0)/Hrho);

    if V > 1e-6
        D = 0.5*rho*V^2*Cd*Aref;
        Dx = -D*Vx/V;
        Dz = -D*Vz/V;
    else
        Dx = 0;
        Dz = 0;
    end

    VzCmd = -min(180,Kv*sqrt(max(2*gTitan*h,0)));
    azCmd = Kv*(VzCmd-Vz);
    axCmd = Kx*(0-x)-Kvx*Vx;

    TxReq = mass*axCmd-Dx;
    TzReq = mass*(azCmd+gTitan)-Dz;
    Treq = hypot(TxReq,TzReq);
    T = min(Treq,Tmax);

    thetaCmd = atan2(TxReq,max(TzReq,1e-6));
    thetaError = atan2(sin(thetaCmd-theta),cos(thetaCmd-theta));
    Mcmd = Iy*(wn^2*thetaError-2*zeta*wn*q);

    if T > 1e-6
        deltaCmd = asin(max(-1,min(1,Mcmd/(T*LTVC))));
    else
        deltaCmd = 0;
    end

    deltaCmd = max(-deltaMax,min(deltaMax,deltaCmd));
    alpha = theta + delta;

    Tx = T*sin(alpha);
    Tz = T*cos(alpha);
    ax = (Tx+Dx)/mass;
    az = (Tz+Dz)/mass-gTitan;
    M = T*sin(delta)*LTVC;
    qDot = M/Iy;
    deltaDot = (deltaCmd-delta)/tauTVC;
    mdot = -T/(Isp*g0);

    if mass + mdot*dt < mDry
        mdot = (mDry-mass)/dt;
    end

    k1 = [Vx; Vz; ax; az; q; qDot; deltaDot; mdot];

    s2 = state + 0.5*dt*k1;
    x = s2(1); h = s2(2); Vx = s2(3); Vz = s2(4);
    theta = s2(5); q = s2(6); delta = s2(7); mass = s2(8);
    V = hypot(Vx,Vz);
    rho = rho0*exp(-max(h,0)/Hrho);

    if V > 1e-6
        D = 0.5*rho*V^2*Cd*Aref;
        Dx = -D*Vx/V;
        Dz = -D*Vz/V;
    else
        Dx = 0; Dz = 0;
    end

    VzCmd = -min(180,Kv*sqrt(max(2*gTitan*h,0)));
    azCmd = Kv*(VzCmd-Vz);
    axCmd = Kx*(0-x)-Kvx*Vx;
    TxReq = mass*axCmd-Dx;
    TzReq = mass*(azCmd+gTitan)-Dz;
    Treq = hypot(TxReq,TzReq);
    T = min(Treq,Tmax);
    thetaCmd = atan2(TxReq,max(TzReq,1e-6));
    thetaError = atan2(sin(thetaCmd-theta),cos(thetaCmd-theta));
    Mcmd = Iy*(wn^2*thetaError-2*zeta*wn*q);

    if T > 1e-6
        deltaCmd = asin(max(-1,min(1,Mcmd/(T*LTVC))));
    else
        deltaCmd = 0;
    end

    deltaCmd = max(-deltaMax,min(deltaMax,deltaCmd));
    alpha = theta+delta;
    Tx = T*sin(alpha);
    Tz = T*cos(alpha);
    ax = (Tx+Dx)/mass;
    az = (Tz+Dz)/mass-gTitan;
    M = T*sin(delta)*LTVC;
    qDot = M/Iy;
    deltaDot = (deltaCmd-delta)/tauTVC;
    mdot = -T/(Isp*g0);

    if mass + mdot*dt < mDry
        mdot = (mDry-mass)/dt;
    end

    k2 = [Vx; Vz; ax; az; q; qDot; deltaDot; mdot];

    s3 = state + 0.5*dt*k2;
    x = s3(1); h = s3(2); Vx = s3(3); Vz = s3(4);
    theta = s3(5); q = s3(6); delta = s3(7); mass = s3(8);
    V = hypot(Vx,Vz);
    rho = rho0*exp(-max(h,0)/Hrho);

    if V > 1e-6
        D = 0.5*rho*V^2*Cd*Aref;
        Dx = -D*Vx/V;
        Dz = -D*Vz/V;
    else
        Dx = 0; Dz = 0;
    end

    VzCmd = -min(180,Kv*sqrt(max(2*gTitan*h,0)));
    azCmd = Kv*(VzCmd-Vz);
    axCmd = Kx*(0-x)-Kvx*Vx;
    TxReq = mass*axCmd-Dx;
    TzReq = mass*(azCmd+gTitan)-Dz;
    Treq = hypot(TxReq,TzReq);
    T = min(Treq,Tmax);
    thetaCmd = atan2(TxReq,max(TzReq,1e-6));
    thetaError = atan2(sin(thetaCmd-theta),cos(thetaCmd-theta));
    Mcmd = Iy*(wn^2*thetaError-2*zeta*wn*q);

    if T > 1e-6
        deltaCmd = asin(max(-1,min(1,Mcmd/(T*LTVC))));
    else
        deltaCmd = 0;
    end

    deltaCmd = max(-deltaMax,min(deltaMax,deltaCmd));
    alpha = theta+delta;
    Tx = T*sin(alpha);
    Tz = T*cos(alpha);
    ax = (Tx+Dx)/mass;
    az = (Tz+Dz)/mass-gTitan;
    M = T*sin(delta)*LTVC;
    qDot = M/Iy;
    deltaDot = (deltaCmd-delta)/tauTVC;
    mdot = -T/(Isp*g0);

    if mass + mdot*dt < mDry
        mdot = (mDry-mass)/dt;
    end

    k3 = [Vx; Vz; ax; az; q; qDot; deltaDot; mdot];

    s4 = state + dt*k3;
    x = s4(1); h = s4(2); Vx = s4(3); Vz = s4(4);
    theta = s4(5); q = s4(6); delta = s4(7); mass = s4(8);
    V = hypot(Vx,Vz);
    rho = rho0*exp(-max(h,0)/Hrho);

    if V > 1e-6
        D = 0.5*rho*V^2*Cd*Aref;
        Dx = -D*Vx/V;
        Dz = -D*Vz/V;
    else
        Dx = 0; Dz = 0;
    end

    VzCmd = -min(180,Kv*sqrt(max(2*gTitan*h,0)));
    azCmd = Kv*(VzCmd-Vz);
    axCmd = Kx*(0-x)-Kvx*Vx;
    TxReq = mass*axCmd-Dx;
    TzReq = mass*(azCmd+gTitan)-Dz;
    Treq = hypot(TxReq,TzReq);
    T = min(Treq,Tmax);
    thetaCmd = atan2(TxReq,max(TzReq,1e-6));
    thetaError = atan2(sin(thetaCmd-theta),cos(thetaCmd-theta));
    Mcmd = Iy*(wn^2*thetaError-2*zeta*wn*q);

    if T > 1e-6
        deltaCmd = asin(max(-1,min(1,Mcmd/(T*LTVC))));
    else
        deltaCmd = 0;
    end

    deltaCmd = max(-deltaMax,min(deltaMax,deltaCmd));
    alpha = theta+delta;
    Tx = T*sin(alpha);
    Tz = T*cos(alpha);
    ax = (Tx+Dx)/mass;
    az = (Tz+Dz)/mass-gTitan;
    M = T*sin(delta)*LTVC;
    qDot = M/Iy;
    deltaDot = (deltaCmd-delta)/tauTVC;
    mdot = -T/(Isp*g0);

    if mass + mdot*dt < mDry
        mdot = (mDry-mass)/dt;
    end

    k4 = [Vx; Vz; ax; az; q; qDot; deltaDot; mdot];

    state = state + dt*(k1+2*k2+2*k3+k4)/6;

    if state(8) < mDry
        state(8) = mDry;
    end

    if state(2) <= 0
        state(2) = 0;
        state(4) = 0;
        X(:,k) = state;
        U(:,k) = [T;delta;thetaCmd;Treq;V];
        break;
    end

    X(:,k) = state;
    U(:,k) = [T;delta;thetaCmd;Treq;V];
end

valid = find(X(8,:) ~= 0);

if isempty(valid)
    valid = 1:N;
end

tv = t(valid);
Xv = X(:,valid);
Uv = U(:,valid);

figure;
plot(tv,Xv(2,:),'LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Altitude (m)');
title('Titan Powered Descent');

figure;
plot(tv,Xv(4,:),'LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Vertical Velocity (m/s)');
title('Titan Vertical Velocity');

figure;
plot(Xv(1,:),Xv(2,:),'LineWidth',1.5);
grid on;
xlabel('Horizontal Position (m)');
ylabel('Altitude (m)');
title('Titan Descent Trajectory');

figure;
plot(tv,Uv(1,:),'LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Thrust (N)');
title('Titan Thrust');

figure;
plot(tv,rad2deg(Uv(2,:)),'LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('TVC Angle (deg)');
title('Titan TVC Gimbal');

figure;
plot(tv,Xv(8,:),'LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Mass (kg)');
title('Titan Propellant Mass Evolution');

disp('Titan GNC and TVC simulation complete.');

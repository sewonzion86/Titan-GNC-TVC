function out = Titan_GNC_TVC_V2(doPlots,seed,opts)
% 6-DOF Titan powered descent: 2-axis TVC, ZEM/ZEV guidance, noisy sensors, optional engine-out.
%   out = Titan_GNC_TVC_V2(doPlots, seed, opts)    opts = struct overriding fields of params()
%   Titan_GNC_TVC_V2(false,12,struct('windScale',3,'engFailIdx',2))
% World: x east, y north, z up. Body x = thrust axis. Quaternion scalar-first, body->world.
% State: pos(3) vel(3) quat(4) rate(3) mass pitchGimbal yawGimbal thrust.
% Vehicle values are assumed (notional). Not flight software.

if nargin < 1, doPlots = true; end
if nargin < 2, seed = 7; end
if nargin < 3, opts = struct(); end
P = params();
fn = fieldnames(opts);
for i = 1:numel(fn)
    if ~isfield(P,fn{i}), error('Titan_GNC_TVC_V2:badOption','Unknown option: %s',fn{i}); end
    P.(fn{i}) = opts.(fn{i});
end
if P.m0*P.massScale <= P.mDry, error('Titan_GNC_TVC_V2:badMass','Wet mass must exceed dry mass.'); end
P.atm = buildAtmosphere(P);
if exist('OCTAVE_VERSION','builtin') ~= 0, randn('state',seed); else, rng(seed); end

dt = P.dt;  ctrlEvery = round(1/(P.ctrlHz*dt));
nSteps = round(P.tEnd/dt);
nTick  = ceil(nSteps/ctrlEvery) + 2;

qUp = [cos(-pi/4); 0; sin(-pi/4); 0];
a0 = P.tilt0;
q0 = quatMul(qUp,[cos(a0/2); 0; sin(a0/2); 0]);
X = [P.x0(1);P.x0(2);P.h0; P.v0; q0; P.w0; P.m0*P.massScale; 0;0;0];
C = initController(X,P);
cmd = [0;0;0;0];

L.t = zeros(1,nTick);  L.X = zeros(17,nTick);  L.qE = zeros(4,nTick);
L.pE = zeros(3,nTick); L.vE = zeros(3,nTick);  L.phase = zeros(1,nTick);
L.eB = zeros(3,nTick); L.cmd = zeros(4,nTick); L.qdyn = zeros(1,nTick);
nl = 0;  status = 'TIMEOUT';  tEndSim = P.tEnd;

for k = 0:nSteps-1
    t = k*dt;
    if mod(k,ctrlEvery) == 0
        S = sensorSample(t,X,k/ctrlEvery,P);
        [C,cmd] = flightSoftware(t,S,C,P);
        nl = nl + 1;
        L.t(nl) = t;  L.X(:,nl) = X;  L.qE(:,nl) = C.qE;  L.pE(:,nl) = C.pE;
        L.vE(:,nl) = C.vE;  L.phase(nl) = C.phase;  L.eB(:,nl) = C.eB;  L.cmd(:,nl) = cmd;
        [~,~,rho] = atmLookup(P.atm,max(X(3),0));
        L.qdyn(nl) = 0.5*rho*P.rhoScale*norm(X(4:6)-windAt(X(3),t,P))^2;
    end
    Xp = X;
    f1 = dyn(t,X,cmd,P);
    f2 = dyn(t+dt/2,X+dt*f1/2,cmd,P);
    f3 = dyn(t+dt/2,X+dt*f2/2,cmd,P);
    f4 = dyn(t+dt,X+dt*f3,cmd,P);
    X = X + dt*(f1+2*f2+2*f3+f4)/6;
    X(7:10) = X(7:10)/norm(X(7:10));
    X(14) = max(X(14),P.mDry);
    X(17) = max(X(17),0);
    X(15:16) = max(min(X(15:16),P.gMax),-P.gMax);
    if X(3) <= 0
        a = Xp(3)/(Xp(3)-X(3));
        X = Xp + a*(X-Xp);
        X(7:10) = X(7:10)/norm(X(7:10));
        status = 'TOUCHDOWN';  tEndSim = t + a*dt;  break;
    end
end

nl = min(nl,nTick);
fn = fieldnames(L);
for i = 1:numel(fn), L.(fn{i}) = L.(fn{i})(:,1:nl); end

R = quat2R(X(7:10));
iBurn = find(L.phase == 1,1);
out.status = status;
out.tFlight = tEndSim;
out.touchSpeed = norm(X(4:6));
out.touchVz = X(6);
out.touchLatSpeed = norm(X(4:5));
out.touchTiltDeg = acos(max(-1,min(1,R(3,1))))*180/pi;
out.landErr = norm(X(1:2)-P.target(:));
out.propUsed = P.m0*P.massScale - X(14);
out.propLeft = X(14) - P.mDry;
out.navErrEnd = norm(L.pE(:,end) - L.X(1:3,end));
if isempty(iBurn)
    out.ignitionAlt = NaN; out.ignitionSpeed = NaN; out.burnTime = 0;
else
    out.ignitionAlt = L.X(3,iBurn); out.ignitionSpeed = -L.X(6,iBurn);
    out.burnTime = tEndSim - L.t(iBurn);
end
out.peakThrust = max(L.X(17,:));
out.peakQdyn = max(L.qdyn);
out.maxGimbalDeg = max(max(abs(L.X(15:16,:))))*180/pi;
out.safeLanding = strcmp(status,'TOUCHDOWN') && out.touchVz > -2.0 && out.touchVz < 0 && ...
                  out.touchLatSpeed < 1.0 && out.touchTiltDeg < 10 && out.landErr < 5.0 && ...
                  out.propLeft > 0;
out.log = L;  out.params = P;

fprintf('Titan powered-descent GNC/TVC simulation: %s\n',status);
fprintf('  Flight time         : %.1f s (burn %.1f s)\n',out.tFlight,out.burnTime);
fprintf('  Ignition            : %.0f m altitude at %.1f m/s\n',out.ignitionAlt,out.ignitionSpeed);
fprintf('  Touchdown           : vertical %.2f m/s, lateral %.2f m/s, tilt %.1f deg\n', ...
        out.touchVz,out.touchLatSpeed,out.touchTiltDeg);
fprintf('  Landing error       : %.2f m from target\n',out.landErr);
fprintf('  Propellant used     : %.1f kg (left %.1f kg)\n',out.propUsed,out.propLeft);
fprintf('  Peak thrust / gimbal: %.0f N / %.1f deg, peak dynamic pressure %.0f Pa\n', ...
        out.peakThrust,out.maxGimbalDeg,out.peakQdyn);
if out.safeLanding, fprintf('  RESULT              : SAFE LANDING\n');
else, fprintf('  RESULT              : NOT a safe landing\n'); end

if doPlots, makePlots(L,P,out); end
end

% =========================================================================
function P = params()
d = pi/180;
P.GM = 8.978e12;  P.Rt = 2574.73e3;  P.g0 = 9.80665;
P.p0 = 146700;  P.Meff = 0.028535;  P.Rgas = 8.314462;
P.Tpts = [0 93.65; 10 83.0; 44 70.43; 60 71.5];
P.rhoScale = 1.0;  P.windScale = 1.0;
P.m0 = 1200;  P.mDry = 1050;  P.massScale = 1.0;  P.Isp = 310;  P.IspScale = 1.0;
P.Tmax = 10000;  P.Tmin = 750;  P.TmaxScale = 1.0;  P.tauT = 0.15;
P.D = 1.07;  P.Aref = pi*(P.D/2)^2;  P.Cd = 0.8;  P.CNa = 1.5;  P.cpOff = 0.12;
P.Cmq = -0.8;  P.Clp = -0.5;
P.Iwet = diag([650 1200 1200]);  P.Idry = diag([600 1100 1100]);
P.armWet = 1.00;  P.armDry = 1.04;
P.tauServo = 0.08;  P.servoRate = 20*d;  P.gMax = 12*d;  P.misalign = [0.2;-0.15]*d;
P.rcsMax = 60;
P.nEng = 4;  P.engRadius = 0.35;  P.engFailIdx = 0;  P.engFailAlt = 100;  P.fdirDelay = 0.5;
P.h0 = 2000;  P.x0 = [-60;30];  P.v0 = [0.8;0.1;-7];  P.w0 = [0;0.005;-0.004];  P.tilt0 = 3*d;
P.target = [0;0];
P.windLow = [0.3;0.1];  P.windHigh = [1.0;0.2];  P.windH = 5000;
P.gustAmp = [0.3;0.2;0.1];  P.gustW = [0.5;0.8;0.3];
P.gyroSigma = 0.01*d;  P.gyroBias = [0.010;-0.008;0.006]*d;  P.gyroBiasResid = 0.002*d;
P.accSigma = 0.01;  P.accBias = [0.004;-0.003;0.005];
P.radarSigma0 = 0.3;  P.radarSigmaK = 0.005;  P.radarMax = 5000;
P.dopSigma0 = 0.03;  P.dopSigmaK = 0.002;
P.trnSigma = 3.0;  P.trnEvery = 50;  P.trnAlt = 1500;
P.navPosErr = [8;8;5];  P.navVelErr = 0.2;  P.navAttErr = 0.1*d;  P.massEstErr = 0.005;
P.kR1 = 5;  P.kR2 = 6.25;  P.kDop = 0.2;  P.kT1 = 1.0;  P.kT2 = 0.3;
P.dt = 0.01;  P.ctrlHz = 50;  P.tEnd = 300;
P.aBrake = 4.0;  P.hMargin = 400;  P.vTouch = 1.0;  P.tgoMin = 1.5;
P.tiltMax = 25*d;  P.tiltLow = 8*d;  P.hTilt = 100;  P.hFinal = 15;
P.vLatMax = 15;  P.kLat = 1.0;  P.kPos = 0.4;  P.latUse = 0.5;
P.wn = 1.2;  P.zeta = 0.9;  P.Kp = P.wn^2;  P.Kd = 2*P.zeta*P.wn;  P.Ki = 0.1;  P.intMax = 1.0;
end

% =========================================================================
function atm = buildAtmosphere(P)
dh = 50;  h = 0:dh:60000;
T = interp1(P.Tpts(:,1)*1e3,P.Tpts(:,2),h);
p = zeros(size(h));  p(1) = P.p0;
for i = 2:numel(h)
    hm = 0.5*(h(i-1)+h(i));
    g = P.GM/(P.Rt+hm)^2;
    Tm = 0.5*(T(i-1)+T(i));
    p(i) = p(i-1)*exp(-g*P.Meff/(P.Rgas*Tm)*dh);
end
atm.dh = dh;  atm.T = T;  atm.p = p;  atm.rho = p*P.Meff./(P.Rgas*T);
end

function [T,p,rho] = atmLookup(atm,h)
x = min(max(h,0),59999)/atm.dh;
i = floor(x)+1;  f = x-(i-1);
T = atm.T(i)*(1-f)+atm.T(i+1)*f;
p = atm.p(i)*(1-f)+atm.p(i+1)*f;
rho = atm.rho(i)*(1-f)+atm.rho(i+1)*f;
end

function g = gravityAt(h,P)
g = P.GM/(P.Rt+max(h,0))^2;
end

function w = windAt(h,t,P)
s = min(max(h/P.windH,0),1);
base = (P.windLow*(1-s) + P.windHigh*s)*P.windScale;
w = [base; 0] + P.gustAmp.*sin(P.gustW*t)*P.windScale;
end

function off = engineOffset(P)
ang = (45 + 90*(P.engFailIdx-1))*pi/180;
off = -P.engRadius*[cos(ang); sin(ang)]/(P.nEng-1);
end

function [frac,off] = engineState(h,P)
frac = 1;  off = [0;0];
if P.engFailIdx > 0 && h < P.engFailAlt
    frac = (P.nEng-1)/P.nEng;  off = engineOffset(P);
end
end

% =========================================================================
function [dX,fb] = dyn(t,X,cmd,P)
v = X(4:6);  q = X(7:10);  w = X(11:13);  m = X(14);
dP = X(15);  dY = X(16);  T = X(17);
R = quat2R(q);
[burn,I] = massProps(m,P);
arm = P.armWet + burn*(P.armDry-P.armWet);
fuel = m > P.mDry + 1e-6;
[frac,off] = engineState(X(3),P);
Teff = T*fuel*frac;

[~,~,rho] = atmLookup(P.atm,X(3));
rho = rho*P.rhoScale;
vrel = v - windAt(X(3),t,P);
V = norm(vrel);
qd = 0.5*rho*V^2;
vb = R'*vrel;
if V > 1e-3, ub = vb/V; else, ub = [0;0;0]; end
Fd = -qd*P.Cd*P.Aref*ub;
Fn = -qd*P.Aref*P.CNa*[0; ub(2); ub(3)];
Vs = max(V,1);
Mdamp = qd*P.Aref*P.D*(P.D/(2*Vs))*[P.Clp*w(1); P.Cmq*w(2); P.Cmq*w(3)];
Maero = cross([P.cpOff;0;0],Fn) + Mdamp;

dPa = dP + P.misalign(1);  dYa = dY + P.misalign(2);
tb = Teff*[cos(dPa)*cos(dYa); sin(dYa); sin(dPa)*cos(dYa)];
Mthr = cross([-arm;off],tb);
Mrcs = [cmd(4);0;0];

Fb = tb + Fd + Fn;
a = R*Fb/m + [0;0;-gravityAt(X(3),P)];
fb = Fb/m;
wd = I\(Maero + Mthr + Mrcs - cross(w,I*w));
Om = [0 -w(1) -w(2) -w(3); w(1) 0 w(3) -w(2); w(2) -w(3) 0 w(1); w(3) w(2) -w(1) 0];
qdot = 0.5*Om*q;
mdot = -Teff/(P.Isp*P.IspScale*P.g0);
dPdot = max(min((cmd(2)-dP)/P.tauServo,P.servoRate),-P.servoRate);
dYdot = max(min((cmd(3)-dY)/P.tauServo,P.servoRate),-P.servoRate);
dTdot = (cmd(1)-T)/P.tauT;
dX = [v; a; qdot; wd; mdot; dPdot; dYdot; dTdot];
end

function [burn,I] = massProps(m,P)
mWet = P.m0*P.massScale;
burn = min(max((mWet-m)/(mWet-P.mDry),0),1);
I = P.Iwet + burn*(P.Idry-P.Iwet);
end

% =========================================================================
function S = sensorSample(t,X,tick,P)
[~,fb] = dyn(t,X,[0;0;0;0],P);
R = quat2R(X(7:10));
h = X(3);
S.gyro = X(11:13) + P.gyroBias + P.gyroSigma*randn(3,1);
S.acc  = fb + P.accBias + P.accSigma*randn(3,1);
slant = h/max(R(3,1),0.3);
S.radarValid = (h < P.radarMax);
S.radar = slant + (P.radarSigma0 + P.radarSigmaK*slant)*randn;
vb = R'*X(4:6);
S.dopp = vb + (P.dopSigma0 + P.dopSigmaK*norm(vb))*randn(3,1);
S.trnNew = (mod(tick,P.trnEvery) == 0) && (h < P.trnAlt);
S.trnPos = X(1:2) + P.trnSigma*randn(2,1);
end

% =========================================================================
function C = initController(X,P)
C.phase = 0;  C.tPhase = 0;
dq = P.navAttErr*randn(3,1);
C.qE = quatMul(X(7:10),[1;dq/2]);
C.pE = X(1:3) + P.navPosErr.*randn(3,1);
C.vE = X(4:6) + P.navVelErr*randn(3,1);
C.bE = P.gyroBias + P.gyroBiasResid*randn(3,1);
C.intY = 0;  C.intZ = 0;  C.Tf = 0;  C.Tprev = 0;
C.mE = max(X(14)*(1 + P.massEstErr*randn), P.mDry);
C.eB = zeros(3,1);
C.engOut = false;  C.tFail = -1;
end

function [C,cmd] = flightSoftware(t,S,C,P)
dtc = 1/P.ctrlHz;
cmd = [0;0;0;0];
fE = 1;  if C.engOut, fE = (P.nEng-1)/P.nEng; end

wE = S.gyro - C.bE;
C.qE = quatIntegrate(C.qE,wE,dtc);
RE = quat2R(C.qE);
aW = RE*S.acc + [0;0;-gravityAt(C.pE(3),P)];
C.pE = C.pE + C.vE*dtc + 0.5*aW*dtc^2;
C.vE = C.vE + aW*dtc;
if S.radarValid
    r = S.radar*RE(3,1) - C.pE(3);
    C.pE(3) = C.pE(3) + P.kR1*dtc*r;   C.vE(3) = C.vE(3) + P.kR2*dtc*r;
end
C.vE = C.vE + P.kDop*(RE*S.dopp - C.vE);
if S.trnNew
    Tg = P.trnEvery*dtc;
    rg = S.trnPos - C.pE(1:2);
    C.pE(1:2) = C.pE(1:2) + P.kT1*Tg*rg;
    C.vE(1:2) = C.vE(1:2) + P.kT2*Tg*rg;
end
C.Tf = C.Tf + (C.Tprev - C.Tf)*dtc/P.tauT;
C.mE = max(C.mE - fE*C.Tf/(P.Isp*P.g0)*dtc, P.mDry);
[burnE,IE] = massProps(C.mE,P);
armE = P.armWet + burnE*(P.armDry-P.armWet);

cmd(4) = max(min(-IE(1,1)*wE(1),P.rcsMax),-P.rcsMax);

hE = C.pE(3);  vE = C.vE;

if C.phase == 0
    hIgn = vE(3)^2/(2*P.aBrake) + P.hMargin;
    if vE(3) < 0 && hE <= hIgn
        C.phase = 1;  C.tPhase = 0;
    else
        C.eB = zeros(3,1);
        return
    end
end
C.tPhase = C.tPhase + dtc;

if P.engFailIdx > 0 && ~C.engOut && hE < P.engFailAlt
    if C.tFail < 0, C.tFail = t; end
    if t - C.tFail >= P.fdirDelay, C.engOut = true; end
end
fE = 1;  if C.engOut, fE = (P.nEng-1)/P.nEng; end

gE = [0;0;-gravityAt(hE,P)];
tgo = max(P.tgoMin, 2*max(hE,0.5)/(abs(vE(3)) + P.vTouch));
rf = [P.target(:); 0];  vf = [0;0;-P.vTouch];
if hE > P.hFinal
    ZEM = rf - (C.pE + vE*tgo + 0.5*gE*tgo^2);
    ZEV = vf - (vE + gE*tgo);
    ac = 6*ZEM/tgo^2 - 2*ZEV/tgo;
else
    ac = [0; 0; 2.0*(-P.vTouch - vE(3))] - gE;
end
tiltLim = P.tiltLow + (P.tiltMax-P.tiltLow)*min(max((hE-P.hFinal)/(P.hTilt-P.hFinal),0),1);
e = P.target(:) - C.pE(1:2);  dist = norm(e);
aLat = P.latUse*max(ac(3),-gE(3))*tan(P.tiltLow);
vdes = min([P.vLatMax, sqrt(2*aLat*dist), P.kPos*dist]);
if dist > 1e-6, vdir = e/dist; else, vdir = [0;0]; end
ac(1:2) = P.kLat*(vdes*vdir - vE(1:2));
[~,~,rhoE] = atmLookup(P.atm,max(hE,0));
Fd = -0.5*rhoE*norm(vE)*P.Cd*P.Aref*vE;
Fw = C.mE*ac - Fd;
Fw(3) = max(Fw(3),1.0);
horiz = norm(Fw(1:2));  maxH = Fw(3)*tan(tiltLim);
if horiz > maxH && horiz > 0, Fw(1:2) = Fw(1:2)*maxH/horiz; end
Tcmd = max(P.Tmin,min(P.Tmax*P.TmaxScale,norm(Fw)/fE));
dDes = Fw/norm(Fw);

aNose = RE*[1;0;0];
qe = qFromTwoVec(aNose,dDes);
eB = RE'*(2*qe(2:4));
C.eB = eB;

intY = max(min(C.intY + eB(2)*dtc, P.intMax),-P.intMax);
intZ = max(min(C.intZ + eB(3)*dtc, P.intMax),-P.intMax);
alY = P.Kp*eB(2) + P.Ki*intY - P.Kd*wE(2);
alZ = P.Kp*eB(3) + P.Ki*intZ - P.Kd*wE(3);
Teff = fE*max(C.Tf,P.Tmin);
ffP = 0;  ffY = 0;
if C.engOut
    off = engineOffset(P);
    ffP = -off(2)/armE;  ffY = -off(1)/armE;
end
dPraw =  IE(2,2)*alY/(armE*Teff) + ffP;
dYraw = -IE(3,3)*alZ/(armE*Teff) + ffY;
if abs(dPraw) < P.gMax, C.intY = intY; end
if abs(dYraw) < P.gMax, C.intZ = intZ; end
dPcmd = max(min(dPraw,P.gMax),-P.gMax);
dYcmd = max(min(dYraw,P.gMax),-P.gMax);

cmd(1:3) = [Tcmd; dPcmd; dYcmd];
C.Tprev = Tcmd;
end

% =========================================================================
function R = quat2R(q)
q = q/norm(q);  a = q(1); b = q(2); c = q(3); d = q(4);
R = [1-2*(c^2+d^2), 2*(b*c-a*d),   2*(b*d+a*c);
     2*(b*c+a*d),   1-2*(b^2+d^2), 2*(c*d-a*b);
     2*(b*d-a*c),   2*(c*d+a*b),   1-2*(b^2+c^2)];
end

function q = quatNorm(q)
q = q/norm(q);
end

function q = quatMul(a,b)
q = [a(1)*b(1)-dot(a(2:4),b(2:4));
     a(1)*b(2:4)+b(1)*a(2:4)+cross(a(2:4),b(2:4))];
q = q/norm(q);
end

function q = quatIntegrate(q,w,dt)
th = norm(w)*dt;
if th < 1e-12, dq = [1;0;0;0]; else, dq = [cos(th/2); sin(th/2)*w/norm(w)]; end
q = quatMul(q,dq);
end

function q = qFromTwoVec(a,b)
a = a/norm(a);  b = b/norm(b);
c = dot(a,b);
if c < -1 + 1e-9
    ax = cross(a,[0;0;1]);
    if norm(ax) < 1e-6, ax = cross(a,[0;1;0]); end
    q = [0; ax/norm(ax)];
else
    q = quatNorm([1+c; cross(a,b)]);
end
end

function c = tiltCos(Q)
c = zeros(1,size(Q,2));
for i = 1:size(Q,2)
    R = quat2R(Q(:,i));  c(i) = R(3,1);
end
end

% =========================================================================
function makePlots(L,P,out)
t = L.t;  X = L.X;  d = 180/pi;
 tilt = acos(max(-1,min(1,tiltCos(X(7:10,:)))))*d;
 tiltE = acos(max(-1,min(1,tiltCos(L.qE))))*d;
mach = zeros(size(t));  rho = mach;
for i = 1:numel(t)
    [Tk,~,rr] = atmLookup(P.atm,max(X(3,i),0));
    rho(i) = rr;  mach(i) = norm(X(4:6,i))/sqrt(1.4*P.Rgas*Tk/P.Meff);
end

figure('Name','Trajectory');
plot3(X(1,:),X(2,:),X(3,:),'LineWidth',1.5); hold on;
plot3(P.target(1),P.target(2),0,'rx','MarkerSize',10,'LineWidth',2);
grid on; axis equal; xlabel('East x (m)'); ylabel('North y (m)'); zlabel('Altitude (m)');
title(sprintf('Titan powered descent (touchdown %.2f m/s, %.1f m from target)',out.touchSpeed,out.landErr));

figure('Name','Dashboard','Position',[100 100 1300 850]);
subplot(3,3,1); plot(t,X(3,:),t,L.pE(3,:),'--','LineWidth',1.2); grid on;
xlabel('t (s)'); ylabel('Altitude (m)'); legend('true','estimated'); title('Altitude');
subplot(3,3,2); plot(t,X(4,:),t,X(5,:),t,X(6,:),'LineWidth',1.2); grid on;
xlabel('t (s)'); ylabel('Velocity (m/s)'); legend('vx','vy','vz'); title('Velocity');
subplot(3,3,3); plot(X(1,:),X(2,:),'LineWidth',1.3); hold on;
plot(P.target(1),P.target(2),'rx','LineWidth',2); grid on; axis equal;
xlabel('x (m)'); ylabel('y (m)'); title('Ground track');
subplot(3,3,4); plot(t,L.cmd(1,:),t,X(17,:),'LineWidth',1.2); grid on;
xlabel('t (s)'); ylabel('Thrust (N)'); legend('command','actual'); title('Thrust');
subplot(3,3,5); plot(t,L.cmd(2,:)*d,t,X(15,:)*d,'LineWidth',1.2); grid on;
xlabel('t (s)'); ylabel('Pitch TVC (deg)'); legend('command','actual'); title('Pitch gimbal');
subplot(3,3,6); plot(t,L.cmd(3,:)*d,t,X(16,:)*d,'LineWidth',1.2); grid on;
xlabel('t (s)'); ylabel('Yaw TVC (deg)'); legend('command','actual'); title('Yaw gimbal');
subplot(3,3,7); plot(t,tilt,t,tiltE,'--','LineWidth',1.2); grid on;
xlabel('t (s)'); ylabel('Tilt from vertical (deg)'); legend('true','estimated'); title('Attitude');
subplot(3,3,8); plot(t,vecnormCols(L.pE-X(1:3,:)),t,vecnormCols(L.vE-X(4:6,:)),'LineWidth',1.2); grid on;
xlabel('t (s)'); ylabel('Error (m, m/s)'); legend('position','velocity'); title('Navigation error');
subplot(3,3,9); plot(t,X(14,:),'LineWidth',1.3); grid on;
xlabel('t (s)'); ylabel('Mass (kg)'); title('Propellant burn');

figure('Name','Environment');
subplot(1,3,1); plot(rho,X(3,:),'LineWidth',1.3); grid on; xlabel('Density (kg/m^3)'); ylabel('Altitude (m)'); title('Titan density along path');
subplot(1,3,2); plot(t,L.qdyn,'LineWidth',1.3); grid on; xlabel('t (s)'); ylabel('Dynamic pressure (Pa)'); title('Dynamic pressure');
subplot(1,3,3); plot(t,mach,'LineWidth',1.3); grid on; xlabel('t (s)'); ylabel('Mach number'); title('Mach number');
end

function n = vecnormCols(A)
n = sqrt(sum(A.^2,1));
end

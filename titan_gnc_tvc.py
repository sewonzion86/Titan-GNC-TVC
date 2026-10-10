"""6-DOF Titan powered descent: 2-axis TVC, ZEM/ZEV guidance, noisy sensors, optional engine-out.

    python titan_gnc_tvc.py                  nominal run + plots
    python titan_gnc_tvc.py --noplot         nominal run, no plots
    python titan_gnc_tvc.py --cases          named stress cases
    python titan_gnc_tvc.py --mc 200         Monte Carlo (parallel) + requirement pass rates

Python API: out = simulate(seed=7, plots=False, windScale=3, engFailIdx=2)
World: x east, y north, z up. Body x = thrust axis. Quaternion scalar-first, body->world.
State: pos(3) vel(3) quat(4) rate(3) mass pitchGimbal yawGimbal thrust.
Port of Titan_GNC_TVC.m (same model; random streams differ, so results match statistically only).
Vehicle values are assumed (notional). Not flight software. Requires numpy; matplotlib for plots.
"""
import argparse
import math
from types import SimpleNamespace

import numpy as np

D2R = math.pi / 180.0


def params(**opts):
    d = D2R
    P = SimpleNamespace(
        GM=8.978e12, Rt=2574.73e3, g0=9.80665,
        p0=146700.0, Meff=0.028535, Rgas=8.314462,
        Tpts=np.array([[0, 93.65], [10, 83.0], [44, 70.43], [60, 71.5]]),
        rhoScale=1.0, windScale=1.0,
        m0=1200.0, mDry=1050.0, massScale=1.0, Isp=310.0, IspScale=1.0,
        Tmax=10000.0, Tmin=750.0, TmaxScale=1.0, tauT=0.15,
        D=1.07, Cd=0.8, CNa=1.5, cpOff=0.12, Cmq=-0.8, Clp=-0.5,
        Iwet=np.array([650.0, 1200.0, 1200.0]), Idry=np.array([600.0, 1100.0, 1100.0]),
        armWet=1.00, armDry=1.04,
        tauServo=0.08, servoRate=20 * d, gMax=12 * d, misalign=np.array([0.2, -0.15]) * d,
        rcsMax=60.0,
        nEng=4, engRadius=0.35, engFailIdx=0, engFailAlt=100.0, fdirDelay=0.5,
        h0=2000.0, x0=np.array([-60.0, 30.0]), v0=np.array([0.8, 0.1, -7.0]),
        w0=np.array([0.0, 0.005, -0.004]), tilt0=3 * d,
        target=np.array([0.0, 0.0]),
        windLow=np.array([0.3, 0.1]), windHigh=np.array([1.0, 0.2]), windH=5000.0,
        gustAmp=np.array([0.3, 0.2, 0.1]), gustW=np.array([0.5, 0.8, 0.3]),
        gyroSigma=0.01 * d, gyroBias=np.array([0.010, -0.008, 0.006]) * d, gyroBiasResid=0.002 * d,
        accSigma=0.01, accBias=np.array([0.004, -0.003, 0.005]),
        radarSigma0=0.3, radarSigmaK=0.005, radarMax=5000.0,
        dopSigma0=0.03, dopSigmaK=0.002,
        trnSigma=3.0, trnEvery=50, trnAlt=1500.0,
        navPosErr=np.array([8.0, 8.0, 5.0]), navVelErr=0.2, navAttErr=0.1 * d, massEstErr=0.005,
        kR1=5.0, kR2=6.25, kDop=0.2, kT1=1.0, kT2=0.3,
        dt=0.01, ctrlHz=50, tEnd=300.0,
        aBrake=4.0, hMargin=400.0, vTouch=1.0, tgoMin=1.5,
        tiltMax=25 * d, tiltLow=8 * d, hTilt=100.0, hFinal=15.0,
        vLatMax=15.0, kLat=1.0, kPos=0.4, latUse=0.5,
        wn=1.2, zeta=0.9, Ki=0.1, intMax=1.0,
    )
    P.Aref = math.pi * (P.D / 2) ** 2
    P.Kp = P.wn ** 2
    P.Kd = 2 * P.zeta * P.wn
    for k, v in opts.items():
        if not hasattr(P, k):
            raise KeyError(f"Unknown option: {k}")
        setattr(P, k, np.array(v, dtype=float) if isinstance(getattr(P, k), np.ndarray) else v)
    if P.m0 * P.massScale <= P.mDry:
        raise ValueError("Wet mass must exceed dry mass.")
    return P


def build_atmosphere(P):
    dh = 50.0
    h = np.arange(0.0, 60000.0 + dh, dh)
    T = np.interp(h, P.Tpts[:, 0] * 1e3, P.Tpts[:, 1])
    p = np.empty_like(h)
    p[0] = P.p0
    for i in range(1, len(h)):
        hm = 0.5 * (h[i - 1] + h[i])
        g = P.GM / (P.Rt + hm) ** 2
        Tm = 0.5 * (T[i - 1] + T[i])
        p[i] = p[i - 1] * math.exp(-g * P.Meff / (P.Rgas * Tm) * dh)
    return SimpleNamespace(dh=dh, T=T, p=p, rho=p * P.Meff / (P.Rgas * T))


def atm_lookup(atm, h):
    x = min(max(h, 0.0), 59999.0) / atm.dh
    i = int(x)
    f = x - i
    return (atm.T[i] * (1 - f) + atm.T[i + 1] * f,
            atm.p[i] * (1 - f) + atm.p[i + 1] * f,
            atm.rho[i] * (1 - f) + atm.rho[i + 1] * f)


def gravity(h, P):
    return P.GM / (P.Rt + max(h, 0.0)) ** 2


def wind_at(h, t, P):
    s = min(max(h / P.windH, 0.0), 1.0)
    base = (P.windLow * (1 - s) + P.windHigh * s) * P.windScale
    return np.append(base, 0.0) + P.gustAmp * np.sin(P.gustW * t) * P.windScale


def engine_offset(P):
    ang = (45 + 90 * (P.engFailIdx - 1)) * D2R
    return -P.engRadius * np.array([math.cos(ang), math.sin(ang)]) / (P.nEng - 1)


def engine_state(h, P):
    if P.engFailIdx > 0 and h < P.engFailAlt:
        return (P.nEng - 1) / P.nEng, engine_offset(P)
    return 1.0, np.zeros(2)


def cross3(a, b):
    return np.array([a[1] * b[2] - a[2] * b[1],
                     a[2] * b[0] - a[0] * b[2],
                     a[0] * b[1] - a[1] * b[0]])


def quat2R(q):
    q = q / np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[1 - 2 * (c * c + d * d), 2 * (b * c - a * d), 2 * (b * d + a * c)],
                     [2 * (b * c + a * d), 1 - 2 * (b * b + d * d), 2 * (c * d - a * b)],
                     [2 * (b * d - a * c), 2 * (c * d + a * b), 1 - 2 * (b * b + c * c)]])


def quat_mul(a, b):
    q = np.empty(4)
    q[0] = a[0] * b[0] - a[1:] @ b[1:]
    q[1:] = a[0] * b[1:] + b[0] * a[1:] + cross3(a[1:], b[1:])
    return q / np.linalg.norm(q)


def quat_integrate(q, w, dt):
    n = np.linalg.norm(w)
    th = n * dt
    if th < 1e-12:
        return q
    dq = np.concatenate(([math.cos(th / 2)], math.sin(th / 2) * w / n))
    return quat_mul(q, dq)


def q_from_two_vec(a, b):
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    c = float(a @ b)
    if c < -1 + 1e-9:
        ax = cross3(a, np.array([0.0, 0.0, 1.0]))
        if np.linalg.norm(ax) < 1e-6:
            ax = cross3(a, np.array([0.0, 1.0, 0.0]))
        return np.concatenate(([0.0], ax / np.linalg.norm(ax)))
    q = np.concatenate(([1 + c], cross3(a, b)))
    return q / np.linalg.norm(q)


def mass_props(m, P):
    mwet = P.m0 * P.massScale
    burn = min(max((mwet - m) / (mwet - P.mDry), 0.0), 1.0)
    return burn, P.Iwet + burn * (P.Idry - P.Iwet)


def dyn(t, X, cmd, P):
    v = X[3:6]
    q = X[6:10]
    w = X[10:13]
    m = X[13]
    dP, dY, T = X[14], X[15], X[16]
    R = quat2R(q)
    burn, I = mass_props(m, P)
    arm = P.armWet + burn * (P.armDry - P.armWet)
    fuel = 1.0 if m > P.mDry + 1e-6 else 0.0
    frac, off = engine_state(X[2], P)
    Teff = T * fuel * frac

    rho = atm_lookup(P.atm, X[2])[2] * P.rhoScale
    vrel = v - wind_at(X[2], t, P)
    V = float(np.linalg.norm(vrel))
    qd = 0.5 * rho * V * V
    vb = R.T @ vrel
    ub = vb / V if V > 1e-3 else np.zeros(3)
    Fd = -qd * P.Cd * P.Aref * ub
    Fn = -qd * P.Aref * P.CNa * np.array([0.0, ub[1], ub[2]])
    Vs = max(V, 1.0)
    Mdamp = qd * P.Aref * P.D * (P.D / (2 * Vs)) * np.array([P.Clp * w[0], P.Cmq * w[1], P.Cmq * w[2]])
    Maero = cross3(np.array([P.cpOff, 0.0, 0.0]), Fn) + Mdamp

    dPa = dP + P.misalign[0]
    dYa = dY + P.misalign[1]
    tb = Teff * np.array([math.cos(dPa) * math.cos(dYa), math.sin(dYa), math.sin(dPa) * math.cos(dYa)])
    Mthr = cross3(np.array([-arm, off[0], off[1]]), tb)
    Mrcs = np.array([cmd[3], 0.0, 0.0])

    Fb = tb + Fd + Fn
    a = R @ Fb / m + np.array([0.0, 0.0, -gravity(X[2], P)])
    fb = Fb / m
    wd = (Maero + Mthr + Mrcs - cross3(w, I * w)) / I
    Om = np.array([[0, -w[0], -w[1], -w[2]],
                   [w[0], 0, w[2], -w[1]],
                   [w[1], -w[2], 0, w[0]],
                   [w[2], w[1], -w[0], 0]])
    dX = np.empty(17)
    dX[0:3] = v
    dX[3:6] = a
    dX[6:10] = 0.5 * Om @ q
    dX[10:13] = wd
    dX[13] = -Teff / (P.Isp * P.IspScale * P.g0)
    dX[14] = min(max((cmd[1] - dP) / P.tauServo, -P.servoRate), P.servoRate)
    dX[15] = min(max((cmd[2] - dY) / P.tauServo, -P.servoRate), P.servoRate)
    dX[16] = (cmd[0] - T) / P.tauT
    return dX, fb


def sensor_sample(t, X, tick, P, rng):
    _, fb = dyn(t, X, np.zeros(4), P)
    R = quat2R(X[6:10])
    h = X[2]
    S = {}
    S["gyro"] = X[10:13] + P.gyroBias + P.gyroSigma * rng.standard_normal(3)
    S["acc"] = fb + P.accBias + P.accSigma * rng.standard_normal(3)
    slant = h / max(R[2, 0], 0.3)
    S["radarValid"] = h < P.radarMax
    S["radar"] = slant + (P.radarSigma0 + P.radarSigmaK * slant) * rng.standard_normal()
    vb = R.T @ X[3:6]
    S["dopp"] = vb + (P.dopSigma0 + P.dopSigmaK * np.linalg.norm(vb)) * rng.standard_normal(3)
    S["trnNew"] = (tick % P.trnEvery == 0) and (h < P.trnAlt)
    S["trnPos"] = X[0:2] + P.trnSigma * rng.standard_normal(2)
    return S


def init_controller(X, P, rng):
    C = SimpleNamespace()
    C.phase = 0
    C.tPhase = 0.0
    dq = P.navAttErr * rng.standard_normal(3)
    C.qE = quat_mul(X[6:10], np.concatenate(([1.0], dq / 2)))
    C.pE = X[0:3] + P.navPosErr * rng.standard_normal(3)
    C.vE = X[3:6] + P.navVelErr * rng.standard_normal(3)
    C.bE = P.gyroBias + P.gyroBiasResid * rng.standard_normal(3)
    C.intY = 0.0
    C.intZ = 0.0
    C.Tf = 0.0
    C.Tprev = 0.0
    C.mE = max(X[13] * (1 + P.massEstErr * rng.standard_normal()), P.mDry)
    C.eB = np.zeros(3)
    C.engOut = False
    C.tFail = -1.0
    return C


def flight_software(t, S, C, P):
    dtc = 1.0 / P.ctrlHz
    cmd = np.zeros(4)
    fE = (P.nEng - 1) / P.nEng if C.engOut else 1.0

    wE = S["gyro"] - C.bE
    C.qE = quat_integrate(C.qE, wE, dtc)
    RE = quat2R(C.qE)
    aW = RE @ S["acc"] + np.array([0.0, 0.0, -gravity(C.pE[2], P)])
    C.pE = C.pE + C.vE * dtc + 0.5 * aW * dtc ** 2
    C.vE = C.vE + aW * dtc
    if S["radarValid"]:
        r = S["radar"] * RE[2, 0] - C.pE[2]
        C.pE[2] += P.kR1 * dtc * r
        C.vE[2] += P.kR2 * dtc * r
    C.vE = C.vE + P.kDop * (RE @ S["dopp"] - C.vE)
    if S["trnNew"]:
        Tg = P.trnEvery * dtc
        rg = S["trnPos"] - C.pE[0:2]
        C.pE[0:2] += P.kT1 * Tg * rg
        C.vE[0:2] += P.kT2 * Tg * rg
    C.Tf += (C.Tprev - C.Tf) * dtc / P.tauT
    C.mE = max(C.mE - fE * C.Tf / (P.Isp * P.g0) * dtc, P.mDry)
    burnE, IE = mass_props(C.mE, P)
    armE = P.armWet + burnE * (P.armDry - P.armWet)

    cmd[3] = min(max(-IE[0] * wE[0], -P.rcsMax), P.rcsMax)

    hE = C.pE[2]
    vE = C.vE

    if C.phase == 0:
        hIgn = vE[2] ** 2 / (2 * P.aBrake) + P.hMargin
        if vE[2] < 0 and hE <= hIgn:
            C.phase = 1
            C.tPhase = 0.0
        else:
            C.eB = np.zeros(3)
            return C, cmd
    C.tPhase += dtc

    if P.engFailIdx > 0 and not C.engOut and hE < P.engFailAlt:
        if C.tFail < 0:
            C.tFail = t
        if t - C.tFail >= P.fdirDelay:
            C.engOut = True
    fE = (P.nEng - 1) / P.nEng if C.engOut else 1.0

    gE = np.array([0.0, 0.0, -gravity(hE, P)])
    tgo = max(P.tgoMin, 2 * max(hE, 0.5) / (abs(vE[2]) + P.vTouch))
    rf = np.array([P.target[0], P.target[1], 0.0])
    vf = np.array([0.0, 0.0, -P.vTouch])
    if hE > P.hFinal:
        ZEM = rf - (C.pE + vE * tgo + 0.5 * gE * tgo ** 2)
        ZEV = vf - (vE + gE * tgo)
        ac = 6 * ZEM / tgo ** 2 - 2 * ZEV / tgo
    else:
        ac = np.array([0.0, 0.0, 2.0 * (-P.vTouch - vE[2])]) - gE
    tilt_lim = P.tiltLow + (P.tiltMax - P.tiltLow) * min(max((hE - P.hFinal) / (P.hTilt - P.hFinal), 0.0), 1.0)
    e = P.target - C.pE[0:2]
    dist = float(np.linalg.norm(e))
    aLat = P.latUse * max(ac[2], -gE[2]) * math.tan(P.tiltLow)
    vdes = min(P.vLatMax, math.sqrt(2 * aLat * dist), P.kPos * dist)
    vdir = e / dist if dist > 1e-6 else np.zeros(2)
    ac[0:2] = P.kLat * (vdes * vdir - vE[0:2])
    rhoE = atm_lookup(P.atm, max(hE, 0.0))[2]
    Fd = -0.5 * rhoE * np.linalg.norm(vE) * P.Cd * P.Aref * vE
    Fw = C.mE * ac - Fd
    Fw[2] = max(Fw[2], 1.0)
    horiz = float(np.linalg.norm(Fw[0:2]))
    maxH = Fw[2] * math.tan(tilt_lim)
    if horiz > maxH and horiz > 0:
        Fw[0:2] *= maxH / horiz
    Tcmd = max(P.Tmin, min(P.Tmax * P.TmaxScale, float(np.linalg.norm(Fw)) / fE))
    dDes = Fw / np.linalg.norm(Fw)

    aNose = RE[:, 0]
    qe = q_from_two_vec(aNose, dDes)
    eB = RE.T @ (2 * qe[1:4])
    C.eB = eB

    intY = min(max(C.intY + eB[1] * dtc, -P.intMax), P.intMax)
    intZ = min(max(C.intZ + eB[2] * dtc, -P.intMax), P.intMax)
    alY = P.Kp * eB[1] + P.Ki * intY - P.Kd * wE[1]
    alZ = P.Kp * eB[2] + P.Ki * intZ - P.Kd * wE[2]
    Teff = fE * max(C.Tf, P.Tmin)
    ffP = ffY = 0.0
    if C.engOut:
        off = engine_offset(P)
        ffP = -off[1] / armE
        ffY = -off[0] / armE
    dPraw = IE[1] * alY / (armE * Teff) + ffP
    dYraw = -IE[2] * alZ / (armE * Teff) + ffY
    if abs(dPraw) < P.gMax:
        C.intY = intY
    if abs(dYraw) < P.gMax:
        C.intZ = intZ
    cmd[0] = Tcmd
    cmd[1] = min(max(dPraw, -P.gMax), P.gMax)
    cmd[2] = min(max(dYraw, -P.gMax), P.gMax)
    C.Tprev = Tcmd
    return C, cmd


def simulate(seed=7, plots=False, verbose=True, **opts):
    P = params(**opts)
    P.atm = build_atmosphere(P)
    rng = np.random.default_rng(seed)

    dt = P.dt
    ctrl_every = round(1 / (P.ctrlHz * dt))
    n_steps = round(P.tEnd / dt)

    q_up = np.array([math.cos(-math.pi / 4), 0.0, math.sin(-math.pi / 4), 0.0])
    a0 = P.tilt0
    q0 = quat_mul(q_up, np.array([math.cos(a0 / 2), 0.0, math.sin(a0 / 2), 0.0]))
    X = np.concatenate((P.x0, [P.h0], P.v0, q0, P.w0, [P.m0 * P.massScale, 0.0, 0.0, 0.0]))
    C = init_controller(X, P, rng)
    cmd = np.zeros(4)

    logs = {k: [] for k in ("t", "X", "qE", "pE", "vE", "phase", "cmd", "qdyn", "eB")}
    status, t_end = "TIMEOUT", P.tEnd

    for k in range(n_steps):
        t = k * dt
        if k % ctrl_every == 0:
            S = sensor_sample(t, X, k // ctrl_every, P, rng)
            C, cmd = flight_software(t, S, C, P)
            logs["t"].append(t)
            logs["X"].append(X.copy())
            logs["qE"].append(C.qE.copy())
            logs["pE"].append(C.pE.copy())
            logs["vE"].append(C.vE.copy())
            logs["phase"].append(C.phase)
            logs["cmd"].append(cmd.copy())
            logs["eB"].append(C.eB.copy())
            rho = atm_lookup(P.atm, max(X[2], 0.0))[2]
            vr = X[3:6] - wind_at(X[2], t, P)
            logs["qdyn"].append(0.5 * rho * P.rhoScale * float(vr @ vr))
        Xp = X
        f1, _ = dyn(t, X, cmd, P)
        f2, _ = dyn(t + dt / 2, X + dt * f1 / 2, cmd, P)
        f3, _ = dyn(t + dt / 2, X + dt * f2 / 2, cmd, P)
        f4, _ = dyn(t + dt, X + dt * f3, cmd, P)
        X = X + dt * (f1 + 2 * f2 + 2 * f3 + f4) / 6
        X[6:10] /= np.linalg.norm(X[6:10])
        X[13] = max(X[13], P.mDry)
        X[16] = max(X[16], 0.0)
        X[14:16] = np.clip(X[14:16], -P.gMax, P.gMax)
        if X[2] <= 0:
            a = Xp[2] / (Xp[2] - X[2])
            X = Xp + a * (X - Xp)
            X[6:10] /= np.linalg.norm(X[6:10])
            status, t_end = "TOUCHDOWN", t + a * dt
            break

    L = {k: np.array(v) for k, v in logs.items()}
    R = quat2R(X[6:10])
    burn_idx = np.flatnonzero(L["phase"] == 1)
    out = dict(status=status, tFlight=t_end)
    out["touchSpeed"] = float(np.linalg.norm(X[3:6]))
    out["touchVz"] = float(X[5])
    out["touchLatSpeed"] = float(np.linalg.norm(X[3:5]))
    out["touchTiltDeg"] = math.degrees(math.acos(max(-1.0, min(1.0, R[2, 0]))))
    out["landErr"] = float(np.linalg.norm(X[0:2] - P.target))
    out["propUsed"] = float(P.m0 * P.massScale - X[13])
    out["propLeft"] = float(X[13] - P.mDry)
    out["navErrEnd"] = float(np.linalg.norm(L["pE"][-1] - L["X"][-1, 0:3]))
    if burn_idx.size:
        i = burn_idx[0]
        out["ignitionAlt"] = float(L["X"][i, 2])
        out["ignitionSpeed"] = float(-L["X"][i, 5])
        out["burnTime"] = t_end - float(L["t"][i])
    else:
        out["ignitionAlt"] = out["ignitionSpeed"] = float("nan")
        out["burnTime"] = 0.0
    out["peakThrust"] = float(L["X"][:, 16].max())
    out["peakQdyn"] = float(L["qdyn"].max())
    out["maxGimbalDeg"] = float(np.degrees(np.abs(L["X"][:, 14:16]).max()))
    out["safeLanding"] = (status == "TOUCHDOWN" and -2.0 < out["touchVz"] < 0 and
                          out["touchLatSpeed"] < 1.0 and out["touchTiltDeg"] < 10 and
                          out["landErr"] < 5.0 and out["propLeft"] > 0)
    out["log"] = L
    out["params"] = P

    if verbose:
        print(f"Titan powered-descent GNC/TVC simulation: {status}")
        print(f"  Flight time         : {out['tFlight']:.1f} s (burn {out['burnTime']:.1f} s)")
        print(f"  Ignition            : {out['ignitionAlt']:.0f} m altitude at {out['ignitionSpeed']:.1f} m/s")
        print(f"  Touchdown           : vertical {out['touchVz']:.2f} m/s, lateral {out['touchLatSpeed']:.2f} m/s, "
              f"tilt {out['touchTiltDeg']:.1f} deg")
        print(f"  Landing error       : {out['landErr']:.2f} m from target")
        print(f"  Propellant used     : {out['propUsed']:.1f} kg (left {out['propLeft']:.1f} kg)")
        print(f"  Peak thrust / gimbal: {out['peakThrust']:.0f} N / {out['maxGimbalDeg']:.1f} deg, "
              f"peak dynamic pressure {out['peakQdyn']:.0f} Pa")
        print("  RESULT              : " + ("SAFE LANDING" if out["safeLanding"] else "NOT a safe landing"))
    if plots:
        make_plots(out)
    return out


def make_plots(out):
    import matplotlib.pyplot as plt
    L, P = out["log"], out["params"]
    t, X = L["t"], L["X"]
    r2d = 180 / math.pi
    tilt = np.degrees(np.arccos(np.clip([quat2R(q)[2, 0] for q in X[:, 6:10]], -1, 1)))
    tilt_e = np.degrees(np.arccos(np.clip([quat2R(q)[2, 0] for q in L["qE"]], -1, 1)))

    ax = plt.figure("Trajectory").add_subplot(projection="3d")
    ax.plot(X[:, 0], X[:, 1], X[:, 2])
    ax.plot([P.target[0]], [P.target[1]], [0], "rx", ms=10)
    ax.set_xlabel("East x (m)")
    ax.set_ylabel("North y (m)")
    ax.set_zlabel("Altitude (m)")
    ax.set_title(f"Titan powered descent ({out['touchSpeed']:.2f} m/s, {out['landErr']:.1f} m from target)")

    fig, axs = plt.subplots(3, 3, figsize=(14, 9), num="Dashboard")
    a = axs.ravel()
    a[0].plot(t, X[:, 2], label="true")
    a[0].plot(t, L["pE"][:, 2], "--", label="estimated")
    a[0].set(title="Altitude", xlabel="t (s)", ylabel="m")
    a[1].plot(t, X[:, 3], label="vx")
    a[1].plot(t, X[:, 4], label="vy")
    a[1].plot(t, X[:, 5], label="vz")
    a[1].set(title="Velocity", xlabel="t (s)", ylabel="m/s")
    a[2].plot(X[:, 0], X[:, 1])
    a[2].plot(*P.target, "rx")
    a[2].set(title="Ground track", xlabel="x (m)", ylabel="y (m)", aspect="equal")
    a[3].plot(t, L["cmd"][:, 0], label="command")
    a[3].plot(t, X[:, 16], label="actual")
    a[3].set(title="Thrust", xlabel="t (s)", ylabel="N")
    a[4].plot(t, L["cmd"][:, 1] * r2d, label="command")
    a[4].plot(t, X[:, 14] * r2d, label="actual")
    a[4].set(title="Pitch gimbal", xlabel="t (s)", ylabel="deg")
    a[5].plot(t, L["cmd"][:, 2] * r2d, label="command")
    a[5].plot(t, X[:, 15] * r2d, label="actual")
    a[5].set(title="Yaw gimbal", xlabel="t (s)", ylabel="deg")
    a[6].plot(t, tilt, label="true")
    a[6].plot(t, tilt_e, "--", label="estimated")
    a[6].set(title="Attitude tilt", xlabel="t (s)", ylabel="deg")
    a[7].plot(t, np.linalg.norm(L["pE"] - X[:, 0:3], axis=1), label="position")
    a[7].plot(t, np.linalg.norm(L["vE"] - X[:, 3:6], axis=1), label="velocity")
    a[7].set(title="Navigation error", xlabel="t (s)")
    a[8].plot(t, X[:, 13])
    a[8].set(title="Mass", xlabel="t (s)", ylabel="kg")
    for i, ax_ in enumerate(a):
        ax_.grid(True)
        if i != 2 and i != 8:
            ax_.legend()
    fig.tight_layout()
    plt.show()


REQ = [
    ("R1 touchdown vz in (-2, 0) m/s", lambda o: -2.0 < o["touchVz"] < 0),
    ("R2 lateral speed < 1 m/s", lambda o: o["touchLatSpeed"] < 1.0),
    ("R3 tilt < 10 deg", lambda o: o["touchTiltDeg"] < 10),
    ("R4 landing error < 5 m", lambda o: o["landErr"] < 5.0),
    ("R5 propellant reserve >= 10 kg", lambda o: o["propLeft"] >= 10),
    ("R6 peak gimbal < 10.8 deg (90% authority)", lambda o: o["maxGimbalDeg"] < 10.8),
    ("R7 terminated by touchdown", lambda o: o["status"] == "TOUCHDOWN"),
    ("R8 ignition altitude > 100 m", lambda o: o["ignitionAlt"] > 100),
    ("R9 final nav position error < 5 m", lambda o: o["navErrEnd"] < 5.0),
]

CASES = [
    ("nominal", {}),
    ("engine out (engine 2, at 100 m)", dict(engFailIdx=2)),
    ("dense air +10%", dict(rhoScale=1.10)),
    ("thin air -10%", dict(rhoScale=0.90)),
    ("wind x3", dict(windScale=3.0)),
    ("Isp -3%", dict(IspScale=0.97)),
    ("mass +3%", dict(massScale=1.03)),
    ("TVC misalign x3", dict(misalign=3 * np.array([0.2, -0.15]) * D2R)),
    ("gyro bias x3", dict(gyroBias=3 * np.array([0.010, -0.008, 0.006]) * D2R)),
    ("offset start 150 m", dict(x0=np.array([-150.0, 100.0]))),
]


def _check(o):
    return [bool(f(o)) for _, f in REQ]


def _case_job(item):
    name, op = item
    try:
        o = simulate(seed=7, verbose=False, **op)
        return name, o
    except Exception as e:
        return name, e


def _mc_job(args):
    k, = args
    rng = np.random.default_rng(5000 + k)
    u = lambda a, b: rng.uniform(a, b)
    g = lambda n=None: rng.standard_normal(n)
    d = D2R
    op = dict(
        windScale=u(0.5, 3.0), rhoScale=u(0.90, 1.10), massScale=u(0.97, 1.03),
        IspScale=u(0.97, 1.03), TmaxScale=u(0.90, 1.00),
        x0=np.array([-60.0, 30.0]) + 40 * g(2),
        v0=np.array([0.8, 0.1, -7.0]) + 0.5 * g(3),
        tilt0=abs(3 * d + 1.5 * d * g()),
        misalign=np.array([0.2, -0.15]) * d + 0.1 * d * g(2),
        gyroBias=np.array([0.010, -0.008, 0.006]) * d * (1 + 0.5 * g(3)),
        accBias=np.array([0.004, -0.003, 0.005]) + 0.002 * g(3),
    )
    try:
        o = simulate(seed=1000 + k, verbose=False, **op)
        return _check(o), o
    except Exception:
        return [False] * len(REQ), None


def run_cases(workers=None):
    from concurrent.futures import ProcessPoolExecutor
    print("=== STRESS CASES ===")
    with ProcessPoolExecutor(workers) as ex:
        results = list(ex.map(_case_job, CASES))
    for name, o in results:
        if isinstance(o, Exception):
            print(f"{name:34s} ERROR: {o}")
            continue
        ok = all(_check(o))
        print(f"{name:34s} {'PASS' if ok else 'FAIL'}  vz {o['touchVz']:.2f}  lat {o['touchLatSpeed']:.2f}  "
              f"err {o['landErr']:.2f} m  prop left {o['propLeft']:.1f} kg")


def wilson(k, n, z=1.96):
    p = k / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, c - h), min(1.0, c + h)


def monte_carlo(n=100, workers=None):
    from concurrent.futures import ProcessPoolExecutor
    print(f"=== MONTE CARLO (N = {n}) ===")
    with ProcessPoolExecutor(workers) as ex:
        results = list(ex.map(_mc_job, [(k,) for k in range(n)]))
    passes = np.array([r[0] for r in results])
    all_ok = passes.all(axis=1)
    lo, hi = wilson(int(all_ok.sum()), n)
    print(f"Overall pass: {all_ok.sum()}/{n} = {100 * all_ok.mean():.1f} %  (95% Wilson CI {100 * lo:.1f} - {100 * hi:.1f} %)")
    print("\nPer-requirement pass rate:")
    for (name, _), col in zip(REQ, passes.T):
        print(f"  {name:48s} {100 * col.mean():5.1f} %")
    keys = ["touchVz", "touchLatSpeed", "touchTiltDeg", "landErr", "propLeft",
            "maxGimbalDeg", "ignitionAlt", "burnTime", "peakQdyn", "navErrEnd"]
    outs = [r[1] for r in results if r[1] is not None]
    print(f"\n{'metric':16s} {'mean':>9s} {'std':>9s} {'p1':>9s} {'p99':>9s}")
    for key in keys:
        x = np.array([o[key] for o in outs], dtype=float)
        x = x[~np.isnan(x)]
        print(f"{key:16s} {x.mean():9.3f} {x.std():9.3f} {np.percentile(x, 1):9.3f} {np.percentile(x, 99):9.3f}")
    return passes


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--noplot", action="store_true")
    ap.add_argument("--cases", action="store_true")
    ap.add_argument("--mc", type=int, default=0)
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()
    if args.cases:
        run_cases(args.workers)
    elif args.mc:
        monte_carlo(args.mc, args.workers)
    else:
        simulate(plots=not args.noplot)

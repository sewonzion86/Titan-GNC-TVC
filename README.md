# Titan-GNC-TVC

MATLAB simulation of a powered descent lander for Titan with GNC and TVC control.

## Overview

This repository models a powered descent and landing scenario for Titan using a simplified guidance, navigation, and control architecture. The simulation includes altitude tracking, thrust modulation, TVC gimbal control, sensor noise, and Monte Carlo verification. It is intended as a notional Titan powered-descent concept model, not flight software.

## Repo status

This repository includes the original scripts plus newer update / v2 variants:

- `Titan_GNC_TVC.m` — original MATLAB implementation
- `Titan_GNC_TVC (1).m` — updated MATLAB v2 script
- `Titan_Verify_MC.m` — Monte Carlo verification script for the updated model
- `titan_gnc_tvc.py` — Python port / v2 update script
- `launcher.m` — convenience runner for the original MATLAB simulation

## MATLAB Files

- `Titan_GNC_TVC.m` — main original simulation script
- `Titan_GNC_TVC (1).m` — v2 / updated simulation model
- `Titan_Verify_MC.m` — verification and Monte Carlo analysis for the updated model
- `launcher.m` — convenience runner that executes the main simulation

## Python Files

- `titan_gnc_tvc.py` — Python equivalent / update script of the Titan GNC/TVC model

## Simulation Output Figures

### Powered descent altitude

![Titan Lander powered descent altitude](figures/titan_lander_powered_descent_altitude.svg)

### Thrust command

![Titan Lander thrust command](figures/titan_lander_thrust_command.svg)

### TVC gimbal

![Titan Lander TVC gimbal](figures/titan_lander_tvc_gimbal.svg)

### Vertical velocity

![Titan Lander vertical velocity](figures/titan_lander_vertical_velocity.svg)

### Conceptual 3D descent path

![Titan Lander conceptual 3D descent path](figures/titan_lander_3d_descent_path.svg)

## Repository Structure

```text
Titan-GNC-TVC/
├── README.md
├── Titan_GNC_TVC.m
├── Titan_GNC_TVC (1).m
├── Titan_Verify_MC.m
├── titan_gnc_tvc.py
├── launcher.m
├── figures/
│   ├── titan_lander_powered_descent_altitude.svg
│   ├── titan_lander_thrust_command.svg
│   ├── titan_lander_tvc_gimbal.svg
│   ├── titan_lander_vertical_velocity.svg
│   └── titan_lander_3d_descent_path.svg
└── main simulation files
```

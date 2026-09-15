# CHRONOS — Nonlinear MPC through LPV models

**Define and tune your controller in MATLAB, then generate C/C++ for embedded deployment.**

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=arielmb94/CHRONOS-MPC)

Changing operating conditions, coupled dynamics, and actuator limits make demanding control problems a natural fit for Model Predictive Control (MPC). CHRONOS helps you build predictive controllers for nonlinear and time-varying systems: define your prediction model, tracking objectives, and operating constraints through a MATLAB API, and CHRONOS assembles and solves the optimization problem.

CHRONOS uses Linear Parameter Varying (LPV) models to bring nonlinear control into a structured convex optimization framework. You supply the prediction-model matrices and update them as operating conditions change; CHRONOS uses the supplied matrices for each solve. This gives you a practical nonlinear MPC workflow without writing solver code or deriving gradients and Hessians for the supported costs and constraints.

---

## Table of Contents

1. [Why CHRONOS?](#why-chronos)
2. [Key Benefits](#key-benefits)
3. [Getting Started](#getting-started)
4. [API Workflow Overview](#api-workflow-overview)
5. [How It Works](#how-it-works)
6. [Contact](#contact)
7. [Citation](#citation)

---

## Why CHRONOS?

- **Robust numerical methods:** Built for demanding control applications, CHRONOS combines well established robust convex optimization methods with a tailored Riccati solver. For fixed state and input dimensions, the solve time scales linearly with the prediction horizon.
- **Focus on control, not optimization:** Put your engineering effort into the controller: choose the model, objectives, and constraints, and let CHRONOS handle optimization assembly, gradients, Hessians, and the numerical solve.
- **MATLAB → C/C++:** Take your controller from MATLAB development to embedded execution. Tune it in MATLAB, then generate C/C++ with MATLAB Coder to bring your design to fast control loops.

---

## Key Benefits

- **MPC-LPV paradigm:** The LPV (or Linear Time-Varying) framework combines exact representations of suitable nonlinear dynamics with the freedom to incorporate partial linearization or hybrid behaviour, while keeping each MPC solve within a structured convex formulation.
- **Custom behaviours:** Shape the controller around your application. Custom linear and quadratic costs and linear constraints let you combine model signals, balance competing objectives, and build advanced control behaviours beyond standard tracking.
- **API-driven:** Express your controller in familiar control-engineering terms. Simple API functions make reference tracking, actuator limits, control-effort and control-rate penalties, and terminal costs easy to configure and tune, with dedicated functions for online updates.
- **Purpose-built for MPC performance:** CHRONOS exploits the specific structure of your dynamics, costs, and constraints throughout the solve. Its integrated problem definition and numerical core enable specialized operations that cut unnecessary operations and arithmetic.
- **Open and adaptable:** Fully accessible, MIT-licensed MATLAB source lets you inspect, debug, and extend CHRONOS to meet your application's needs.

---

## Getting Started

1. Clone the repository.
2. Add CHRONOS to your MATLAB path.
3. Browse the Tutorials for step‑by‑step guides.
4. Try one of the Examples to see CHRONOS in action.

## API Workflow Overview

**Define once → build once → update model/data → solve → apply the first control action.**

The snippet below illustrates the API workflow. For a complete application, see the [Two Tank example](examples/Two%20Tank/readme.md).

```matlab
% Define once: prediction model, constraints, and tracking objective
mpc = init_mpc(N);
mpc = init_mpc_dynamics(mpc, A, B, []);
mpc = init_mpc_output(mpc, C, D, []);
mpc = init_mpc_control_rate_cnstr(mpc, du_min, du_max);
mpc = init_mpc_Tracking_cost(mpc, Qe);

% Build once: initialize the fixed-size solver workspace
mpc = build_chronos_mpc(mpc, s_prev, u_prev, [], []);

% At each control sample: supply updated model matrices and current inputs
mpc = update_mpc_dynamics(mpc, A, B, []);

% Compute the next control action
[u0, mpc, iter] = mpc_solve(mpc, s_prev, u_prev, r_in, [], [], [], []);

% Apply u0 to the plant and retain mpc for the next control sample
```

---

## How It Works



### MPC problem formulation and custom extensions

CHRONOS solves a finite-horizon MPC problem with ingredients:

$$
\begin{aligned}
\min_{x,u} J= \quad & \frac{1}{2}\sum_{k=0}^{N}(r_k-y_k)^TQ_{e,k}(r_k-y_k) \\
& + \sum_{k=0}^{N-1}\left(\frac{1}{2}u_k^TR_{u,k}u_k+r_{u,k}^Tu_k\right) \\
& + \frac{1}{2}\sum_{k=0}^{N-1}\Delta u_k^TR_{du,k}\Delta u_k \\
& + (x_{N,\mathrm{ref}}-x_N)^TP(x_{N,\mathrm{ref}}-x_N) \\
{} & {} \\
\text{s.t.} \quad & x_{k+1}=A_kx_k+B_ku_k+B_{d,k}d_k,\quad k=0,\dots,N-1 \\
{} & {} \\
& x_{\min}(k)\leq x(k)\leq x_{\max}(k),\quad k=1,\dots,N \\
& u_{\min}(k)\leq u(k)\leq u_{\max}(k),\quad k=0,\dots,N-1 \\
& \Delta u_{\min}(k)\leq \Delta u(k)\leq \Delta u_{\max}(k),\quad k=0,\dots,N-1 \\
& y_{\min}(k)\leq y(k)\leq y_{\max}(k),\quad k=0,\dots,N
\end{aligned}
$$

$N$ is the prediction horizon, $y_k=C_kx_k+D_ku_k+D_{d,k}d_k$, and $\Delta u_k=u_k-u_{k-1}$.

Additionally, CHRONOS gives you the flexibility to add custom terms to the cost function. For this, you can define a user vector $z(k)$ as a time-varying linear combination of the state, control action, preceding control action, and known input disturbance.

$$
z_k = C_{z,k}x_k+D_{z,k}u_k+D_{su,z,k}u_{k-1}+D_{dz,k}dz_k.
$$
CHRONOS then allows you to add quadratic and linear costs to this vector:
$$
\begin{aligned}
J_{\mathrm{custom}} &= \sum_{k=0}^{N}\left(\frac{1}{2}z_k^TQ_{z,k}z_k+q_{z,k}^Tz_k\right), \\
\end{aligned}
$$


Similarly you can define custom constraints by defining a user vector $h(k)$ in the same way:

$$
\begin{aligned}
h_k &= C_{h,k}x_k+D_{h,k}u_k+D_{su,h,k}u_{k-1}+D_{dh,k}dh_k, \\
\end{aligned}
$$

and constrain it componentwise: 

$$
\begin{aligned}
h_{\min}(k)&\leq h(k)\leq h_{\max}(k),\quad k=0,\dots,N
\end{aligned}
$$

---

## Contact

If you are interested in using CHRONOS in a professional or industrial setting and need help in the process, or wish to discuss potential collaborations, please feel free to get in touch via [LinkedIn](https://www.linkedin.com/in/ariel-medero-borrell). We welcome technical inquiries and are open to exploring tailored applications.



---

## Citation

If you use CHRONOS in your academic work, please cite:

M. Borrell, A. (2025). CHRONOS: solver for receding horizon control  of parameter varying convex systems, Online:  https://github.com/arielmb94/CHRONOS-MPC

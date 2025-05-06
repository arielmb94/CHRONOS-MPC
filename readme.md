# CHRONOS: solver for reCeding Horizon contROl  of parameter varyiNg cOnvex Systems 

CHRONOS is a Model Predictive Control (MPC) solver tailored for Linear Parameter Varying (LPV) systems. Built from the ground up in MATLAB, it leverages mature convex optimization techniques to deliver a robust, high‑performance control engine—while keeping optimization details completely under the hood. Thanks to the MPC‑LPV paradigm, CHRONOS brings the power of nonlinear MPC to real-time applications with the efficiency of convex solvers. You can tackle complex, nonlinear, time‑varying systems without ever writing a single line of optimization code.

---

## Table of Contents

1. [Why CHRONOS?](#why-chronos)
2. [Key Benefits](#key-benefits)
3. [Getting Started](#getting-started)
4. [Quick API Example](#quick-api-example)
5. [MPC Problem Formulation](#mpc-problem-formulation)
6. [Contact](#contact)
7. [Citation](#citation)

---

## Why CHRONOS?

- **Solid theory, simplified workflow**: Underneath the hood, CHRONOS uses a log‑barrier interior point method with tailored enhancements for MPC. You get convergence and robustness without writing a line of solver code yourself.
- **Focus on control, not optimization**: All gradient, Hessian and constraint assembly is automated. Define your model, cost and constraints—and CHRONOS handles the rest.
- **Seamless MATLAB→C/C++**: Write and tune your MPC in MATLAB, then generate production‑ready code in a few clicks.

---

## Key Benefits

- **MPC‑LPV paradigm**: Solve nonlinear MPC problems efficiently by leveraging LPV models. CHRONOS uses LPV models to obtain a structured convex optimization problem that’s fast, and realiable while being able to handle natively time‑varying dynamics.
- **Custom behaviours**: CHRONOS already implements the cost terms and constraints typically used in MPC, e.g. tracking error penalty, control action penalty, state constraints, constraints on the control change in between samples and more. In addition, CHRONOS allows to set user-defined penalty terms and inequalities, so that you can implement custom extensions to the classical MPC definition.
- **Speed tweaks**: Single‑layer interior‑point, computationally efficient feasibility line search, adjustable iteration limits, and warm‑start strategies—so you solve faster, with graceful approximation.
- **API‑driven**: Intuitive API functions allow you to define, initialize and update during runtime every parameter of your MPC problem—without worrying about gradients, Hessians, or solver internals. 
- **Code Generation**: CHRONOS has been completely programmed from scratch using only basic Matlab operators and functions. This allows to use Matlab Coder to export your tuned CHRONOS controller to C/C++, enabling deployment to embedded systems or high‑performance applications.

---

## Getting Started

1. Clone the repository.
2. Add CHRONOS to your MATLAB path.
3. Browse the Tutorials for step‑by‑step guides.
4. Try one of the Examples to see CHRONOS in action.

## Quick API Examples
 
 ```matlab
% Create mpc problem structure
 mpc = init_mpc(N, N_h_ctr);

% Initialize system dynamics
mpc = init_mpc_system(mpc, A, B, Bd, C, D, Dd);

% Control inputs variation constraints
 mpc = init_mpc_delta_u_cnstr(mpc, du_min, u_max);

% Tracking error penalty
mpc = init_mpc_Tracking_cost(mpc, Qe);

% Update dynamics model during runtime
mpc = update_mpc_sys_dynamics(mpc, A, B, []);

% lunch MPC iteration
[u, x0] = mpc_solve(mpc, x0, x, u_prev, ref , d, x_refN, dz, dh);
```
---

## MPC Problem Formulation

CHRONOS solves a finite-horizon MPC problem of the form:

$$
\begin{aligned}
\min_{x,u} J= \quad & (x_{refN}-x_N)^T P (x_{refN}-x_N) + \sum_{i=0}^{N-1} (r_i - y_i)^T Q_e (r_i - y_i) + \sum_{i=0}^{N_{\text{ctr}}-1} \left( u_i^T R_u u_i + r_u^T u_i \right) + \sum_{i=0}^{N_{\text{ctr}}-1}\Delta u_i^T dR_u \Delta u_i + \sum_{i=1}^{N-1} \left( z_i^T Q_z z_i + q_z^T z_i \right) \\
{} & {} \\
\text{s.t.} \quad & x_{i+1} = A x_i + B u_i + D_i d_i,\quad i = 0, \dots, N \\
{} & {} \\
& x_i \in [x_{\min}, x_{\max}],\quad i = 1, \dots, N+1 \\
& u_i \in [u_{\min}, u_{\max}],\quad i = 0, \dots, N_{\text{ctr}} - 1 \\
& \Delta u_i \in [\Delta u_{\min}, \Delta u_{\max}],\quad i = 0, \dots, N_{\text{ctr}} - 1 \\
& y_i \in [y_{\min}, y_{\max}],\quad i = 1, \dots, N \\
& h_i \in [h_{\min}, h_{\max}],\quad i = 1, \dots, N \\
& (x_{refN}-x_N)^T P (x_{refN}-x_N) \leq 1
\end{aligned}
$$

Notes:
* $N$ is the prediction horizon, $N_{ctr}$ the control horizon.
* $y$: tracking output
$$y = Cx+Du+D_dd$$ 
* $z$: custom cost signal
$$z = C_zx+D_zu+D_{dz}d_z$$
* $h$: custom constraint signal
$$h = C_hx+D_hu+D_{dh}d_z$$ 
* $d_i$, $d_z$ and $d_h$ are vectors of known disturbance signals on the tracking output and on the custom cost and constraint signals respectevely. These signals are passed to CHRONOS during runtime execution of the MPC.
* Terminal constraints and terminal cost built-in for stability and feasibility guarantees.

---

## Contact

If you are interested in using CHRONOS in a professional or industrial setting and need help in the process, or wish to discuss potential collaborations, please feel free to get in touch via [LinkedIn](https://www.linkedin.com/in/ariel-medero-borrell). We welcome technical inquiries and are open to exploring tailored applications.



---

## Citation

If you use CHRONOS in your academic work, please cite:

M. Borrell, A. (2025). CHRONOS: solver for receding horizon control  of parameter varying convex systems, Online:  https://github.com/arielmb94/CHRONOS-MPC


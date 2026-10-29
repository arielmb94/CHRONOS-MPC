# LPV-MPC Scheduling Trajectory Estimation in CHRONOS

In order to solve the sampled quasi-Linear Parameter-Varying (qLPV) MPC optimization problem in **CHRONOS**, an estimate for the future scheduling parameter trajectory along the prediction horizon $N_p$,

$$\hat{P}_k := \text{col}\Big(\rho(k \mid k),\, \rho(k+1 \mid k),\, \dots,\, \rho(k+N_p-1 \mid k)\Big) \in \mathbb{R}^{n_\rho N_p},$$

is required at each sampling instant $k$. Once $\hat{P}_k$ is defined, the state-space prediction matrices become known over the horizon, preserving a convex Quadratic Programming (QP) structure that CHRONOS solves using fast block-structured Riccati recursions.

Depending on your system's degree of nonlinearity and available computational budget, CHRONOS provides **four main approaches** to estimate $\hat{P}_k$. This folder contains interactive MATLAB Live Scripts (`.mlx`) demonstrating each method using a **Two-Tank** system as an intuitive working benchmark.

---

## Decoupled Pipeline Overview

In CHRONOS, the scheduling trajectory estimation is completely decoupled from the QP solver:
1. **Estimate $\hat{P}_k$:** Call one of the four `compute_schedul_*` functions.
2. **Build 3D Horizon Matrices:** Pass $\hat{P}_k$ to `project_matrix` using your affine matrix function handles (`compute_A`, `compute_B`, `compute_Bd`).
3. **Update & Solve:** Call `update_mpc_dynamics` followed by `mpc_solve`.

---

## The 4 Available Methods

### 1) The Frozen Approach (`01_Frozen_Approach.mlx`)
**[Read Detailed Documentation](./01_frozen_approach/readme.md)**

The simplest approach evaluates the exact nonlinear scheduling function $f_\rho(x)$ only once at the current measured state $x(k)$ and freezes (repeats) this vector across the entire prediction horizon $N_p$. Despite its lower prediction accuracy during fast nonlinear transients, it is widely used in practical applications due to its minimal computational cost and implementation simplicity.

```matlab
% 1. Estimate frozen scheduling trajectory (1 evaluation of f_rho)
Pk = compute_schedul_frozen(mpc, x_k, my_sched_fun);

% 2. Assemble 3D matrices across the horizon and solve QP once
A_lpv = project_matrix(compute_A, Pk, n_rho, Np);
mpc   = update_mpc_dynamics(mpc, A_lpv, B_lpv, []);
[u_k, mpc, iter] = mpc_solve(mpc, x_k, u_prev, y_ref, x_ter_ref, [], [], []);
```

---

### 2) The LPV-RTI / Warm-Start Approach (`02_LPV_RTI_Approach.mlx`)
**[Read Detailed Documentation](./02_LPV_RTI_approach/readme.md)**

Inspired by the Real-Time Iteration (RTI) paradigm, this approach achieves a highly accurate approximation with low computational overhead by exploiting the predicted state trajectory stored from the previous MPC solution (the warm-start vector). Instead of freezing the current state, each future state prediction is passed through $f_\rho(x)$ to populate $\hat{P}_k$.

> **Remark:** Because $f_\rho(x)$ is evaluated $N_p$ times per sampling instant, the computational cost of this method scales with the algebraic complexity of the nonlinear scheduling function $f_\rho(x)$.

```matlab
% 1. Map previous predicted state trajectory (warm-start) through f_rho
Pk = compute_schedul_lpv_rti(mpc, x_k, my_sched_fun, n_rho);

% 2. Assemble 3D matrices across the horizon and solve QP once
A_lpv = project_matrix(compute_A, Pk, n_rho, Np);
mpc   = update_mpc_dynamics(mpc, A_lpv, B_lpv, []);
[u_k, mpc, iter] = mpc_solve(mpc, x_k, u_prev, y_ref, x_ter_ref, [], [], []);
```

---

### 3) The LPV-SQP / Iterative Approach (`03_LPV_SQP_Approach.mlx`)
**[Read Detailed Documentation](./03_LPV_SQP_approach/readme.md)**

Also known as the iterative qLPV method, this approach implements a successive fixed-point refinement loop within the same sampling interval. It initializes a base trajectory, solves the MPC, extracts the newly predicted state trajectory to re-evaluate $\hat{P}_k^l$, and resolves the QP repeatedly until convergence ($\Vert \hat{P}_k^{l} - \hat{P}_k^{l-1} \Vert_\infty < \epsilon_{\text{cis}}$) or a maximum iteration budget $l_{\max}$ is reached.

```matlab
% 1. Iteratively refine Pk and re-solve the QP up to n_iter times
Pk = compute_schedul_lpv_sqp(mpc, x_k, u_prev, y_ref, x_ter_ref, [], ...
                               my_sched_fun, compute_A, compute_B, [], n_rho, n_iter);

% 2. Update dynamics with the converged trajectory and compute final control action
A_lpv = project_matrix(compute_A, Pk, n_rho, Np);
mpc   = update_mpc_dynamics(mpc, A_lpv, B_lpv, []);
[u_k, mpc, iter] = mpc_solve(mpc, x_k, u_prev, y_ref, x_ter_ref, [], [], []);
```

---

### 4) The Taylor Recursive Approach (`04_Taylor_Approach.mlx`)
**[Read Detailed Documentation](./04_Taylor_recursive_approach/readme.md)**

Based on a first-order Taylor series expansion of the scheduling function, this method extracts the base $\rho$ values from the time-shifted trajectory of the previous sampling instant and corrects them using the analytical Jacobian $\sigma_k = \frac{\partial f_\rho}{\partial x}$ multiplied by the state deviation $\Delta x$. To prevent numerical divergence under sharp nonlinearities, every extrapolated parameter is strictly clipped between physical bounds $[\underline{\rho}, \overline{\rho}]$.

```matlab
% 1. Extrapolate Pk via 1st-order Taylor expansion with saturation bounds
Pk = compute_schedul_taylor_recursive(mpc, x_k, n_rho, my_sched_fun, my_jacob_fun, rho_min, rho_max);

% 2. Assemble 3D matrices across the horizon and solve QP once
A_lpv = project_matrix(compute_A, Pk, n_rho, Np);
mpc   = update_mpc_dynamics(mpc, A_lpv, B_lpv, []);
[u_k, mpc, iter] = mpc_solve(mpc, x_k, u_prev, y_ref, x_ter_ref, [], [], []);
```

---

## Comparison Summary & Which Method Should I Use?

| Method | Function Name | $f_\rho(x)$ Evals / Sample | QP Solves / Sample | Prediction Accuracy | Computational Burden |
| :--- | :--- | :---: | :---: | :---: | :---: |
| **1. Frozen** | `compute_schedul_frozen` | $1$ | $1$ | Low (in fast transients) | Minimal |
| **2. LPV-RTI** | `compute_schedul_iterative_fast` | $N_p$ | $1$ | High | Low |
| **3. LPV-SQP** | `compute_schedul_iterative` | $1 + l_{\max}(N_p - 1)$ | $l_{\max} + 1$ | Very High | High |
| **4. Taylor** | `compute_schedul_recursive` | $1$ (+ $1$ Jacobian) | $1$ | Moderate to High | Very Low |

### Practical Recommendations
* **Start with Frozen or LPV-RTI:** They offer an exceptional trade-off, providing high stability and accuracy for the computational price of a single QP solve per sample, with zero tuning overhead.
* **Upgrade to LPV-SQP (Iterative):** Recommended when your plant is strongly nonlinear, operates near stability boundaries (e.g., pendulum swing-up), or experiences large reference steps/disturbances that invalidate the previous warm-start trajectory.
* **Choose Taylor Recursive:** Ideal when $f_\rho(x)$ is smooth ($\mathcal{C}^2$ with bounded derivatives) and algebraically expensive to evaluate $N_p$ times, allowing you to replace $N_p - 1$ nonlinear function evaluations with fast matrix-vector products.
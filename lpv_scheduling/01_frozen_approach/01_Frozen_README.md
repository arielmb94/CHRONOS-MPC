# 1. The Frozen (Gain-Scheduled) Approach

The **Frozen Approach** (often referred to as standard gain-scheduled MPC) is the most direct strategy to handle quasi-Linear Parameter-Varying (qLPV) models in **CHRONOS**. It converts the nonlinear prediction problem into a Linear Time-Invariant (LTI) model at each sampling instant by holding the current scheduling parameter constant across the entire prediction horizon.

---

## Conceptual Overview & Visual Scheme

At each sampling instant $k$, the controller measures the current state $x(k)$ and computes the exact scheduling parameter vector $\rho(k) = f_\rho(x(k))$. Because future states $x(k+1 \mid k), \dots, x(k+N_p-1 \mid k)$ are unknown prior to solving the optimization problem, the Frozen method simply assumes that the operating point remains unchanged over the next $N_p$ steps:

$$\hat{\rho}(k+j \mid k) = \rho(k), \quad \forall j \in \{0, 1, \dots, N_p-1\}.$$

<!-- [FIGURE PLACEHOLDER] -->
<p align="center">
  <img src="./figures/01_frozen_diagram.png" alt="Frozen Scheduling Trajectory Diagram" width="650"/>
  <br>
  <em>Figure 1: Conceptual illustration of the Frozen approach keeping the scheduling parameter constant over the prediction horizon.</em>
</p>

---

## Mathematical Formulation

Consider a discrete-time nonlinear plant embedded into the state-space qLPV representation:

$$x(k+1) = A\big(\rho(k)\big)x(k) + B\big(\rho(k)\big)u(k),$$

where $\rho(k) = f_\rho(x(k)) \in \mathcal{P} \subset \mathbb{R}^{n_\rho}$ is the vector of bounded scheduling variables.

### 1. Scheduling Trajectory Construction
At time step $k$, given the measured state $x(k) \in \mathbb{R}^{n_x}$, the scheduling proxy is evaluated once:

$$\rho(k) = f_\rho\big(x(k)\big).$$

The estimated scheduling trajectory vector $\hat{P}_k \in \mathbb{R}^{n_\rho N_p}$ is then assembled by repeating $\rho(k)$ vertically $N_p$ times:

$$\hat{P}_k = \text{col}\Big(\rho(k),\, \rho(k),\, \dots,\, \rho(k)\Big) = \mathbf{1}_{N_p} \otimes f_\rho\big(x(k)\big).$$

### 2. Prediction Model & Associated Uncertainty
Substituting $\hat{P}_k$ into the MPC equality constraints yields a stage-invariant prediction model along the horizon $\forall j \in \{0, \dots, N_p-1\}$:

$$x(k+j+1 \mid k) = A\big(\rho(k)\big)x(k+j \mid k) + B\big(\rho(k)\big)u(k+j \mid k).$$

While the first prediction step ($j=0$) is exact because $\hat{\rho}(k \mid k) = \rho(k)$, any movement of the state trajectory away from $x(k)$ introduces a scheduling mismatch $\Delta \rho(k+j \mid k) = \rho(k+j) - \rho(k)$ for $j \ge 1$. Consequently, the accumulated state prediction error $W_k = X_k - \hat{X}_k$ grows along the horizon:

$$W_k = \Big[\mathcal{A}(P_k) - \mathcal{A}(\hat{P}_k)\Big]x(k) + \Big[\mathcal{B}(P_k) - \mathcal{B}(\hat{P}_k)\Big]U_k,$$

which can degrade performance during rapid nonlinear transients if the sampling period $T_s$ or horizon $N_p$ is large.

---

## Implementation in CHRONOS

### Function Signature
```matlab
Pk = compute_schedul_frozen(mpc, x_k, f_rho)
```
* **`mpc`**: Initialized CHRONOS structure (uses `mpc.N` for the horizon length $N_p$).
* **`x_k`**: Current measured state vector ($n_x \times 1$).
* **`f_rho`**: Function handle `@(x)` returning the scheduling vector $\rho \in \mathbb{R}^{n_\rho \times 1}$.
* **`Pk`**: Output column vector of dimension $(n_\rho \cdot N_p) \times 1$.

### Complete Closed-Loop Snippet
```matlab
% Define scheduling function handle and affine matrix generator
my_sched_fun = @(x) [1 / sqrt(max(x(1), 1e-4)); 
                     1 / sqrt(max(x(2), 1e-4))];

compute_A = @(rho) eye(2) + Ts * [-k_tank * rho(1),  0; 
                                   k_tank * rho(1), -k_tank * rho(2)];

% Inside the closed-loop control loop at instant k:
Pk    = compute_schedul_frozen(mpc, x_prev, my_sched_fun);
A_lpv = project_matrix(compute_A, Pk, n_rho, mpc.N);
mpc   = update_mpc_dynamics(mpc, A_lpv, [], []);

[u_opt, mpc, iter] = mpc_solve(mpc, x_prev, u_prev, y_ref, x_ref, [], [], []);
```

---

## Pros, Cons & When to Use

* **Advantages:**
  * Requires only **one** evaluation of $f_\rho(x)$ and **one** QP solve per sampling instant.
  * Completely independent of previous MPC solutions (immune to warm-start corruption after large unmeasured disturbances).
* **Limitations:**
  * Ignores future parameter variations over the horizon ($\hat{\rho}(k+j \mid k) = \rho(k)$), which may fail to stabilize aggressive maneuvers (such as the Furuta pendulum swing-up).
* **Best suited for:** Systems with slow-to-moderate dynamics relative to the sampling rate, mildly nonlinear plants, or ultra-fast embedded loops with strict CPU budgets.
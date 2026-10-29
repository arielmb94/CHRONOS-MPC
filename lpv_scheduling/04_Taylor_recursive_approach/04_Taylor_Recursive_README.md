# 4. The Taylor-Based Recursive Extrapolation Approach

The **Taylor Recursive Approach** provides an analytical mechanism to update the scheduling parameter trajectory over the prediction horizon while retaining a **single QP solution** and **a single nonlinear function evaluation** per sampling instant. It corrects the scheduling trajectory from the previous time step using a first-order Taylor series expansion evaluated at the current measured state.

---

## Conceptual Overview & Visual Scheme

Assuming the scheduling proxy $f_\rho(x)$ is continuously differentiable ($\mathcal{C}^2$ with bounded derivatives), we can avoid re-evaluating the nonlinear function $f_\rho(\cdot)$ at every prediction step $j = 1, \dots, N_p-1$. Instead, at instant $k$:
1. The current scheduling parameter $\rho(k) = f_\rho(x(k))$ and the local analytical Jacobian matrix $\sigma_k = \left.\frac{\partial f_\rho}{\partial x}\right|_{x(k)}$ are evaluated **once** at the measured state $x(k)$.
2. For each future step $j$, the algorithm measures the spatial deviation $\Delta x(k+j)$ between the current state $x(k)$ and the previously predicted state $x^\star(k+j \mid k-1)$.
3. The previous scheduling parameter $\rho(k+j \mid k-1)$ is linearly corrected by $\sigma_k \Delta x(k+j)$ and projected onto the admissible physical set $\mathcal{P} = [\underline{\rho}, \overline{\rho}]$ via a saturation operator.

<!-- [FIGURE PLACEHOLDER] -->
<p align="center">
  <img src="./figures/04_taylor_recursive_diagram.png" alt="Taylor Recursive Extrapolation Diagram" width="650"/>
  <br>
  <em>Figure 1: First-order Taylor tangent correction and component-wise saturation between physical bounds.</em>
</p>

---

## Mathematical Formulation

Assume that $f_\rho : \mathbb{R}^{n_x} \mapsto \mathbb{R}^{n_\rho}$ is continuously differentiable with respect to the state vector $x$.

### 1. Local Jacobian Evaluation
At sampling instant $k$, given the measured state $x(k)$, the exact initial parameter is $\hat{\rho}(k \mid k) = f_\rho(x(k))$, and the analytical Jacobian matrix $\sigma_k \in \mathbb{R}^{n_\rho \times n_x}$ is evaluated at $x(k)$:

$$\sigma_k := \left.\frac{\partial f_\rho(x)}{\partial x}\right|_{x(k)} = \begin{bmatrix} \frac{\partial f_{\rho,1}}{\partial x_1} & \cdots & \frac{\partial f_{\rho,1}}{\partial x_{n_x}} \\ \vdots & \ddots & \vdots \\ \frac{\partial f_{\rho,n_\rho}}{\partial x_1} & \cdots & \frac{\partial f_{\rho,n_\rho}}{\partial x_{n_x}} \end{bmatrix}_{x = x(k)}.$$

### 2. Spatial State Deviation & First-Order Extrapolation
Let $x^\star(k+j \mid k-1)$ and $\rho(k+j \mid k-1)$ denote the predicted state and scheduling vectors saved from the previous sampling instant $k-1$. The spatial state deviation along the horizon is defined as:

$$\Delta x(k+j) := x(k) - x^\star(k+j \mid k-1), \quad \forall j \in \{1, \dots, N_p-1\}.$$

*(Note: For the terminal horizon slot where $k-1$ has no further shift, the last entry of the previous trajectory is repeated).*

Using a first-order Taylor expansion around the previous trajectory, the updated scheduling parameters for $j = 1, \dots, N_p-1$ are computed recursively as:

$$\hat{\rho}(k+j \mid k) = \text{sat}_{\mathcal{P}}\Big(\rho(k+j \mid k-1) + \sigma_k \Delta x(k+j)\Big),$$

where $\text{sat}_{\mathcal{P}}(\cdot)$ is the component-wise saturation (clipping) operator enforcing the physical bounds $\mathcal{P} = \big\{\rho \in \mathbb{R}^{n_\rho} : \underline{\rho} \le \rho \le \overline{\rho}\big\}$:

$$\text{sat}_{\mathcal{P}}(v)_i = \min\Big(\overline{\rho}_i,\, \max\big(\underline{\rho}_i,\, v_i\big)\Big), \quad \forall i \in \{1, \dots, n_\rho\}.$$

### 3. Why Saturation Clipping is Essential
> **Important Note on Smoothness & Clipping:** The Taylor expansion intrinsically assumes that $f_\rho(x)$ is smooth enough to be locally approximated by linear tangents ($\mathcal{C}^2$ boundedness). For highly exponential or rational nonlinearities (such as $1/\sqrt{x}$ near zero or $1/D(q)$ during a pendulum swing-up), $\sigma_k$ can become very large during fast transients, causing the linear tangent to overshoot wildly. The strict clipping operator $\text{sat}_{\mathcal{P}}(\cdot)$ prevents the extrapolated parameters from leaving the validity region of the qLPV embedding and protects the QP solver from numerical divergence.

---

## Implementation in CHRONOS

### Function Signature
```matlab
Pk = compute_schedul_recursive(mpc, x_k, n_rho, f_rho, jacob_rho, rho_min, rho_max)
```
* **`mpc`**: CHRONOS structure (stores previous scheduling trajectory and states).
* **`x_k`**: Current measured state vector ($n_x \times 1$).
* **`n_rho`**: Number of scheduling parameters ($n_\rho$).
* **`f_rho`**: Function handle `@(x)` returning $\rho \in \mathbb{R}^{n_\rho \times 1}$.
* **`jacob_rho`**: Function handle `@(x)` returning the Jacobian matrix $\sigma_k \in \mathbb{R}^{n_\rho \times n_x}$.
* **`rho_min`, `rho_max`**: Lower ($\underline{\rho}$) and upper ($\overline{\rho}$) physical bound vectors ($n_\rho \times 1$).

### Complete Closed-Loop Snippet
```matlab
% 1. Define scheduling function, analytical Jacobian, and physical bounds
my_sched_fun = @(x) [1 / sqrt(max(x(1), 1e-4)); 
                     1 / sqrt(max(x(2), 1e-4))];

my_jacob_fun = @(x) [-0.5 / (max(x(1), 1e-4)^1.5),  0; 
                      0,                           -0.5 / (max(x(2), 1e-4)^1.5)];

rho_min = my_sched_fun([1.00; 1.00]); % Minimum rho (at max water height)
rho_max = my_sched_fun([0.01; 0.01]); % Maximum rho (at min water height)

% 2. Inside the closed-loop control loop at instant k:
Pk    = compute_schedul_taylor_recursive(mpc, x_prev, n_rho, my_sched_fun, my_jacob_fun, rho_min, rho_max);
A_lpv = project_matrix(compute_A, Pk, n_rho, mpc.N);
mpc   = update_mpc_dynamics(mpc, A_lpv, [], []);

[u_opt, mpc, iter] = mpc_solve(mpc, x_prev, u_prev, xf, x_ref, [], [], []);
```

---

## Pros, Cons & When to Use

* **Advantages:**
  * Extremely fast online execution: requires only **1 evaluation of $f_\rho(x)$**, **1 evaluation of $\sigma_k$**, and $N_p - 1$ lightweight matrix-vector products, plus **1 QP solve**.
  * Avoids evaluating complex nonlinear functions $N_p$ times per sample.
* **Limitations:**
  * Requires the user to provide the analytical Jacobian matrix $\sigma_k(x)$ and explicit parameter bounds $[\underline{\rho}, \overline{\rho}]$.
  * Linear tangent extrapolation can produce spikes when $f_\rho(x)$ has sharp curvature or non-smooth regions.
* **Best suited for:** Systems with smooth, continuously differentiable scheduling functions ($\mathcal{C}^2$) and long prediction horizons $N_p$ where evaluating $f_\rho(x)$ repeatedly is computationally prohibitive.
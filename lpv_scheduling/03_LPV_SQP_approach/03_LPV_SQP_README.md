# 3. The LPV-SQP (Iterative Refinement) Approach

The **LPV-SQP Approach** (also known in the literature as the iterative function-based qLPV method) enforces strict consistency between the predicted state trajectory $\hat{X}_k$ and the scheduling parameter trajectory $\hat{P}_k$ **within the same sampling instant $k$**. Instead of relying solely on the previous sample's prediction, it performs a fixed-point refinement loop by repeatedly updating the LPV model and re-solving the QP.

---

## Conceptual Overview & Visual Scheme

At a single sampling instant $k$:
1. An initial trajectory guess $\hat{P}_k^0$ is generated using the Frozen approach at $x(k)$.
2. At sub-iteration $l$, the LPV model matrices are assembled from $\hat{P}_k^l$ and the QP is solved, yielding an updated optimal state trajectory $\hat{X}_k^l$.
3. Each predicted state in $\hat{X}_k^l$ is passed through $f_\rho(x)$ to build a refined scheduling trajectory $\hat{P}_k^{l+1}$.
4. Steps 2 and 3 are repeated until the change between consecutive scheduling trajectories falls below a convergence tolerance $\epsilon_{\text{cis}}$, or the maximum prescribed number of iterations $l_{\max}$ is reached.

<!-- [FIGURE PLACEHOLDER] -->
<p align="center">
  <img src="./figures/03_lpv_sqp_diagram.png" alt="LPV-SQP Iterative Refinement Diagram" width="650"/>
  <br>
  <em>Figure 1: Fixed-point iterative refinement loop of the LPV-SQP approach within a single sampling period.</em>
</p>

---

## Mathematical Formulation

Let $l \in \{0, 1, \dots, l_{\max}\}$ denote the inner sub-iteration index at sampling instant $k$.

### 1. Initialization ($l = 0$)
The algorithm initializes the scheduling trajectory using the exact measurement at $x(k)$ repeated across the horizon (Frozen guess):

$$\hat{P}_k^0 = \text{col}\Big(f_\rho\big(x(k)\big),\, f_\rho\big(x(k)\big),\, \dots,\, f_\rho\big(x(k)\big)\Big).$$

### 2. Successive Fixed-Point Refinement ($l = 1, \dots, l_{\max}$)
At each sub-iteration $l$:
1. **Assemble LPV Matrices:** Evaluate the 3D affine prediction matrices $\mathcal{A}(\hat{P}_k^{l-1})$ and $\mathcal{B}(\hat{P}_k^{l-1})$ along the horizon.
2. **Solve QP:** Solve the convex MPC optimization problem:
   
   $$\big(U_k^{\star, l},\, \hat{X}_k^l\big) = \arg\min_{U_k} J_k\big(x(k), U_k, \hat{P}_k^{l-1}\big),$$
   
   where $\hat{X}_k^l = \text{col}\big(x^l(k+1 \mid k),\, \dots,\, x^l(k+N_p \mid k)\big)$ is directly extracted from the CHRONOS primal solution vector (`mpc.x`) without needing an external forward simulator.
3. **Update Scheduling Trajectory:** Map the newly predicted states through the scheduling function:
   
   $$\hat{P}_k^l = \begin{bmatrix} f_\rho\big(x(k)\big) \\ f_\rho\big(x^l(k+1 \mid k)\big) \\ \vdots \\ f_\rho\big(x^l(k+N_p-1 \mid k)\big) \end{bmatrix}.$$

### 3. Stopping Criteria
The inner loop terminates as soon as the infinity norm (or Euclidean norm) of the difference between consecutive scheduling estimates satisfies the convergence threshold $\epsilon_{\text{cis}}$:

$$\left\Vert \hat{P}_k^l - \hat{P}_k^{l-1} \right\Vert_\infty < \epsilon_{\text{cis}},$$

or when $l = l_{\max}$. In embedded real-time applications, bounding the loop by $l_{\max}$ (typically $3$ to $5$ iterations) guarantees a deterministic worst-case execution time.

---

## Implementation in CHRONOS

### Function Signature
```matlab
Pk = compute_schedul_iterative(mpc, x_k, u_prev, y_ref, x_ter_ref, d_k, ...
                               f_rho, compute_A, compute_B, compute_Bd, n_rho, n_iter)
```
* **`mpc`**: Initialized CHRONOS structure.
* **`x_k`, `u_prev`, `y_ref`, `x_ter_ref`, `d_k`**: Current state, previous input, output reference, terminal reference, and measured disturbances required to solve the inner QPs.
* **`f_rho`**: Function handle `@(x)` returning $\rho \in \mathbb{R}^{n_\rho \times 1}$.
* **`compute_A`, `compute_B`, `compute_Bd`**: Function handles `@(rho)` returning the stage matrices (pass `[]` if a matrix is constant and does not depend on $\rho$).
* **`n_rho`**: Number of scheduling parameters ($n_\rho$).
* **`n_iter`**: Maximum number of inner SQP iterations ($l_{\max}$).

### Complete Closed-Loop Snippet
```matlab
n_iter = 5; % Maximum number of inner refinement iterations

% Inside the closed-loop control loop at instant k:
Pk = compute_schedul_lpv_sqp(mpc, x_prev, u_prev, xf, x_ref, [], ...
                               my_sched_fun, compute_A, [], [], n_rho, n_iter);

% Finalize dynamics with the converged trajectory and compute applied control
A_lpv = project_matrix(compute_A, Pk, n_rho, mpc.N);
mpc   = update_mpc_dynamics(mpc, A_lpv, [], []);

[u_opt, mpc, iter] = mpc_solve(mpc, x_prev, u_prev, xf, x_ref, [], [], []);
```

---

## Pros, Cons & When to Use

* **Advantages:**
  * Delivers the **highest prediction accuracy** among all four methods.
  * Self-corrects the entire horizon trajectory within the **current** sampling instant, even after abrupt reference steps or large external disturbances.
* **Limitations:**
  * Computationally expensive: requires up to $l_{\max} + 1$ full QP solutions per sampling interval.
* **Best suited for:** Highly nonlinear systems operating under aggressive transients, systems with longer sampling times $T_s$ that can afford multiple QP solves, or offline/high-accuracy benchmarking.
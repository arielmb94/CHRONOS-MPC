# LPV-MPC Scheduling Trajectory Estimation in CHRONOS

In order to solve the sampled LPV MPC optimization in CHRONOS, an accurate estimate for the future scheduling trajectory, $\hat{P}_k := \text{col}\Big(\rho(k),\rho(k+1|k),\dots,\rho(k+N_p-1|k)\Big)$, is required at each sampling instant. 

Depending on your system's nonlinearity and the available computational budget, CHRONOS provides four main approaches to compute this trajectory. This folder contains interactive MATLAB Live Scripts (`.mlx`) demonstrating each method using a **Two Tank** system as an intuitive working example.

### The 4 Available Methods

**1) The Frozen (Gain-Scheduling) Approach:** `01_Frozen_Approach.mlx`
The simplest approach computes the scheduling variables by evaluating the exact nonlinear scheduling function $f_\rho(x)$ at the current measured state $x(k)$, and simply freezing (repeating) this vector across the entire prediction horizon. Despite its inherent lower prediction accuracy during fast transients, this method is widely used in practical applications due to its minimal computational cost and implementation simplicity.

**2) The Warm-Start Evaluation Approach (Iterative Fast):** `02_Warm_Start_Approach.mlx`
A highly accurate approximation with lower computational overhead can be achieved by exploiting the predicted state trajectory from the previous MPC solution (the warm-start vector). Instead of freezing the current state, each future state prediction is passed through the nonlinear scheduling function $f_\rho(x)$ to populate $\hat{P}_k$.
> **Remark:** This method is subjected to the complexity of the nonlinear scheduling function $f_{\rho}(x)$. Consequently, the computational cost of this algorithm might increase for costly nonlinear functions.

**3) The Iterative (SQP) Approach:** `03_SQP_Approach.mlx`
Also known as the iterative qLPV method, this approach implements a successive refinement loop and relies on iteratively solving the MPC multiple times per sample. It computes a base trajectory, solves the MPC, uses the optimal solution to compute a new state prediction, and evaluates a new scheduling trajectory, repeating until convergence ($\Vert \hat{P}_k^{l} - \hat{P}_k^{l-1} \Vert < \epsilon_{\text{cis}}$).
This method yields extremely accurate predictions but is computationally heavy.

**4) The Recursive Extrapolation Approach:** `04_Recursive_Approach.mlx`
Based on a first-order Taylor expansion of the scheduling function, it extracts the base $\rho$ values from the shifted optimal trajectory of the previous time step and corrects them using the analytical Jacobian $\sigma_k$ multiplied by the spatial deviation $\Delta x$. To prevent numerical divergence for highly exponential nonlinearities, the algorithm enforces strict clipping between physical limits.

---
### Which method should I use?
In practical terms, we encourage you to start your implementation using the **Constant (Frozen) Approach** or the **Warm-Start Evaluation** method. They offer an exceptional balance, providing sufficiently high prediction accuracy with the computational cost of a single standard QP solve, despite being easy to set up. 

If your system is highly nonlinear, operates under strong disturbances, or experiences massive reference steps that invalidate the warm-start, upgrade to the **Iterative (SQP)** method to ensure robust performance. The **Recursive** method is also a great option when your system is continuously differentiable with bounded derivatives, providing high accuracy with low computational overhead.
## Scheduling Trajectory Estimates

In order to solve the sampled LPV MPC optimisation, an accurate estimate for the future scheduling trajectory, $\hat{P}_k := \text{col}\Big(\rho(k),\rho(k+1|k),\dots,\rho(k+N_p-1|k)\Big)$, is required at each sampling instant. Depending on the system's nonlinearity and the available computational budget, four main approaches can be implemented:

**1) The frozen (gain-scheduling) approach:**
The simplest approach computes the scheduling variables by evaluating the exact nonlinear scheduling function $f_\rho(x)$ at the current measured state $x(k)$, and simply freezing (repeating) this vector across the entire prediction horizon. Accordingly, the future scheduling trajectory estimate is taken as $\hat{P}_k = \text{col}(\rho(k), \cdots, \rho(k))$. Despite its inherent lower prediction accuracy during fast transients, this method is widely used in practical applications due to its minimal computational cost and implementation simplicity.

**2) The warm-start evaluation approach (Iterative Fast):**
A highly accurate approximation with lower computational overhead can be achieved by exploiting the predicted state trajectory from the previous MPC solution (the warm-start vector). Instead of freezing the current state, each future state prediction is passed through the nonlinear scheduling function $f_\rho(x)$ to populate $\hat{P}_k$. The vector is computed as:

$$ \hat{P}_k = \text{col}\Big(f_\rho(x(k)), f_\rho(x_0(k+1|k-1)), \dots, f_\rho(x_0(k+N_p-1|k-1))\Big) $$

The only strict requirement for this method is that the warm-start vector must be properly time-shifted at the end of each control step.

> **Remark:** This method is subjected to the complexity of the nonlinear scheduling function $f_{\rho}(x)$. Consequently, the computational cost of this algorithm might increase for costly nonlinear functions.

**3) The iterative (SQP) approach:**
Also known as the iterative qLPV method, this approach implements a successive refinement loop and relies on iteratively solving the MPC multiple times per sample. The procedure is as follows:
* First, $\hat{P}_k^0$ is computed using the frozen gain-scheduling approach, and the MPC is solved.
* Next, the optimal solution $U_k^\star$ is used to compute the state prediction $X_k^0$, and a new scheduling trajectory is evaluated using $\hat{P}_k^1 = f_\rho(X_k^0)$.
* This process is repeated until convergence is established, i.e., $\Vert \hat{P}_k^{l}-\hat{P}_k^{l-1}\Vert < \epsilon_{\text{cis}}$.

This method yields extremely accurate predictions but is computationally heavy since it requires solving the optimization problem multiple times per sampling instant. This computational burden may be an issue for embedded real-time systems with strict sampling times.

**4) The recursive extrapolation approach:**
This approach is based on a first-order Taylor expansion of the scheduling function. It extracts the base $\rho$ values from the shifted optimal trajectory of the previous time step ($P_{k-1}^*$) and corrects them using the analytical Jacobian $\sigma_k$ multiplied by the spatial deviation $\Delta x$. Specifically, the approach admits the following features:
* The method relies on the assumption that $f_\rho(\cdot)$ is continuously differentiable, with bounded derivatives.
* At each instant, the local first-order derivative is computed as $\sigma_k := \left.\frac{\partial f_\rho(x)}{\partial x}\right\vert_{\scriptscriptstyle x(k)}$.
* Incremental trajectories are computed based on the last sampled trajectories, and available data at the sampling instant $k$ is used to correct biased predictions from the previous instant.
* The extrapolation follows the recursive law: 

$$ \rho(k+j\mid k) \,=\, \rho_{\text{base}}(k+j) + \sigma_k \Delta x $$

* To prevent numerical divergence—such as tangents deviating drastically for highly exponential nonlinearities—the algorithm enforces strict clipping, bounding every extrapolated $\rho$ between its physical limits.

---
> **Final Remark:** In terms of numerical complexity, the frozen and warm-start evaluation approaches yield an MPC scheme with the complexity of a single QP, requiring only basic function evaluations. The recursive method computes a linear vector operation before solving a single QP. Conversely, the iterative method introduces the associated complexity of an SQP, requiring multiple QPs to be solved per sampling instant.
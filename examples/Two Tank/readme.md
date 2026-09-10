# Two Tank System Example

### Folder structure

In this folder you will find 3 files with the following objectives:

* *two_tank_init.m*: script to define the MPC problem using the CHRONOS init functions.
* *sim_two_tank_lpv.m*: script to simulate the Two Tank system in closed-loop using the CHRONOS mpc solver, at each iteration we use the CHRONOS update functions to adapt its internal Linear Parameter Varying model to the instantaneous water height level.
* *sim_two_tank_lti.m*: this script uses the same controller definition, nonlinear plant, reference, and solver iteration budget as *sim_two_tank_lpv.m*, but skips the online prediction-model update. The supplied simulation therefore isolates the effect of keeping the prediction model frozen at its initial operating point. In this scenario, that model mismatch produces a noticeable tracking offset.

### Example introduction

The two tanks system is a classical example used in control lectures. It consists of two inter-connected cylindrical tanks, where the first tank receives a steady flow of water, which can be regulated, and the second tank has a downstream discharge of water at its base. The dynamics of the two tank system is the following:

$$ \dot h_1 = u/A_b-\sqrt{2gh_1}/A_b $$
$$ \dot h_2 = \sqrt{2gh_1}/A_b-\sqrt{2gh_2}/A_b $$

where $h_1$ and $h_2$ are the water heights of each tank, $u$ is the controlled water massflow into tank 1 and $A_b$ is the tank area, equal for both tanks. The regulation objective is to control the water level on the second tank, e.g. our tracking target is:

$$ y = h_2$$

### From non-linear to linear time varying system description

The CHRONOS solver is designed to solve Nonlinear MPC by making use of Linear Parameter Varying (LPV) models. In order to transform the non-linear system into a LPV one, we do a linear embeddeding by dividing and multiplying the square root terms by the respective tank height:

$$ \sqrt{2gh_i} := \frac{\sqrt{2gh_i}}{h_i}h_i  $$

Substituting the linear embeddings on the non-linear dynamics equation, we arrive at the following state-space LPV description of the Two Tank system dynamics:

```math
\begin{bmatrix}
\dot h_1\\
\dot h_2
\end{bmatrix}
=
\begin{bmatrix}
-\frac{\sqrt{2gh_1}}{h_1A_b} & 0\\
\frac{\sqrt{2gh_1}}{h_1A_b} & -\frac{\sqrt{2gh_2}}{h_2A_b}
\end{bmatrix}
\begin{bmatrix}h_1\\h_2\end{bmatrix}
+\begin{bmatrix}1/A_b\\0\end{bmatrix}u
```

Note that if we expand the state-space LPV model, we recover the exact nonlinear dynamics of the two-tank system. This highlights a key advantage of the LPV representation: it captures the full nonlinear behavior of the system while casting it in a form compatible with convex optimization. As a result, CHRONOS can solve the nonlinear MPC problem exactly, using fast, reliable, and well-established convex optimization algorithms.

### MPC Definition

The controller uses a prediction horizon of $N=10$ samples. Its state, input,
and tracked output are

```math
s_k=\begin{bmatrix}h_{1,k}\\h_{2,k}\end{bmatrix},\qquad
u_k=\text{inlet flow},\qquad y_k=\begin{bmatrix}0&1\end{bmatrix}s_k=h_{2,k}.
```

At every sampling instant, CHRONOS solves the configured problem

```math
\begin{aligned}
\min_{s,u}\quad
J={}&(s_{ref,N}-s_N)^T P(s_{ref,N}-s_N) \\
&+\frac{1}{2}\sum_{k=1}^{N}(r_k-y_k)^TQ_e(r_k-y_k)
+\frac{1}{2}\sum_{k=0}^{N-1}\Delta u_k^T R_{du}\Delta u_k
\\
\text{subject to}\quad
&s_{k+1}=A_{d,k}s_k+B_du_k, && k=0,\ldots,N-1,\\
&0.01\le s_k\le 1, && k=1,\ldots,N,\\
&0\le u_k\le 10, && k=0,\ldots,N-1,\\
&-0.1\le\Delta u_k\le0.1, && k=0,\ldots,N-1.
\end{aligned}
```


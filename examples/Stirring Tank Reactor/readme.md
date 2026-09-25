# Continuous Stirred-Tank Reactor Example

### Folder structure

In this folder you will find the following files:

* *stirring_tank_init.m*: script to define the MPC problem using the CHRONOS init functions.
* *stirring_tank_sim_lpv.m*: script to simulate the reactor in closed loop while updating the LPV prediction model from the measured state.

### Example introduction

The Continuous Stirred-Tank Reactor (CSTR) model and parameters are taken from
[1]. Its nonlinear dynamics are

$$ \dot c =(1-c)/\theta_f - k c e^{-M/v} $$
$$ \dot v = (x_f-v)/\theta_f + kce^{-M/v}-\alpha u(v-x_c) $$

where $c$ is the product concentration, $v$ is the reactor temperature, $u$ is
the coolant flow rate, and $\theta_f$, $k$, $M$, $x_f$, $x_c$, and $\alpha$ are
reactor parameters. The controller regulates $c$ while tracking the
corresponding temperature $v$, so the tracking objective is the full state:

```math
s=\begin{bmatrix}c\\v\end{bmatrix}.
```

CHRONOS uses $y=s$ by default. Because this example tracks the full state,
`stirring_tank_init.m` does not need to define a separate output model.

### From nonlinear to LPV dynamics

The nonlinear terms can first be embedded in the following LPV model:

```math
\begin{bmatrix}
\dot c\\
\dot v
\end{bmatrix}
=
\begin{bmatrix}
-1/\theta_f-ke^{-M/v} & 0\\
ke^{-M/v} & -1/\theta_f
\end{bmatrix}
\begin{bmatrix}c\\v\end{bmatrix}
+
\begin{bmatrix}0\\-\alpha(v-x_c)\end{bmatrix}u
+
\begin{bmatrix}1/\theta_f\\x_f/\theta_f\end{bmatrix}
1
```

After discretization, CHRONOS supports models of the form

$$ x^{+}=Ax+Bu+B_d d $$

so constant or affine terms that do not fit in $A$ or $B$ can be represented
through $B_d d$. For this first embedding, the known input is simply $d=1$.

This direct embedding reproduces the nonlinear equations, but the resulting LPV
model is not controllable. The coolant input affects the temperature equation,
while the first row of $A$ has no term multiplying $v$. Once the scheduling
variables are fixed, the model therefore has no path from $u$ to $v$ and then
from $v$ to the concentration $c$.

The nonlinear model does contain this connection: the temperature appears in
the reaction term $-kce^{-M/v}$ of the concentration equation. To make that
dependence explicit, we approximate only this term with a first-order Taylor
expansion around the current operating point $(c^o,v^o)$:

$$
q^o=\frac{k c^o M e^{-M/v^o}}{(v^o)^2},\qquad
-kce^{-M/v}\approx-k e^{-M/v^o}c-q^o(v-v^o).
$$

The new term $-q^o v$ introduces the missing temperature-to-concentration
coupling, while $q^o v^o$ is the offset required for the approximation to match
the nonlinear term at the expansion point. Substituting this approximation gives
the new LPV model

```math
\begin{bmatrix}\dot c\\ \dot v\end{bmatrix}=
\begin{bmatrix}
-1/\theta_f-ke^{-M/v^o} & -q^o\\
ke^{-M/v^o} & -1/\theta_f
\end{bmatrix}
\begin{bmatrix}c\\v\end{bmatrix}
+\begin{bmatrix}0\\-\alpha(v^o-x_c)\end{bmatrix}u
+\begin{bmatrix}1/\theta_f&q^o\\x_f/\theta_f&0\end{bmatrix}
\begin{bmatrix}1\\v^o\end{bmatrix}.
```

At each control sample, $(c^o,v^o)$ is set to the measured reactor state. The
physical states $c$ and $v$ are kept unchanged, and the linearization point is
passed as the known input $d=[1\;v^o]^T$. This illustrates the flexibility of
the MPC-LPV approach: most of the dynamics retain their LPV embedding, while
only the term that causes the controllability problem is locally linearized.
CHRONOS supports this hybrid model by allowing $A$, $B$, $B_d$, and $d$ to be
updated online.

### MPC Definition

The controller tracks the full reactor state using the default output $y=s$:

```math
s_k=\begin{bmatrix}c_k\\v_k\end{bmatrix},\qquad
r_k=\begin{bmatrix}c_{ref,k}\\v_{ref,k}\end{bmatrix}
```

over a horizon of $N=15$ samples:

```math
\begin{aligned}
\min_{s,u}\quad
&\frac{1}{2}\sum_{k=1}^{N}(r_k-s_k)^TQ_e(r_k-s_k)\\
\text{subject to}\quad
&s_{k+1}=A_{d,k}s_k+B_k u_k+B_{d,k}d_k,
&&k=0,\ldots,N-1,\\
&0\le s_k\le1, &&k=1,\ldots,N,\\
&0\le u_k\le1, &&k=0,\ldots,N-1,\\
&-0.1\le\Delta u_k\le0.1, &&k=0,\ldots,N-1.
\end{aligned}
```

### References

[1] Nonhoff, M., Köhler, J., & Müller, M. A. (2024). [Online convex optimization for constrained control of nonlinear systems](https://arxiv.org/abs/2412.00922). arXiv:2412.00922.

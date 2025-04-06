# Stirring Tank System Example

### Folder structure

In this folder you will find the following files:

* *stirring_tank_init.m*: script to define the MPC problem using the CHRONOS init functions.
* *stirring_tank_sim_lpv.m*: script to simulate the Stirring Tank system in closed-loop using the CHRONOS mpc solver, at each iteration we use the CHRONOS update functions to adapt its internal Linear Time Varying model to the instantaneous Stirring Tank states.

### Example introduction

We borrowed the Continous Stirring Tank Reactor (CSTR) example and parameters from [1].  The nonlinear dynamics of the CSTR are given by:

$$ \dot c =(1-c)/\theta_f - k c e^{-M/v} $$
$$ \dot v = (x_f-v)/\theta_f + kce^{-M/v}-\alpha u(v-x_c) $$

where $c$ is the product concentration, $v$ is the CSTR temperature, $u$ is the coolant flow rate and $\theta_f$, $k$, $M$, $x_f$, $x_c$ and $\alpha_f$ are reactor parameters. The control target is to regulate the product concentration $c$, however, given that the product concentration evolution is regulated through the CSTR temperature $v$, both $c$ and $v$ must be controlled. Thus, our tracking objective is the full state vector:

$$ y = \left [\begin{array}{c} c\\\ v \end{array} \right ] $$

### From non-linear to linear time varying system description

Given that the states $c$, $v$ and the control input $u$ appear linearly on the model equation, we can rewrite the CSTR dynamics as the following Linear Time Varying (LTV) model:

$$ \left [\begin{array}{c} \dot c\\\ \dot v \end{array} \right ] =
\left [\begin{array}{cc}  -1/\theta_f-ke^{-M/v} & 0\\\ 
ke^{-M/v} &-1/\theta_f\end{array} \right ]
\left [\begin{array}{c}c\\\  v \end{array} \right ] + 
\left [\begin{array}{c} 0\\\ -\alpha(v-x_c) \end{array} \right ]u + 
\left [\begin{array}{c} 1/\theta_f\\\ x_f/\theta_f \end{array} \right ]
1$$

Note that we used the capability of the CHRONOS solver to work with systems of the form:

$$ x^{+}=Ax+Bu+B_dd $$

to place on the disturbance matrix $B_d$ the terms which cannot be directly grouped on the $A$ and $B$ matrices of the LTV model, using a fake known disturbance input set to 1.

 There is a problem however with the previous LTV model, it is not controllable. This highligh one of the challenges with translating a non-linear model into a LTV one, if not done carefully we might built LTV models with negetative characteristics, e.g. lack of controllability or observability, which are not true of the nonlinear system.

 To find a controllable LTV representation lets replace the term $-kce^{-M/v}$ from the $\dot c$ differential equation by its first order taylor expansion with repect the CSTR temperature $v$:

$$ -kce^{-M/v} \approx  -kce^{-M/{v^o}} - \frac{kce^{-M/{v^o}}}{{v^o}^2}(v-{v^o}) $$
 
In practice, we will update the Taylor expansion point such that we pick $v^o = v(k)$. Replacing the Taylor expansion on the nonlinear model equation for $\dot c$ we then obtain the following LTV model:

$$ \left [\begin{array}{c} \dot c\\\ \dot v \end{array} \right ] =
\left [\begin{array}{cc}  -1/\theta_f-ke^{-M/v^o} & - \frac{kce^{-M/v^o}}{{v^o}^2}\\\ 
ke^{-M/v} &-1/\theta_f\end{array} \right ]
\left [\begin{array}{c}c\\\  v \end{array} \right ] + 
\left [\begin{array}{c} 0\\\ -\alpha(v-x_c) \end{array} \right ]u + 
\left [\begin{array}{cc} 1/\theta_f & \frac{kce^{-M/v^o}}{{v^o}^2}\\\ x_f/\theta_f & 0 \end{array} \right ]
\left [\begin{array}{c}1\\\  v^o \end{array} \right ]$$
 
Note that we made use again of the input disturbance matrix $B_d$, this time to handle the term $`\frac{kce^{-M/{v^o}}}{{v^o}^2}{v^o}`$ from the Taylor expansion. Now we have a controllable LTV model that captures almost exactly the nonlinear dynamics of the CSTR system. We can then use the CHRONOS solver to define and solve the MPC problem for controlling the CSTR system using fast and realiable convex optimization algorithms.

### MPC Definition

In order to control the height of the second tank we solve at each iteration the following MPC problem using the CHRONOS solver:

$$\min_{u,x}J = (r-x_N)^TP(r-x_N) + \sum_{i=1}^{N-1} (r-x_i)^TQ_{e}(r-x_i) + \sum_{i=0}^{N_{ctr}-1}\Delta u_0^TdR_u\Delta u_0$$

s.t.

$$ x^+=A(c,v)x+B(v)u+B_d(v)d$$
$$ c,v \in [0,\\\ 1]$$
$$ u \in [0,\\\ 1] $$

### References

[1] Nonhoff, M., Köhler, J., & Müller, M. A. (2024). Online convex optimization for constrained control of nonlinear systems. arXiv preprint arXiv:2412.00922.


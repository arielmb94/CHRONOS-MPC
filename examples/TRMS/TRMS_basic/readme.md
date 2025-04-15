# Twin Rotor MIMO System (TRMS) Basic MPC Example

### Folder structure

In this folder you will find the following files:

* *TRMS_init.m*: script to define the MPC problem using the CHRONOS init functions.
* *TRMS_sim_lpv.m*.m: script to simulate the TRMS in closed-loop using the CHRONOS MPC solver, at each iteration we use the CHRONOS update functions to adapt its internal Linear Parameter Varying model to the instantaneous TRMS states.
* *qLPV_TRMS_SS.m*: computes the LPV model of the TRMS based on the current values of the state vector. The LPV model is extracted from the nonlinear model provided in [1].
* *compute_ref.m*: computes the state references based on the angle set points for the TRMS and from the equilibrium equations as given in [1].

### Example introduction

The TRMS is a simplified representation of a helicotper, it has no translation, however, it can rotate freely on the horizontal and vertical frames. It counts with a main rotor to control the vertical angle enabling the TRMS to pitch and a tail rotor to control the horizontal angle, which allows changes on the TRMS yaw angle. Each rotor operated by a dedicated DC motor. 

Existing multiple representations of the TRMS dynamics available in the literature, we borrowed the nonlinear model and parameters identified in [1]. The nonlinear model presented in [1] captures very well the nonlinear dynamics of the vertical and horzizontal TRMS dynamics, their couplings and the effect of friction forces, which are represented by discontinuous equations taking into account the differences between negative and positive displacement directions. The quality of the identification work carried by the authors in [1] is demonstrated by the perfect matching between simulated model response and real data from the physical TRMS behaviour. 

Without detailing the full nonlinear equations of the model terms, the LPV model extracted from the nonlinear model given in [1] can be represented as:

$$ \dot x=
\left [\begin{array}{cccccc}  
a_{11}(\rho) & 0 & 0 & 0 & 0 & 0 \\\ 
a_{21}(\rho) & a_{22}(\rho) & a_{23}(\rho) & a_{24}(\rho) & a_{25}(\rho) & a_{26}(\rho) \\\
0 & a_{32} & 0 & 0 & 0 & 0 \\\
0 & 0 & 0 & a_{44}(\rho) & 0 & 0 \\\
0 & a_{52}(\rho) & 0 & a_{54}(\rho) & a_{55} & a_{56}(\rho) \\\
0 & 0 & 0 & 0 & a_{65} & 0
\end{array} \right ]
x + 
\left [\begin{array}{cc} 
b_{11} & 0 \\\
0 & b_{22}(\rho) \\\
0 & 0 \\\
0 & b_{42} \\\
0 & 0 \\\
0 & 0 
\end{array} \right ]u$$

The state vector is $x = [\omega_h,\Omega_h,\theta_h,\omega_v,\Omega_v,\theta_v]^T $, where

* $\omega_h$: angular speed of the tail rotor DC fan
* $\Omega_h$: TRMS angular speed on the horizontal frame
* $\theta_h$: TRMS horizontal angle
* $\omega_v$: angular speed of the main rotor DC fan
* $\Omega_v$: TRMS angular speed on the vertical frame
* $\theta_v$: TRMS vertical angle

The input vector is  $u = [u_h,u_v]^T$, where

* $u_h$: DC voltage applied to the tail rotor fan
* $u_v$: DC voltage applied to the main rotor fan

Finally, note that parameter varying terms have been made explicit by showing on the state-space model their depedency on the varying parameter vector $\rho$. On the LPV model developed for the TRMS, the varying parameter vector $\rho$ is formed by $\rho = [\omega_h,\Omega_h,\theta_h,\omega_v,\theta_v]^T$. The file *qLPV_TRMS_SS.m* allows to compute the LPV model matrices for a given value of the varying parameters $\rho$.

### MPC Definition

The control objective is to control the TRMS horizontal and vertical angles. This is achieved by solving at each time step the following MPC problem using the CHRONOS solver:

$$\min_{u,x}J = (r-x_N)^TP(r-x_N) + \sum_{i=1}^{N-1} (r-x_i)^TQ_{e}(r-x_i) + \sum_{i=0}^{N_{ctr}-1}\Delta u_i^TdR_u\Delta u_i$$

s.t.

$$ x^+=A(\rho)x+B(\rho)u$$
$$  \left [\begin{array}{c} 
-2.9 \\\
-1.0 \\\
-1.7 \\\
-1.6 \\\
-0.6 \\\
-0.5
\end{array} \right ]
\leq x \leq
\left [\begin{array}{c} 
2.9 \\\
1.0 \\\
1.2 \\\
1.6 \\\
0.6 \\\
1.0 
\end{array} \right ]$$
$$  \left [\begin{array}{c} 
-2.5 \\\
-2.0 
\end{array} \right ]
\leq u \leq
\left [\begin{array}{c} 
2.5 \\\
2.0 
\end{array} \right ]$$

Due to the highly nonlinear dynamics of the TRMS and the couplings between its vertical and horizontal motions, it is required to provide tracking references for all states to obtain good control performance. The state references are computed as follows:
* $\theta_h^{ref}$, $\theta_v^{ref}$: set points provided to the MPC
* $\Omega_h^{ref}$, $\Omega_v^{ref}$: computed from the angle error value and a time constant $\tau$ as:

$$ \Omega_i = \frac{\theta_i^{ref}-\theta_i}{\tau} $$ 

* $\omega_h^{ref}$, $\omega_v^{ref}$: computed from the equilibrium equations, e.g. $\dot\omega_i=0$, of the nonlinear differential equations model. The equilibrium equations were obtained from [1].


### References

[1] Rotondo, D., Nejjari, F., & Puig, V. (2013). Quasi-LPV modeling, identification and control of a Twin Rotor MIMO System. Control Engineering Practice, 21(6), 829-846.


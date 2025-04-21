# Twin Rotor MIMO System (TRMS) Custom Reference MPC Example

### Folder structure

In this folder you will find the following files:

* *TRMS_refMPC_init.m*: script to define the MPC problem using the CHRONOS init functions.
* *TRMS_refMPC_sim_lpv.m*.m: script to simulate the TRMS in closed-loop using the CHRONOS mpc solver, at each iteration we use the CHRONOS update functions to adapt its internal Linear Parameter Varying model to the instantaneous TRMS states.
* *qLPV_TRMS_refMPC_SS.m*: computes the LPV model of the TRMS based on the current values of the state vector. The LPV model is extracted from the nonlinear model provided in [1].

### TRMS introduction

The TRMS is a simplified representation of a helicotper, it has no translation, however, it can rotate freely on the horizontal and vertical frames. It counts with a main rotor to control the vertical angle enabling the TRMS to pitch and a tail rotor to control the horizontal angle, which allows changes on the TRMS yaw direction. Each rotor operated by a dedicated DC motor. 

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

The state vector is $x = [\omega_h,\Omega_h,\theta_h,\omega_v,\Omega_v,\theta_v]^T $, with

* $\omega_h$: angular speed of the tail rotor DC fan
* $\Omega_h$: TRMS angular speed on the horizontal frame
* $\theta_h$: TRMS horizontal angle
* $\omega_v$: angular speed of the main rotor DC fan
* $\Omega_v$: TRMS angular speed on the vertical frame
* $\theta_v$: TRMS vertical angle

The input vector is  $u = [u_h,u_v]^T$, with

* $u_h$: DC voltage applied to the tail rotor fan
* $u_v$: DC voltage applied to the main rotor fan

Finally, note that parameter varying terms have been made explicit by showing on the state-space model their depedency on the varying parameter vector $\rho$. On the LPV model developed for the TRMS, the varying parameter vector $\rho$ is formed by $\rho = [\omega_h,\Omega_h,\theta_h,\omega_v,\theta_v]^T$.

### Example introduction

Different to the basic MPC example, on this example we want to stop computing the references for the rotor fan speeds based on equilibrium equations. Instead, we will make use of CHRONOS cost function custom performance vectors $z$ in order to make the MPC solver compute appropiate references for the rotor fan speed states internally.

Let's consider the extended control action vector $u = [u_h,u_v,\omega_h^{ref},\omega_v^{ref}]^T$ with:

* $\omega_h^{ref}$: the tail rotor fan speed setpoint to be computed by the MPC
* $\omega_v^{ref}$: the main rotor fan speed setpoint to be computed by the MPC

To accomodate the new MPC actions we will decouple the rotor fan speeds states from the TRMS horizontal and vertical dynamics. In their place, terms associated to the rotor fan speed will now be moved to be terms multupliying the newly introduced setpoint actions. Following this, let's redifine the LPV model as:

$$ \dot x=
\left [\begin{array}{cccccc}  
a_{11}(\rho) & 0 & 0 & 0 & 0 & 0 \\\ 
0 & a_{22}(\rho) & a_{23}(\rho) & 0 & a_{25}(\rho) & a_{26}(\rho) \\\
0 & a_{32} & 0 & 0 & 0 & 0 \\\
0 & 0 & 0 & a_{44}(\rho) & 0 & 0 \\\
0 & a_{52}(\rho) & 0 & 0 & a_{55} & a_{56}(\rho) \\\
0 & 0 & 0 & 0 & a_{65} & 0
\end{array} \right ]
x + 
\left [\begin{array}{cccc} 
b_{11} & 0 & 0 & 0 \\\
0 & b_{22}(\rho) & a_{21}(\rho) & a_{24}(\rho) \\\
0 & 0 & 0 & 0 \\\
0 & b_{42} & 0 & 0 \\\
0 & 0 & 0 & a_{54}(\rho) \\\
0 & 0  & 0 & 0
\end{array} \right ]u$$

Finally, we must specify to CHRONOS that it must minimize the difference between the setpoint actions $\omega_i^{ref}$ and the rotor fan speeds states $\omega_i$, e.g. minimize $\omega_i^{ref}-\omega_i$. For this, we will use CHRONOS custom performance vectors

$$ z = C_zx + D_zu+D_{zd}d $$

in order to create a new term on the MPC cost function. The performance vectors indicating the fan speed errors to be minimized by CHRONOS can be defined as:

$$z = \omega_i^{ref}-\omega_i=
\left [\begin{array}{cccccc} 
-1 & 0 & 0 & 0 & 0 & 0 \\\
0 & 0 & 0 & -1 & 0 & 0
\end{array} \right ]x+
\left [\begin{array}{cccc} 
0 & 0 & 1 & 0  \\\
0 & 0 & 0 & 1
\end{array} \right ]u
$$

Once the performance vector $z$ has been defined, CHRONOS will autimatically compute the gradient and Hessian values associated to the customly defined cost function and take them into account during the MPC optimization problem.

### MPC Definition

The control objective is to control the TRMS horizontal and vertical angles while computing setpoints for the rotors fan speeds states. This is achieved by solving at each time step the following MPC problem using the CHRONOS solver:

$$\min_{u,x}J = (r-x_N)^TP(r-x_N) + \sum_{i=1}^{N-1} (r-x_i)^TQ_{e}(r-x_i) + \sum_{i=0}^{N_{ctr}-1}\Delta u_i^TdR_u\Delta u_i + \sum_{i=1}^{N-1} z_i^TQ_{z}z_i$$

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

* $\omega_h^{ref}$, $\omega_v^{ref}$: computed internally by the CHRONOS MPC solver.


### References

[1] Rotondo, D., Nejjari, F., & Puig, V. (2013). Quasi-LPV modeling, identification and control of a Twin Rotor MIMO System. Control Engineering Practice, 21(6), 829-846.


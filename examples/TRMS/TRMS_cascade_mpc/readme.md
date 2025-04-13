# Twin Rotor MIMO System (TRMS) Custom Reference MPC Example

### Folder structure

In this folder you will find the following files:

* *TRMS_cascade_mpc_init.m*: script to define the MPC problem using the CHRONOS init functions.
* *TRMS_cascade_mpc_lpv.m*.m: script to simulate the TRMS in closed-loop using the CHRONOS mpc solver, at each iteration we use the CHRONOS update functions to adapt its internal Linear Parameter Varying model to the instantaneous TRMS states.
* *qLPV_TRMS_cascade_mpc_SS.m*: computes the LPV model of the TRMS based on the current values of the state vector. The LPV model is extracted from the nonlinear model provided in [1].

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

In this example, we explore the performance benefits of using a full cascade of MPC controllers. Although it might seem counterintuitive to apply MPC in the inner control loops, this architecture offers several key advantages:

* MPC controllers compute a whole sequence of control actions but only the first action is often used. By cascading MPC controllers, the inner layers can profit from the whole sequence of setpoints computed by the outer layer MPC.
* Inner loops often deal with actuator dynamics, which may be nonlinear and difficult to handle with traditional PID or linear controllers. Implementing nonlinear MPC strategies, such as LPV MPC as used in CHRONOS, can significantly enhance control performance in these scenarios.
* Decomposing a large MIMO control problem into smaller MPC subproblems can reduce computational complexity. Each subproblem is simpler and can be solved faster than a monolithic MPC formulation. Additionally, the prediction horizons of each layer can be tuned independently, offering further flexibility and potentially faster overall response compared to a single, centralized MPC controller.



For the TRMS MPC cascade architecture, we will create an outer layer MPC controlling the horizontal and vertical dynamics of the TRMS. The output of this MPC will be setpoints for the tail and main rotors fan speeds respectively. Then, the voltage to be applied to each rotor DC motor will be computed by a dedicated SISO MPC.

The model for the MIMO MPC will then have a reduced state vector:

* $\Omega_h$: TRMS angular speed on the horizontal frame
* $\theta_h$: TRMS horizontal angle
* $\Omega_v$: TRMS angular speed on the vertical frame
* $\theta_v$: TRMS vertical angle

and its control inputs will be:

* $\omega_h^{ref}$: the tail rotor fan speed setpoint to be computed by the MPC
* $\omega_v^{ref}$: the main rotor fan speed setpoint to be computed by the MPC

resulting on the following model:

$$ \dot x=
\left [\begin{array}{cccc}  
a_{22}(\rho) & a_{23}(\rho) & a_{25}(\rho) & a_{26}(\rho) \\\
a_{32} & 0 & 0 & 0 \\\
a_{52}(\rho) & 0 & a_{55} & a_{56}(\rho) \\\
0 & 0 & a_{65} & 0
\end{array} \right ]
x + 
\left [\begin{array}{cc} 
a_{21}(\rho) & a_{24}(\rho) \\\
0 & 0 \\\
0 & a_{54}(\rho) \\\
0 & 0
\end{array} \right ]\omega^{ref}+
\left [\begin{array}{c}
b_{22}(\rho) \\\
0 \\\
0 \\\
0
\end{array} \right ]u_v $$

Note that the voltage to the main rotor $u_v$ appears in the model as a measured disturbance, not as a control input. Making use of the fact that CHRONOS accepts models of the form

$$ x^+ = Ax+Bu+B_dd $$

allows us to account for known or measurable perturbations to our system, as the case of the coupling term between main rotor voltage and TRMS horizontal dynamics.

The tail rotor MPC has a single order LPV state-space model with state:

* $\omega_h$: angular speed of the tail rotor DC fan

and control action:

* $u_h$: DC voltage applied to the tail rotor fan

resulting on the LPV model:

$$ \dot x_h =
a_{11}(\rho) x_h + 
b_{11} u_h $$

Similarly, the model used for main rotor MPC is a single order model with state:

* $\omega_v$: angular speed of the main rotor DC fan

and control input:

* $u_v$: DC voltage applied to the main rotor fan

resultin on:

$$ \dot x_v =
a_{44}(\rho) x_v + 
b_{42}u_v $$

### MPC Definition

For the outer MIMO MPC, its objective is to compute appropiate setpoints for the rotor fans speed $\omega^{ref}$. This is achieved by solving at each time step the following MPC problem using the CHRONOS solver:

$$\min_{u,x}J = (r-x_N)^TP(r-x_N) + \sum_{i=1}^{N-1} (r-x_i)^TQ_{e}(r-x_i) + \sum_{i=0}^{N_{ctr}-1}\Delta {\omega_i^{ref}}^TdR_u\Delta \omega_i^{ref} $$

s.t.

$$ x^+=A(\rho)x+B(\rho)u$$
$$  \left [\begin{array}{c} 
-1.0 \\\
-1.7 \\\
-0.6 \\\
-0.5
\end{array} \right ]
\leq x \leq
\left [\begin{array}{c} 
1.0 \\\
1.2 \\\
0.6 \\\
1.0 
\end{array} \right ]$$
$$  \left [\begin{array}{c} 
-2.9 \\\
-1.6 
\end{array} \right ]
\leq \omega^{ref} \leq
\left [\begin{array}{c} 
2.9 \\\
1.6 
\end{array} \right ]$$

The tail rotor SISO MPC uses the full sequence of actions computed by the outer MIMO MPC to solve the following MPC at each iteration:

$$\min_{u,x}J = (\omega_{h_N}^{ref}-x_{h_N})^TP(\omega_{h_N}^{ref}-x_{h_N}) + \sum_{i=1}^{N-1} (\omega_{h_i}^{ref}-x_{h_i})^TQ_{e}(\omega_{h_i}^{ref}-x_{h_i}) + \sum_{i=0}^{N_{ctr}-1}\Delta u_{h_i}^TdR_u\Delta u_{h_i} $$

s.t.

$$ x_h \in [-2.9,2.9] $$
$$ u_h \in [-2.5, 2.5]$$

Very similarly, the SISO MPC for the main rotor is defined in CHRONOS as:

$$\min_{u,x}J = (\omega_{v_N}^{ref}-x_{v_N})^TP(\omega_{v_N}^{ref}-x_{v_N}) + \sum_{i=1}^{N-1} (\omega_{v_i}^{ref}-x_{v_i})^TQ_{e}(\omega_{v_i}^{ref}-x_{v_i}) + \sum_{i=0}^{N_{ctr}-1}\Delta u_{v_i}^TdR_u\Delta u_{v_i} $$


$$ x_v \in [-1.6,1.6] $$
$$ u_v \in [-2.0, 2.0] $$


### References

[1] Rotondo, D., Nejjari, F., & Puig, V. (2013). Quasi-LPV modeling, identification and control of a Twin Rotor MIMO System. Control Engineering Practice, 21(6), 829-846.


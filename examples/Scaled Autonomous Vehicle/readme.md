# Scaled Autonomous Vehicle Example

### Folder structure

In this folder you will find the following files:

* *main_SAV_Simulator.m*: main script, loads the simulator parameters, run the vehicle simulation and generate plots.
* *sav_lat_controller_init.m*: script defining the CHRONOS-MPC for the control of the vehicle lateral dynamics.
* *sav_lat_mpc_solve.m*: function used online call the CHRONOS-MPC solver for the vehicle's lateral dynamics control.
* *update_BM.m*: function to generate the LPV matrices of the vehicle frozen at the given velocity.
* *SAV Simulator*: folder containing all the files needed for the simulation of the Scaled Autonomous Vehicle (SAV) plarform. The main contained files being *CAR_SIM_MPC.slx* and *simNLPVCalpha.m*, the simulink file containing the SAV simulator and the function implementing the nonlinear dynamics of the vehicle.

### SAV Simulator Environment

The example provided here is concerned with the Autonomous Trajectory Tracking for the Scaled Autonomous Vechile (SAV) platform. The SAV platform is a 1:15 scale vehicle used in GIPSA-Lab, Grenoble (France), for research on autonomous vehicle control. For full details on the platform please check [1], [2]. The SAV car has two electric motors on each of the rear wheels for propulsion and uses a servo motor to steer the front wheels. The platform uses a motion capture system to obtain its position and velocities, acting as the main sensor.

To speed up the development of control strategies on the platform, the simulator environment provided in this example was created. The simulator captures to a high level of fidelity the key aspects of the SAV platform:

* Coupled Lateral / Longitudinal dynamics with nonlinear tyre models. All parameters obtained from identification process of the SAV car dynamics using extensive datasets.
* Actuator dynamics and actuator controllers mimic the real issues and implementation as done in the SAV platfrom. The steering servo motor is modeled as a slow second order system with 0.18 seconds input time-delay. The rear-wheel electric burshless motors + ESC are modeled jointly as a DC motor each, with the control system for each motor modeled as implemented exactly on the real platform: large input delay, high noise encoders with poor numerical resolution and random count skips, single current sensor for both motors, PI observer for accurate real-time estimation and Hinfinity controller for each motor.
* The software stack is accuretly replicated. The simulator uses the same methods as in the SAV platform software to obtain the velocities from the global position given by the motion capture system, the software differential used to distribute the speed between the rear wheels is present on the simulator and the reference generator for the lateral control based on the Pure-Pursuit method is implemented as in the SAV platform itself.

All in all, the key characteristics and details of the SAV platform are captured on the simulator shared on this example. Thanks to the accuracy of the simulator, the results obtained in simulation carried almost exactly to the real platform.

### Autonomous Trajectory Tracking: Lateral Control

For trajectory tracking we need to pass to the MPC controller an adequate reference. The reference generation for the lateral controller is carried as follows on the SAV platform:

1. The car position is captured by the motion capture system.
2. The car position is extrapolated at a look-ahead point based on the current speed $v$ and at that point the vehicle position is compared with the desired reference trajectory. The look-ahead distance computed as $L=t_p*v$, where $t_p$ is the look-ahead time in seconds.
3. Based on the lateral and angle errors at the look-ahead point, the Pure Pursuit algorithm is then used to generate a desired yaw rate reference (how fast the vehicle should turn) that the vehicle should follow in order to track the trajectory.

The lateral MPC controller task is then to ensure that the vehicle follows the desired yaw rate signal in order to generate an adequate steering angle command to the servo actuator of the SAV car. The model used on the lateral MPC controller is the well known single-track bicycle model:

$$ \left [\begin{array}{c} \dot v_y\\\ \ddot \psi \end{array} \right ] =
\left [\begin{array}{cc}  -\frac{C_f+C_r}{m  v_x} & -v_x-\frac{C_fl_f-C_rl_r}{mv_x}\\\ 
-\frac{C_fl_f-C_rl_r}{I_zv_x} & -\frac{C_fl_f^2+C_rl_r^2}{I_zv_x} \end{array} \right ]
\left [\begin{array}{c}v_y\\\ \dot \psi \end{array} \right ] + 
\left [\begin{array}{c} \frac{C_f}{m}\\\ \frac{C_fl_f}{I_z} \end{array} \right ]\delta $$

where $v_x$ and $v_y$ are the vehicle longitudinal and lateral speeds respectevely, $\psi$ is the yaw rate, $C_f$ and $C_r$ are the front and rear wheels cornering stiffness respectevely, $l_f$ and $l_r$ are the distances from the car Center of Mass to the front and rear wheels, $m$ is the vehicle mass, $I_z$ the vehicle inertia on the z-axis and $\delta$ is the steering angle. A few things to note:

* The bicycle model depends on the longitudinal speed $v_x$ which is not constant, making the model a pure LPV model. Moreover, from identification results it was found on the SAV platfrom that a cornering stiffness value set as a polynomial of $v_x$ increased the model accuracy; making the parameters $C_f$ and $C_r$ also varying.
* The model used on the simulator does not assume linear as tire forces and small angle approximatios as the bicycle model used on the MPC. Moreover, the simulator model takes into the coupling between the longitudinal and lateral dynamics.
* The simulator takes into account the actuator dynamics of the steering servo motor considering the following model: $$\frac{\delta_{real}(s)}{\delta_{command}(s)}=\frac{\omega_n^2}{s^2+2\xi\omega_ns+\omega_n^2}e^{-t_ds}$$ 

    the time delay being $t_d=0.18$ seconds. The actuator model is not considered in the model used withing the MPC formulation.

As can be seen from these points, despite the bicycle model being known to be an accurate representation of the vehicle lateral dynamics, there are phenomenas that are not captured and unmodeled dynamics that have been ignored in the MPC model formulation. These represent important sources of uncertainty the MPC lateral controller must deal with.

Using the Bicycle model to model the car lateral dynamics and the reference signal provided by the Pure-Pursuit algorithm we can create an MPC controller using CHRONOS-MPC. The MPC problem in this example was defined as follows:

$$\min_{u,x}J = \sum_{i=1}^{N-1} (r-y_i)^TQ_{e}(r-y_i) + \sum_{i=0}^{N_{ctr}-1}u_i^TR_u u_i$$

s.t.

$$ x^+=A(v_x)x+B(v_x)u$$
$$ u \in [-0.7,\\\ 0.7] $$

with

$$ y = Cx=[0 \ 1]\left [\begin{array}{c}v_y\\\ \dot \psi \end{array} \right ] $$


### References

[1] Medero Borrell, A. (2023). LPV lateral control of autonomous and automated vehicles. Universitat Politècnica de Catalunya. https://doi.org/https://dx.doi.org/10.5821/dissertation-2117-412168

[2] Ariel M. Borrell, Vicenç Puig, Olivier Sename, Fixed-structure parameter-dependent state feedback controller: A scaled autonomous vehicle path-tracking application, Control Engineering Practice, https://doi.org/10.1016/j.conengprac.2024.105911.


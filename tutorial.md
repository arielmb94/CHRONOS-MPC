# Tutorial: Implementing LPV-MPC in CHRONOS

This tutorial will guide you through the process of setting up and simulating a Linear Parameter-Varying Model Predictive Controller (LPV-MPC) using CHRONOS. We will use the **Two Tank** system as an intuitive working example.

> **Note:** The code snippets provided in this tutorial cover the core setup. You can find the complete, ready-to-run scripts for this and other examples in our "examples" folder repository.

---

### 1. The Nonlinear Model & LPV Representation

The Two Tank system consists of two connected water tanks. The input is the water flow into the first tank, and the states are the water heights $h_1$ and $h_2$. The nonlinear dynamics are governed by Torricelli's law:

$$ \dot{h}_1 = \frac{1}{A_b} u - \frac{\sqrt{2g}}{A_b} \sqrt{h_1} $$
$$ \dot{h}_2 = \frac{\sqrt{2g}}{A_b} \sqrt{h_1} - \frac{\sqrt{2g}}{A_b} \sqrt{h_2} $$

To use the CHRONOS LPV solver, we must embed the nonlinearities inside bounded scheduling parameters ($\rho$). By rewriting $\sqrt{h}$ as $h \cdot (1/\sqrt{h})$, we extract our scheduling variables:

* $\rho_1 = \frac{1}{\sqrt{h_1}}$
* $\rho_2 = \frac{1}{\sqrt{h_2}}$

Now, the system can be written in the pseudo-linear affine form $\dot{x} = A(\rho)x + Bu$, where the matrix $A$ varies dynamically depending on the current water levels.

---

### 2. The Initialization Script (`init.m`)

The initialization script is responsible for building the base CHRONOS MPC structure. You must initialize the system with the dynamic matrices "frozen" at the starting condition. 

First, we define the physical parameters, the sampling time, and the initial conditions. Then, we initialize the main `mpc` structure by defining the prediction horizon $N$:

```matlab
%% Parameters
% Tank area and gravity value
Ab = 1;
g = 9.81;

% Initial condition
h1 = 0.45;
h2 = 0.45;

% Sampling time
Ts = 0.01;

%% Create MPC object
N = 10;         % Prediction Horizon
mpc = init_mpc(N);
```

Next, we evaluate the continuous-time dynamics at the exact initial states. This provides the solver with a valid baseline LTI system. Do not worry about the matrix varying yet; the simulation loop will handle the dynamic parameter updates:

```matlab
%% LTI system
% Get LPV model frozen at current state vector
A = [-sqrt(2*g)*sqrt(h1)/(Ab*h1) 0;
     sqrt(2*g)*sqrt(h1)/(Ab*h1) -sqrt(2*g)*sqrt(h2)/(Ab*h2)];
B = [1/Ab; 0];
```

With the continuous matrices ready, we apply a Forward Euler discretization and inject them into the CHRONOS structure using `init_mpc_dynamics`. We also specify the output matrix $C$ to define which variable the MPC should track (in this case, the water height of the second tank):

```matlab
% Initialize system dynamics
% System discretized with forward Euler discretization:
% x+ = (I+Ts*A)*x+Ts*B*u+Ts*Bd*d
mpc = init_mpc_dynamics(mpc, eye(2)+Ts*A, Ts*B, 0);

% Tracking objective is the water height on the second tank
C = [0 1];
mpc = init_mpc_output(mpc, C, [], []);
```

Finally, we set up the standard MPC ingredients. We define the physical bounds for the states and control inputs, and initialize the terminal costs using the discrete Linear Quadratic Regulator (dLQR) method to guarantee stability:

```matlab
%% Constraints and Costs (Standard CHRONOS setup)
% State constraints
x_min = 0.01*ones(mpc.nx, 1);
x_max = 1*ones(mpc.nx, 1);
mpc = init_mpc_state_cnstr(mpc, x_min, x_max);

% Control input constraints
u_min = 0*ones(mpc.nu, 1);
u_max = 10*ones(mpc.nu, 1);
mpc = init_mpc_u_cnstr(mpc, u_min, u_max);

% Terminal ingredients
Qx = diag([30 30]);         
Ru = 1;                     
[mpc] = init_mpc_ter_ingredients_dlqr(mpc, Qx, Ru, 0);
```

---

### 3. The Simulation Script (`sim_lpv.m`)

In your main simulation loop, the LPV framework evaluates the scheduling variables and rebuilds the dynamic matrices at each step. 

#### 3.1. Common Setup (Required for all methods)
Regardless of the estimation method you choose, you must define the base scheduling function and the interface handle that constructs the affine matrix for a single time step.

```matlab
n = mpc.N;
n_rho = 2;

% 1. Evaluate rho based on current states
my_sched_fun = @(x) [1/sqrt(max(x(1), 1e-4)); 
                     1/sqrt(max(x(2), 1e-4))];

% 2. Build the discrete A matrix for a given rho vector
k_tank = sqrt(2*g)/Ab;
compute_A = @(rho) eye(2) + Ts * [-k_tank*rho(1), 0; 
                                   k_tank*rho(1), -k_tank*rho(2)];
```
> **Note:** Notice that in this example, we only created a `compute_A` function because the qLPV behavior is embedded entirely within the $A$ matrix. Depending on your system's dynamics, you might also need to create analogous `compute_B` and/or `compute_Bd` functions (and you might not need `compute_A` either).

At the core of the simulation loop, the workflow is always: **Estimate Trajectory ($P_k$) $\rightarrow$ Build 3D Matrices $\rightarrow$ Update Solver**.

```matlab
    % ... inside the simulation loop ...
    
    % Step 1: Estimate scheduling trajectory Pk (Choose a method below)
    % Pk = ... 
    
    % Step 2: Use the traj_mat interface to generate the 3D array
    A_lpv = traj_mat(compute_A, Pk, n_rho, n);
    % Same for B_lpv and/or Bd_lpv if necessary
    
    % Step 3: Update CHRONOS problem dynamics and solve
    mpc = update_mpc_dynamics(mpc, A_lpv, [], []);  % Note that, in our case, B_lpv and Bd_lpv are empty arrays
    % If you have B_lpv and/or Bd_lpv, you might call the function as: update_mpc_dynamics(mpc, A_lpv, B_lpv, Bd_lpv);
    [u_prev, iter, mpc] = mpc_solve(mpc, x_prev, u_prev, xf, x_ref, [], [], []);
```

#### 3.2. Choosing a Trajectory Estimation Method
CHRONOS provides four distinct algorithms to compute $P_k$. You must select the one that fits your computational budget and nonlinear complexity.

##### Method 1: Constant (Frozen) Approach
The simplest method. It freezes the current measured state and assumes the parameters will remain constant over the entire prediction horizon. No extra setup is required.
```matlab
    Pk = compute_schedul_frozen(mpc, x_prev, my_sched_fun);
```

##### Method 2: Warm-Start Evaluation (Iterative Fast)
Exploits the predicted state trajectory from the previous optimal MPC solution. It evaluates the nonlinear function along the predicted horizon in a single shot. Highly accurate and cheap.
```matlab
    Pk = compute_schedul_iterative_fast(mpc, x_prev, my_sched_fun, n_rho);
```

##### Method 3: Successive Refinement (SQP-like)
Iteratively solves the MPC problem to refine the trajectory until convergence (or maximum iterations). It yields the highest accuracy but demands more computational power. Requires declaring `n_iter`.
```matlab
% Setup before the loop:
n_iter = 5;

% Inside the loop:
Pk = compute_schedul_iterative(mpc, x_prev, u_prev, xf, x_ref, [], my_sched_fun, compute_A, [], [], n_rho, n_iter);
% If you have compute_B and/or compute_Bd, place it right after compute_A, as follows:
% Pk = compute_schedul_iterative(mpc, x_prev, u_prev, xf, x_ref, [], my_sched_fun, compute_A, compute_B, compute_Bd, n_rho, n_iter);

```

##### Method 4: Taylor Extrapolation (Recursive)
Uses a first-order Taylor expansion to project the scheduling trajectory using the analytical Jacobian. To prevent numerical divergence, it requires strict physical bounds.
```matlab
% Setup before the loop:
% Define the Jacobian matrix function
my_jacob_fun = @(x) [-0.5/(max(x(1), 1e-4)^1.5), 0; 
                      0, -0.5/(max(x(2), 1e-4)^1.5)];

% Define the operational physical bounds for clipping
rho_min = my_sched_fun([1.0; 1.0]); 
rho_max = my_sched_fun([0.01; 0.01]); 

% Inside the loop:
Pk = compute_schedul_recursive(mpc, x_prev, n_rho, my_sched_fun, my_jacob_fun, rho_min, rho_max);
```

---

### Which method should I use?
In practical terms, we encourage you to start your implementation using the **Constant (Frozen) Approach** or the **Warm-Start Evaluation (Iterative Fast)** method. They offer an exceptional balance, providing sufficiently high prediction accuracy with the computational cost of a single standard QP solve, despite being easy to set up. If your system is highly nonlinear, operates under strong disturbances, or experiences massive reference steps that invalidate the warm-start, upgrade to the **Successive Refinement (SQP)** method to ensure robust performance. The **Taylor Extrapolation (Recursive)** is also a good option when your system is continuously differentiable, with bounded derivatives, providing high accuracy with low computational overhead.
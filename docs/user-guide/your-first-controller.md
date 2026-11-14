# Your first controller

[Solving and updating the MPC](solving-and-updating-the-mpc.md) described the
online control loop. This final chapter puts that loop into practice in a
complete example, using the modeling methods from
[From nonlinear to LPV models](from-nonlinear-problems-to-lpv-mpc.md) and the
initialization workflow from
[Configuring the MPC controller](configuring-the-mpc-controller.md).

We will define the two-tank control problem, configure its costs and
constraints in CHRONOS, and use the resulting controller in a closed-loop
simulation. The complete script at the end brings these steps together.

```{contents}
:local:
:depth: 2
```

(first-controller-problem)=
## 1. The control problem

We want to regulate the water level in the second tank of a two-tank system. The inlet flow $u$ feeds tank 1. Water leaving tank 1 enters tank 2, which also drains. With a common tank area $A_b$ and gravitational acceleration $g$, the nonlinear continuous model is

$$
\begin{aligned}
\dot x_1 &= \frac{u}{A_b} - \frac{\sqrt{2gx_1}}{A_b},\\
\dot x_2 &= \frac{\sqrt{2gx_1}}{A_b} - \frac{\sqrt{2gx_2}}{A_b}.
\end{aligned}
$$

Here $x_1$ and $x_2$ are the water heights. Assume both heights are measured or estimated at each sampling instant. We regulate the water level by adjusting the inlet flow into tank 1.

The system constraints limit the water level in each tank, the inlet flow, and its rate of change. We also constrain the difference between the two water levels so that they remain close even during transients.

For related examples, see the repository's
[First Step tutorial](https://github.com/arielmb94/CHRONOS-MPC/tree/docs/tutorials/01%20-%20First%20Step)
and [Two Tank example](https://github.com/arielmb94/CHRONOS-MPC/tree/docs/examples/Two%20Tank).

### The MPC problem setting

Over a prediction horizon of $N$ steps, we choose the inlet flow to track
the target level in tank 2 while limiting changes in flow. This gives the
following MPC problem:

$$
\min_{x,u}\; J=\frac12\sum_{j=1}^{N}Q_e(r_j-y_j)^2
+\frac12\sum_{j=0}^{N-1}R_{du}(\Delta u_j)^2
$$

subject to

$$
x_{j+1}=A_jx_j+Bu_j,\qquad j=0,\ldots,N-1,
$$

$$
\begin{aligned}
0.01&\le x_j\le\mathbf{1},
&&j=1,\ldots,N,\\
-0.2&\le x_{1,j}-x_{2,j}\le0.2,
&&j=1,\ldots,N,\\
0&\le u_j\le10, &&j=0,\ldots,N-1,\\
-0.3&\le\Delta u_j\le0.3, &&j=0,\ldots,N-1.
\end{aligned}
$$

The state and output are

$$
x=\begin{bmatrix}x_1\\x_2\end{bmatrix},
\qquad y=Cx=x_2,
\qquad C=\begin{bmatrix}0&1\end{bmatrix}.
$$

There are two states, one control input, and one tracked output. Tank 1
remains part of the prediction model and constraints, while the tracking
cost acts only on tank 2. At a steady state, equal connecting and outlet
flows require $x_1=x_2$. No known-input or custom-cost channels are needed.

#### The cost ingredients

The objective has two terms. The tracking cost penalizes the difference
between the target water level $r_j$ and the predicted level in tank 2,
$y_j=x_{2,j}$. We select this output with $C=[0\;1]$ and set its weight
with $Q_e$.

The control-rate cost penalizes changes in inlet flow,
$\Delta u_j=u_j-u_{j-1}$, with weight $R_{du}$. For the first action,
$\Delta u_0=u_0-u_{\mathrm{prev}}$, where $u_{\mathrm{prev}}$ is the
previously applied flow. Together, the two weights balance tracking the
target against making smoother changes in flow. Below, we use $Q_e=200$
and $R_{du}=1$.

#### The LPV prediction model

Using the algebraic manipulation method from [From nonlinear to LPV models](from-nonlinear-problems-to-lpv-mpc.md), we multiply each term by $x_i/x_i$ to factor the nonlinear water-level dynamics:

$$
\alpha_i(x_i)=\frac{\sqrt{2gx_i}}{A_bx_i},
\qquad
A_c(x)=\begin{bmatrix}
-\alpha_1&0\\
\alpha_1&-\alpha_2
\end{bmatrix},
\qquad
B_c=\begin{bmatrix}1/A_b\\0\end{bmatrix}.
$$

Then $\dot x=A_c(x)x+B_cu$ reproduces the plant equations when the
coefficients and state use the same heights.
CHRONOS uses a discrete prediction model. With a sampling period $T_s$,
forward Euler gives

$$
x_{j+1}=Ax_j+Bu_j,
\qquad A=I_2+T_sA_c(x),
\qquad B=T_sB_c.
$$

The index $j$ denotes a stage within the prediction horizon. For
initialization, a frozen LPV model is sufficient: we evaluate $A_c$ at the
initial heights and repeat it over all $N$ transitions. In the simulation,
we instead schedule the continuous state matrices across the horizon using
the current measurement and the previously predicted states. The continuous
dynamics API performs the forward-Euler conversion at initialization and
each update. The resulting discrete matrices remain fixed during each solve.

#### Constraints

We configure three types of limits directly through the CHRONOS API:

- **State limits:** both water levels have the same operating range,
  $0.01\le x_{1,j},x_{2,j}\le1$.
- **Control limits:** the inlet flow must remain between 0 and 10.
- **Control-rate limits:** the flow may increase or decrease by at most
  0.3 per sample, so $-0.3\le\Delta u_j\le0.3$.

The remaining constraint couples the two states: their water levels should
differ by no more than 0.2. We express this as a custom constraint signal,

$$
h_j=x_{1,j}-x_{2,j}.
$$

In terms of the states, this is

$$
h_j=C_hx_j = \begin{bmatrix}1&-1\end{bmatrix}
\begin{bmatrix}x_{1,j}\\ x_{2,j}\end{bmatrix}
$$

Thus, we pass $C_h=[1\;-1]$ and the bounds $h_{\min}=-0.2$ and
$h_{\max}=0.2$ to the custom-constraint initializer. Here, $h$ denotes the
height difference; the individual tank levels are $x_1$ and $x_2$.

(first-controller-configuration)=
## 2. Configuring the MPC

We now translate the formulation into CHRONOS ingredients. Define the
dynamics first, add the costs and constraints, and build the controller
once before entering the online loop. Each step below identifies the
corresponding part of the mathematical problem.

### Initialization

Start MATLAB in your local CHRONOS repository folder, which contains
`mpc_definition_api` and `solver`. Add both folders and their subfolders to
the MATLAB path before running the controller code:

```matlab
addpath(genpath('mpc_definition_api'));
addpath(genpath('solver'));
```

Start by defining the model parameters, an initial state guess, the sampling period, and the prediction horizon length. Use [`init_chronos_mpc`](../api/init/init_chronos_mpc.md) to create the `mpc` structure that holds the controller data.

```matlab
% Plant parameters
Ab = 1;
g = 9.81;

% Starting state guess
x1 = 0.45;
x2 = 0.45;

% Sampling period
Ts = 0.01;

% Prediction horizon
N = 25;

mpc = init_chronos_mpc(N);
```

### Define the prediction dynamics

Evaluate $A_c(x)$ at the initial heights and pass the continuous coefficients
to [`init_mpc_dyn_continuous_forward`](../api/init/init_mpc_dyn_continuous_forward.md).
It applies the forward-Euler equations above using `Ts`.

```matlab
alpha1 = sqrt(2*g*x1)/(Ab*x1);
alpha2 = sqrt(2*g*x2)/(Ab*x2);
A_continuous = [-alpha1 0; alpha1 -alpha2];
B_continuous = [1/Ab; 0];

mpc = init_mpc_dyn_continuous_forward(mpc, A_continuous, ...
                                      B_continuous, [], Ts);
```

`A_continuous` is 2-by-2 and `B_continuous` is 2-by-1, establishing the two
states and single input. Supplying one matrix repeats it throughout the
horizon. The `[]` argument declares that the dynamics have no known input;
omitting the optional mass matrix selects $M=I$.

### Select the tracked output

The formulation uses $y=x_2$. Replace the default full-state output with
$C=[0\;1]$ using
[`init_mpc_output`](../api/init/init_mpc_output.md) after defining the
dynamics and before defining the tracking cost.

```matlab
C = [0 1];
mpc = init_mpc_output(mpc, C, [], []);
```

The row `C = [0 1]` selects the second tank level as the output to track.
The empty `D` and `Dd` arguments leave the output dependent only on the
state, as in $y=x_2$.

### Define the tracking cost and the control-rate cost

Supply $Q_e$ to
[`init_mpc_Tracking_cost`](../api/init/init_mpc_Tracking_cost.md).
The weight is scalar because the only tracked output is $x_2$.

```matlab
Qe = 200;
mpc = init_mpc_Tracking_cost(mpc, Qe);
```

Supply $R_{du}$ to
[`init_mpc_ControlRate_cost`](../api/init/init_mpc_ControlRate_cost.md)
to set the second term of $J$. A scalar is sufficient because there is
only one control input.

```matlab
Rdu = 1;
mpc = init_mpc_ControlRate_cost(mpc, Rdu);
```

### Define state and control limits

Set the operating range for both tank levels to 0.01 through 1 with
[`init_mpc_state_cnstr`](../api/init/init_mpc_state_cnstr.md). These are the
state bounds in the formulation above.

```matlab
x_min = 0.01;
x_max = 1;

mpc = init_mpc_state_cnstr(mpc, x_min, x_max);
```

We pass scalar values because both states share the same limits. CHRONOS automatically expands them to all $n_x$ state components.

The inlet flow must stay between 0 and 10. Register these hard bounds with [`init_mpc_control_cnstr`](../api/init/init_mpc_control_cnstr.md).

```matlab
mpc = init_mpc_control_cnstr(mpc, 0, 10);
```

The inlet flow may change by no more than 0.3 per sample. Supply the lower and
upper increment bounds to [`init_mpc_control_rate_cnstr`](../api/init/init_mpc_control_rate_cnstr.md).

```matlab
mpc = init_mpc_control_rate_cnstr(mpc, -0.3, 0.3);
```

### Define the height-difference constraint

Use [`init_mpc_custom_cnstr`](../api/init/init_mpc_custom_cnstr.md) to limit
the difference between the two tank levels. The row `Ch = [1 -1]` forms
$h=x_1-x_2$, and `h_min` and `h_max` set its allowed range.
We set both violation penalties, `qv_min` and `qv_max`, to 150 with the
aim of treating these limits as effectively hard in this example. The
constraint remains soft, so this penalizes violations without guaranteeing
that they are eliminated.

```matlab
Ch = [1 -1];
h_min = -0.2;
h_max = 0.2;
qv_min = 150;
qv_max = 150;

mpc = init_mpc_custom_cnstr(mpc, Ch, [], [], [], ...
                            h_min, h_max, qv_min, qv_max);
```

The empty `Dh`, `Dxuh`, and `Ddh` arguments leave only the state-dependent
term. The custom limits apply to predicted stages 1 through `N`, including
the terminal state. No `dh_in` forecast is needed because this signal has
no known-input term.


### Build the controller

All ingredients are now defined. Supply the initial measured state and
a guess for the inlet flow to
[`build_chronos_mpc`](../api/init/build_chronos_mpc.md) to prepare the
controller and all of its ingredients for its first solve.

```matlab
x_prev = [x1; x2];
u_prev = 3.7;
mpc = build_chronos_mpc(mpc, x_prev, u_prev);
```

The builder prepares the fixed solver workspace and an initial prediction. The value 3.7 is inside the input bounds. Keep the returned
`mpc` for the online loop. Changing the horizon, dimensions, or enabled
ingredients requires redefining and rebuilding the controller.

(first-controller-online)=
## 3. Online solve

At each sample, update the target and prediction model, solve the current
problem, and apply the first control action. We will repeat those steps
for five seconds while simulating the continuous nonlinear plant, then
plot the heights and actuator signals.

### Prepare the reference and result arrays

Start with a commanded height of 0.7, then change it to 0.25 after 2.5
seconds. This is the target for tank 2. Smooth the command with a
first-order reference filter before supplying it to the controller:

$$
r_f^+=r_f+\frac{T_s}{\tau}(r_{\mathrm{command}}-r_f),
\qquad \tau=0.1.
$$

The filter turns the abrupt command change into a gradual target change,
which suits the limited rate of change of the inlet flow. In the code,
`xf` holds the filtered target and starts at the initial level of tank 2.
We also allocate arrays to record the heights, inlet flow, flow increments,
and filtered target for plotting.

```matlab
tsim = 5;
Sim_samples = round(tsim/Ts);
time = (0:Sim_samples-1)*Ts;
time_state = (0:Sim_samples)*Ts;

r = 0.7*ones(1,Sim_samples);
r(time > 2.5) = 0.25;
tau = 0.1;
xf = x_prev(2);

x_log = zeros(2,Sim_samples+1);
x_log(:,1) = x_prev;
u_log = zeros(1,Sim_samples);
du_log = zeros(1,Sim_samples);
rf_log = zeros(1,Sim_samples);
```

The log arrays belong to the experiment, rather than the controller
workspace. `x_log` includes the initial state and the state after each
sample, so `time_state` has one more entry than `time`.

### Update the LPV prediction model and call the solver

We will explore two ways to update the prediction model. The frozen LPV method uses the current measurement to update the LPV matrices and holds them constant over the prediction horizon. The scheduled LPV method uses both the current state measurement and the states predicted by the previous MPC solve, stored in `mpc.x`. Before the first solve, we use the initial prediction prepared by the builder.

#### Frozen LPV method

This method uses only the current state measurement to update the matrices
and holds them constant over the horizon. We therefore pass a single
continuous state matrix to
[`update_mpc_dyn_continuous_forward`](../api/online/update_mpc_dyn_continuous_forward.md).
Only $A_c$ changes with height, but this full-model updater also requires
the constant `B_continuous` on every call. The controller has no known-input
channels, so we pass `[]` for that coefficient.

```matlab
a1 = sqrt(2*g*x_prev(1))/(Ab*x_prev(1));
a2 = sqrt(2*g*x_prev(2))/(Ab*x_prev(2));

A_continuous = [-a1 0; a1 -a2];

mpc = update_mpc_dyn_continuous_forward(mpc, A_continuous, ...
                                        B_continuous, [], Ts);
```

#### Scheduled LPV method

This method builds one continuous state matrix for each transition of the
horizon. We collect these matrices in `A_continuous` and pass them to
[`update_mpc_dyn_continuous_forward`](../api/online/update_mpc_dyn_continuous_forward.md)
with the constant `B_continuous`, `[]` for the absent known-input coefficient,
and `Ts`. The updater converts each stage to the discrete prediction model.

```matlab
A_continuous = zeros(2,2,mpc.N);

a1 = sqrt(2*g*x_prev(1))/(Ab*x_prev(1));
a2 = sqrt(2*g*x_prev(2))/(Ab*x_prev(2));

A_continuous(:,:,1) = [-a1 0; a1 -a2];

for j = 2:mpc.N

    x_k = mpc.x(:,j);

    a1 = sqrt(2*g*x_k(1))/(Ab*x_k(1));
    a2 = sqrt(2*g*x_k(2))/(Ab*x_k(2));

    A_continuous(:,:,j) = [-a1 0; a1 -a2];
end

mpc = update_mpc_dyn_continuous_forward(mpc, A_continuous, ...
                                        B_continuous, [], Ts);
```

### Call the solver

After updating the LPV prediction model, call
[`mpc_solve`](../api/online/mpc_solve.md) with the current tank levels
`x_prev`, the previously applied inlet flow `u_prev`, and the filtered
target `r_in`. The returned `u_k` is the flow to apply for this sample.
Keep the returned `mpc` for the next model update and solve.

```matlab
[u_k, mpc] = mpc_solve(mpc, x_prev, u_prev, ...
                                r_in, [], [], [], []);
```



### Complete simulation loop

We now combine the model update and solve steps in a loop over the full simulation.

We use the scheduled LPV method for the model update.

The complete loop below updates the reference, updates the model, calls [`mpc_solve`](../api/online/mpc_solve.md) and advances the plant.

```matlab
A_continuous = zeros(2,2,mpc.N);

for k = 1:Sim_samples
    % 1. Update the desired height of tank 2.
    xf = xf + (Ts/tau)*(r(k)-xf);
    r_in = xf;

    % 2. Schedule the LPV model across the horizon.
    a1 = sqrt(2*g*x_prev(1))/(Ab*x_prev(1));
    a2 = sqrt(2*g*x_prev(2))/(Ab*x_prev(2));
    A_continuous(:,:,1) = [-a1 0; a1 -a2];

    for j = 2:mpc.N

        x_k = mpc.x(:,j);

        a1 = sqrt(2*g*x_k(1))/(Ab*x_k(1));
        a2 = sqrt(2*g*x_k(2))/(Ab*x_k(2));

        A_continuous(:,:,j) = [-a1 0; a1 -a2];
    end

    mpc = update_mpc_dyn_continuous_forward(mpc, A_continuous, ...
                                            B_continuous, [], Ts);

    % 3. Solve and retain the returned controller state.
    [u_k, mpc] = mpc_solve(mpc, x_prev, u_prev, ...
                                r_in, [], [], [], []);
    u_log(k) = u_k;
    du_log(k) = u_k-u_prev;
    rf_log(k) = xf;

    % 4. Hold the action while the continuous nonlinear plant evolves.
    [~, x_ode] = ode45(@(t,x) [(u_k-sqrt(2*g*x(1)))/Ab; ...
                    (sqrt(2*g*x(1))-sqrt(2*g*x(2)))/Ab], ...
                    [0 Ts], x_prev);
    x_prev = x_ode(end,:)';
    u_prev = u_k;
    x_log(:,k+1) = x_prev;
end
```

The solver predicts `N` input actions and returns the first one as `u_k`.
The loop holds that action constant while `ode45` integrates the nonlinear
plant over `[0,Ts]` using its default tolerances. The final integrated
state becomes the measurement for the next sample.

### Plot the water heights

Compare the second tank's height with the filtered target supplied to the
solver, and plot tank 1 to see how it feeds tank 2. Plot the command as
well to see how the reference filter shapes the requested change. Use
`time_state` for the state history so that the final plant state is included.

```matlab
figure;
plot(time, r, ':k', time, rf_log, '--k', ...
     time_state, x_log(1,:), time_state, x_log(2,:));
grid on;
xlabel('Time (s)');
ylabel('Water height');
legend('Command', 'Filtered target', 'x1', 'x2');
```

Only tank 2 is tracked directly. Tank 1 supplies the flow that changes
tank 2, so its transient response can differ. Inspect both histories
against the soft operating range from 0.01 to 1.

### Plot the height difference

Compare the simulated height difference with the custom limits. This plot
shows the plant response associated with $h=x_1-x_2$, alongside the
desired lower and upper bounds.

```matlab
h_diff_log = x_log(1,:) - x_log(2,:);

figure;
plot(time_state, h_diff_log, ...
     time_state, h_min*ones(size(time_state)), ':r', ...
     time_state, h_max*ones(size(time_state)), ':r');
grid on;
xlabel('Time (s)');
ylabel('Height difference');
legend('x1-x2', 'Lower limit', 'Upper limit');
```

The custom limits apply to the predicted signal and can be relaxed at a
cost. Use this plot to inspect the actual height difference rather than
assuming it must remain strictly inside the interval.

### Plot the input and its increments

Use separate panels for the inlet flow and its per-sample changes. This
makes it easier to compare the flow with its bounds from 0 to 10 and the
increments with their bounds from -0.3 to 0.3.

```matlab
figure;
subplot(2,1,1);
stairs(time, u_log);
grid on;
ylabel('Inlet flow');

subplot(2,1,2);
stairs(time, du_log);
grid on;
xlabel('Time (s)');
ylabel('Input increment');
```

Each entry of `du_log` is the new inlet flow minus the previously applied
flow. For the first sample, that previous flow is the `u_prev` used to
build the controller; later samples use the action from the preceding step.

## Complete script

The following two code blocks bring the setup and closed-loop simulation
together. Paste them in order into one MATLAB script: run the offline
configuration once, then run the online simulation.

### Offline configuration

```matlab
% Plant and controller parameters
Ab = 1;
g = 9.81;
x1 = 0.45;
x2 = 0.45;
Ts = 0.01;
N = 25;

% Initialize the controller
mpc = init_chronos_mpc(N);

% Define the prediction dynamics
alpha1 = sqrt(2*g*x1)/(Ab*x1);
alpha2 = sqrt(2*g*x2)/(Ab*x2);
A_continuous = [-alpha1 0; alpha1 -alpha2];
B_continuous = [1/Ab; 0];
mpc = init_mpc_dyn_continuous_forward(mpc, A_continuous, ...
                                      B_continuous, [], Ts);

% Select the tracked output
C = [0 1];
mpc = init_mpc_output(mpc, C, [], []);

% Set the costs
Qe = 200;
mpc = init_mpc_Tracking_cost(mpc, Qe);
Rdu = 1;
mpc = init_mpc_ControlRate_cost(mpc, Rdu);

% Set state, input, and input-rate limits
x_min = 0.01;
x_max = 1;
mpc = init_mpc_state_cnstr(mpc, x_min, x_max);
mpc = init_mpc_control_cnstr(mpc, 0, 10);
mpc = init_mpc_control_rate_cnstr(mpc, -0.3, 0.3);

% Set the height-difference limits
Ch = [1 -1];
h_min = -0.2;
h_max = 0.2;
qv_min = 150;
qv_max = 150;
mpc = init_mpc_custom_cnstr(mpc, Ch, [], [], [], ...
                            h_min, h_max, qv_min, qv_max);

% Build the controller
x_prev = [x1; x2];
u_prev = 3.7;
mpc = build_chronos_mpc(mpc, x_prev, u_prev);
```

### Online simulation

```matlab
% Prepare the reference and result arrays
tsim = 5;
Sim_samples = round(tsim/Ts);
time = (0:Sim_samples-1)*Ts;
time_state = (0:Sim_samples)*Ts;
r = 0.7*ones(1,Sim_samples);
r(time > 2.5) = 0.25;
tau = 0.1;
xf = x_prev(2);

x_log = zeros(2,Sim_samples+1);
x_log(:,1) = x_prev;
u_log = zeros(1,Sim_samples);
du_log = zeros(1,Sim_samples);
rf_log = zeros(1,Sim_samples);
A_continuous = zeros(2,2,mpc.N);

for k = 1:Sim_samples
    % Update the filtered reference
    xf = xf + (Ts/tau)*(r(k)-xf);
    r_in = xf;

    % Schedule the LPV model across the horizon
    a1 = sqrt(2*g*x_prev(1))/(Ab*x_prev(1));
    a2 = sqrt(2*g*x_prev(2))/(Ab*x_prev(2));
    A_continuous(:,:,1) = [-a1 0; a1 -a2];

    for j = 2:mpc.N
        x_k = mpc.x(:,j);
        a1 = sqrt(2*g*x_k(1))/(Ab*x_k(1));
        a2 = sqrt(2*g*x_k(2))/(Ab*x_k(2));
        A_continuous(:,:,j) = [-a1 0; a1 -a2];
    end

    mpc = update_mpc_dyn_continuous_forward(mpc, A_continuous, ...
                                            B_continuous, [], Ts);

    % Solve and apply the first control action
    [u_k, mpc] = mpc_solve(mpc, x_prev, u_prev, ...
                           r_in, [], [], [], []);
    u_log(k) = u_k;
    du_log(k) = u_k-u_prev;
    rf_log(k) = xf;

    % Simulate the plant over one sample interval
    [~, x_ode] = ode45(@(t,x) [(u_k-sqrt(2*g*x(1)))/Ab; ...
                    (sqrt(2*g*x(1))-sqrt(2*g*x(2)))/Ab], ...
                    [0 Ts], x_prev);
    x_prev = x_ode(end,:)';
    u_prev = u_k;
    x_log(:,k+1) = x_prev;
end

% Plot the water heights
figure;
plot(time, r, ':k', time, rf_log, '--k', ...
     time_state, x_log(1,:), time_state, x_log(2,:));
grid on;
xlabel('Time (s)');
ylabel('Water height');
legend('Command', 'Filtered target', 'Tank 1', 'Tank 2');

% Plot the height difference
h_diff_log = x_log(1,:) - x_log(2,:);
figure;
plot(time_state, h_diff_log, ...
     time_state, h_min*ones(size(time_state)), ':r', ...
     time_state, h_max*ones(size(time_state)), ':r');
grid on;
xlabel('Time (s)');
ylabel('Height difference');
legend('Tank height difference', 'Lower limit', 'Upper limit');

% Plot the inlet flow and input increments
figure;
subplot(2,1,1);
stairs(time, u_log);
grid on;
ylabel('Inlet flow');
subplot(2,1,2);
stairs(time, du_log);
grid on;
xlabel('Time (s)');
ylabel('Input increment');
```

## Next step

You have now assembled a complete controller workflow: prepare an LPV model,
configure and build the MPC, then update, solve, and apply the first control
action at each sample. The simulation also records the plant response so you
can compare it with the reference and operating limits.

Use the [API reference](../api/index.md) for full function arguments when
adapting this workflow to your plant. For further applications, explore the
[repository examples](https://github.com/arielmb94/CHRONOS-MPC/tree/docs/examples).

# Your first controller

Build a controller that regulates the water heights in two connected tanks. This walkthrough takes you from the plant equations to a
closed-loop simulation, explaining both the one-time controller definition
and the operations performed at every control sample.

This simplifies the repository's Two Tank example by tracking both heights
and using one model matrix across the horizon at each sample. You can
copy the MATLAB blocks below, in order, into one script. They contain the
complete setup and simulation; you do not need to call `two_tank_init` first.

```{contents}
:local:
:depth: 2
```

## 1. Prepare MATLAB

Start MATLAB in the CHRONOS repository root, the folder containing
`mpc_definition_api`, `solver`, and `examples`. Add the controller functions
and their helpers to the path:

```matlab
addpath(genpath('mpc_definition_api'));
addpath(genpath('solver'));
```

The simulation uses `ode45` to integrate the continuous plant. No code
generation is needed to follow it.

## 2. Understand the plant and the control objective

The inlet flow `u` feeds tank 1. Water leaving tank 1 feeds tank 2, which
also drains. For the simplified model used in this example,

$$
\begin{aligned}
\dot h_1 &= \frac{u-\sqrt{2gh_1}}{A_b},\\
\dot h_2 &= \frac{\sqrt{2gh_1}-\sqrt{2gh_2}}{A_b}.
\end{aligned}
$$

Here `h1` and `h2` are water heights, `Ab` is the common tank area, and `g`
is gravitational acceleration. The example absorbs the discharge scaling
into the square-root flow terms. The inlet is the only actuator: changing
it affects tank 1 first, then tank 2 through the connecting flow.

We choose

$$
x=\begin{bmatrix}h_1\\h_2\end{bmatrix},\qquad y=x,\qquad C=I_2.
$$

Both heights are tracked toward the same desired level. Assume both heights
are measured or estimated at every sample. There are two states, one input,
two tracked outputs, and no known-input
or custom-signal channels in this controller.

```matlab
Ab = 1;
g = 9.81;
h1 = 0.45;
h2 = 0.45;
Ts = 0.01;
N = 25;
```

`Ts=0.01` seconds means a new action every 10 ms. With `N=25`, each solve
looks 0.25 seconds ahead. These are the example's starting choices, not a
general tuning rule. A longer horizon sees further ahead but increases the
work per solve.

## 3. Define the prediction model offline

### Express the nonlinear dynamics as an LPV model

For positive heights, write

$$
\alpha_i(h_i)=\frac{\sqrt{2gh_i}}{A_bh_i},\qquad
A_c(h)=\begin{bmatrix}-\alpha_1&0\\\alpha_1&-\alpha_2\end{bmatrix},
\qquad B_c=\begin{bmatrix}1/A_b\\0\end{bmatrix}.
$$

Then $\dot x=A_c(h)x+B_cu$ reproduces the plant equations when the
coefficients are evaluated at the same heights as the state. This is an
LPV embedding: the nonlinear dependence is carried by the matrix
coefficients. During each optimization, the supplied coefficients are fixed.
The resulting prediction problem therefore approximates the nonlinear
trajectory; an exact embedding does not make a fixed-schedule solve an exact
nonlinear MPC solve.

The coefficients divide by height, so this model requires strictly positive
scheduling heights. The example operates away from empty tanks. Its soft
state bounds, introduced below, do not themselves guarantee this domain.

### Discretize and register the model

CHRONOS expects a discrete model. Forward Euler gives

$$x_{k+1}=A_dx_k+B_du_k,\qquad A_d=I+T_sA_c,\quad B_d=T_sB_c.$$

In the code, `B_discrete` denotes the discrete control matrix; it is not
CHRONOS's optional disturbance matrix `Bd`.

```matlab
mpc = init_chronos_mpc(N);
mpc.max_iter = 3;

alpha1 = sqrt(2*g*h1)/(Ab*h1);
alpha2 = sqrt(2*g*h2)/(Ab*h2);
A_continuous = [-alpha1 0; alpha1 -alpha2];
B_continuous = [1/Ab; 0];
A_discrete = eye(2) + Ts*A_continuous;
B_discrete = Ts*B_continuous;

mpc = init_mpc_dynamics(mpc, A_discrete, B_discrete, []);
```

Initialize dynamics first so CHRONOS knows the state and input dimensions.
The nominal matrices use the initial heights and are repeated across the
horizon. The final `[]` declares no disturbance term. Initializing dynamics
also defines the default output `y=x`, with `C=eye(2)` and no direct input
or disturbance dependence. No output initialization call is needed.

The three-iteration limit follows the example's small computational budget.
It caps Newton iterations, not elapsed time, and does not guarantee convergence
within each sample. Keep that distinction in mind when inspecting results.

References: [`init_chronos_mpc`](../api/init/init_chronos_mpc.md) and
[`init_mpc_dynamics`](../api/init/init_mpc_dynamics.md).

## 4. Choose operating limits offline

We want to keep both tanks between 0.01 and 1, keep the inlet between 0 and
10, and prevent it from changing by more than 0.1 per sample:

```matlab
mpc = init_mpc_state_cnstr(mpc, 0.01*ones(2,1), ones(2,1));
mpc = init_mpc_control_cnstr(mpc, 0, 10);
mpc = init_mpc_control_rate_cnstr(mpc, -0.1, 0.1);
```

| Choice | Purpose | CHRONOS treatment |
| --- | --- | --- |
| Bounds on both heights | Express the desired operating range | Soft bounds on predicted states `x1` through `xN` |
| Bounds on inlet flow | Represent actuator limits and exclude reverse pumping | Hard bounds on controls `u0` through `u(N-1)` |
| Bounds on input increments | Limit abrupt actuator changes | Hard bounds on `u0-u_prev` and subsequent input differences |

State bounds can be relaxed at a cost when the predicted trajectory cannot
meet them. Since no explicit relaxation penalties are supplied here, the
builder resolves their defaults from the configured tracking weight and
`mpc.qv`. Control and rate bounds have no such relaxation and must be mutually
compatible with the previous input.

The rate limits are increments per sample, without division by `Ts`. Keeping
these same numbers while changing the sampling period changes the permitted
rate of movement per second. No output bounds are needed here because both tracked heights
are already included in the state bounds.

References: [state bounds](../api/init/init_mpc_state_cnstr.md),
[control bounds](../api/init/init_mpc_control_cnstr.md), and
[control-rate bounds](../api/init/init_mpc_control_rate_cnstr.md).

## 5. Choose the costs offline

### Track both tanks and discourage unnecessary input movement

```matlab
Qe = 50*eye(2);
Rdu = 1;
mpc = init_mpc_Tracking_cost(mpc, Qe);
mpc = init_mpc_ControlRate_cost(mpc, Rdu);
```

Tracking penalizes the squared errors of both heights with equal weights.
The rate cost encourages smooth input movement even when the hard rate
limit would allow a larger change. Its weight trades that preference against
tracking performance. `Qe` is a 2-by-2 matrix for the two outputs, while
`Rdu` is scalar for the single input. Both are constant over the horizon.

We do not add a control-effort cost in this example: the controller may use
the inlet flow needed for tracking, subject to its magnitude and rate limits.
The rate cost and rate constraint cause CHRONOS to maintain the preceding
control internally; you only need to supply the last applied input at solve
time.

At a steady state, equal connecting and outlet flows require `h1=h2`.
For a requested height `r`, we therefore use `[r; r]` as the tracking target.
This walkthrough uses tracking and input-rate costs over the 25-step horizon,
with no separate terminal ingredient cost.

The configured objective, before adding state-bound relaxation penalties, is

$$
J=\frac12\sum_{j=1}^{N}(r_j-x_j)^T(50I_2)(r_j-x_j)
+\frac12\sum_{j=0}^{N-1}(\Delta u_j)^2.
$$

Here `r_j` is the two-component height target. Tracking starts at stage 1
because the current heights are already measured and the output has no
direct dependence on the current input. The same tracking weight applies
through stage `N`, including both final predicted heights.

References: [tracking cost](../api/init/init_mpc_Tracking_cost.md) and
[rate cost](../api/init/init_mpc_ControlRate_cost.md).

## 6. Build once

```matlab
x_prev = [h1; h2];
u_prev = 3.7;
mpc = build_chronos_mpc(mpc, x_prev, u_prev);
```

`x_prev` is the initial state and `u_prev` represents the inlet flow applied
before the first solve. Use the actual previous actuator value in an
application, since it determines the first input increment.

The builder allocates the fixed solver workspace, prepares cost derivatives
and constraint data, and creates an initial prediction and positive slack
values. With no optional rollout reference supplied, this setup starts the
rollout from the previous input. The initial value 3.7 lies within the input
bounds and permits a zero initial input change.

Building prepares the controller; it does not compute the first optimized
action. Keep this structure for every subsequent solve. Changing `N`, signal
dimensions, or which ingredients are enabled requires redefining and rebuilding.

Reference: [`build_chronos_mpc`](../api/init/build_chronos_mpc.md).

## 7. Prepare the online experiment

Simulate five seconds with a command of 0.7, followed by 0.25 after 2.5
seconds. Smooth that command through a first-order reference filter:

$$r_f^+=r_f+\frac{T_s}{\tau}(r_{\mathrm{command}}-r_f).$$

This avoids presenting the full jump to a controller with tight rate limits
and a small iteration budget. It shapes the requested trajectory; it does
not guarantee constraint feasibility.

```matlab
tsim = 5;
Sim_samples = round(tsim/Ts);
time = (0:Sim_samples-1)*Ts;
r = 0.7*ones(1,Sim_samples);
r(time > 2.5) = 0.25;
tau = 0.1;
xf = x_prev(2);

x_log = zeros(2,Sim_samples+1);
x_log(:,1) = x_prev;
u_log = zeros(1,Sim_samples);
du_log = zeros(1,Sim_samples);
rf_log = zeros(1,Sim_samples);
iter_log = zeros(1,Sim_samples);
```

The log arrays record the experiment; they are not part of the controller
workspace.

## 8. Update, solve, and apply at every sample

The complete loop below has four steps. The explanation immediately after
it shows how each step connects to the offline choices.

```matlab
for k = 1:Sim_samples
    % 1. Update the desired heights.
    xf = xf + (Ts/tau)*(r(k)-xf);
    r_in = [xf; xf];

    % 2. Freeze one model at the current measured heights for this horizon.
    a1 = sqrt(2*g*x_prev(1))/(Ab*x_prev(1));
    a2 = sqrt(2*g*x_prev(2))/(Ab*x_prev(2));
    A_lpv = eye(2) + Ts*[-a1 0; a1 -a2];
    mpc = update_mpc_dynamics(mpc, A_lpv, [], []);

    % 3. Solve and retain the returned controller state.
    [u_k, mpc, iter] = mpc_solve(mpc, x_prev, u_prev, ...
                                r_in, [], [], [], []);
    u_log(k) = u_k;
    du_log(k) = u_k-u_prev;
    rf_log(k) = xf;
    iter_log(k) = iter;

    % 4. Hold the action while the continuous nonlinear plant evolves.
    [~, x_ode] = ode45(@(t,h) [(u_k-sqrt(2*g*h(1)))/Ab; ...
                    (sqrt(2*g*h(1))-sqrt(2*g*h(2)))/Ab], ...
                    [0 Ts], x_prev);
    x_prev = x_ode(end,:)';
    u_prev = u_k;
    x_log(:,k+1) = x_prev;
end
```

### Supply one target for both heights

`r_in=[xf;xf]` is a 2-by-1 tracking reference. CHRONOS repeats this column
across the horizon, so this example does not provide advance knowledge of
the future command step. Pass `[]` for the terminal reference because no
terminal ingredient cost is configured.
The final three empty arguments are the unused known inputs `d_in`, `dz_in`,
and `dh_in`.

### Freeze the model across each prediction horizon

At every sample, compute one 2-by-2 `A_lpv` from the measured heights in
`x_prev`. Passing this single matrix to `update_mpc_dynamics` repeats it
across all `N` transitions. No predicted states or matrix-page loop are needed.

This is frozen LPV: the model is constant across the horizon and throughout
the solve, then reevaluated from the new measurement at the next sample.
The predicted states still evolve under that model; only its coefficients
are frozen. Thus the controller adapts between samples without forecasting
how the coefficients will change within a horizon.

Only `A` changes with height. `B_discrete` remains constant and there is no
disturbance model, so `update_mpc_dynamics(mpc,A_lpv,[],[])` retains the
existing input matrix. Costs and limits also remain unchanged. The solver
uses the new dynamics and current state automatically; no manual gradient or
constraint assembly is needed.

### Keep the solver state and apply only the first action

The solve predicts `N` actions but returns the first one as `u_k`. Keeping
the returned `mpc` preserves its updated predictions and internal iterate
for the next sample. `iter` is the number of Newton iterations, not a
convergence flag; frequently reaching three means the configured cap is
being used.

### Advance the plant independently of the prediction model

The simulation holds `u_k` constant while `ode45` integrates the nonlinear
plant over `[0,Ts]`, using default tolerances. Its final state becomes the
next measurement. This plant simulation is separate from the forward-Euler
model used inside MPC.

For a physical system, replace the `ode45` block with applying `u_k`, waiting
for the next sample, and reading or estimating the new state. Set `u_prev`
to the input actually applied, including any actuator-side modification.

References: [`update_mpc_dynamics`](../api/online/update_mpc_dynamics.md) and
[`mpc_solve`](../api/online/mpc_solve.md).

## 9. Inspect the result

```matlab
figure;
subplot(3,1,1);
plot(time, r, ':k', time, rf_log, '--k', ...
     time, x_log(1,1:end-1), time, x_log(2,1:end-1));
grid on;
ylabel('Water height');
legend('Command', 'Filtered target', 'h1', 'h2');

subplot(3,1,2);
plot(time, u_log, time, du_log);
grid on;
ylabel('Inlet flow / increment');
legend('u', 'Delta u');

subplot(3,1,3);
stairs(time, iter_log);
grid on;
xlabel('Time (s)');
ylabel('Newton iterations');
```

Compare both heights against the filtered target, since both are tracked.
Tank 1 can move differently during a transient because
it supplies the flow that changes tank 2. Inspect the input and its increments
against their configured bounds, and inspect both heights against the soft
operating range. The final plant state is retained in `x_log(:,end)`.

To explore the effect of scheduling, start a fresh run from the offline
setup and omit the online model-update call. The controller then keeps its
initial prediction model while the simulated plant remains nonlinear. This
is the same kind of comparison as the repository's LTI simulation. Compare
the resulting tracking error rather than assuming either run has converged
at every sample.

You now have the complete controller lifecycle: define the model and
preferences once, build its workspace, and repeatedly update the model,
solve, and apply one action. The [definition guide](defining-the-mpc-problem.md)
explains additional ingredients, and [Solving and updating online](online-updates.md)
shows how to change weights, limits, and other signals between solves.

The original example tracks only tank 2, uses a stage-varying LPV
schedule, and includes a terminal cost with a 10-step horizon. Its files provide a more detailed variant of this walkthrough:
[`two_tank_init.m`](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Two%20Tank/two_tank_init.m),
[`sim_two_tank_lpv.m`](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Two%20Tank/sim_two_tank_lpv.m),
[`sim_two_tank_lti.m`](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Two%20Tank/sim_two_tank_lti.m),
and the [example README](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Two%20Tank/readme.md).

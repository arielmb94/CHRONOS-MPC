# From nonlinear to LPV models

CHRONOS does not attempt to work directly with nonlinear systems. Instead, it relies on the user providing the MPC problem ingredients (dynamics, costs, constraints, etc.) using Linear Parameter Varying (LPV) or Linear Time Varying (LTV) formulations. In exchange for this effort, the tradeoff for users can be summarized in the following *contract*:

```{note}
Give CHRONOS your dynamics and cost/constraint signals in its supported
LPV/LTV forms. In return, CHRONOS assembles and solves a structured MPC
problem using implementations designed for maximum computational efficiency.

For a fixed LPV/LTV schedule, every optimum is global. Stage-varying models and custom cost/constraint signals still let you capture rich model and performance behavior.
```

This guide explains how to represent nonlinear models and signals in LPV
form. Subsequent guide pages show how to use these representations to
define an MPC problem with the CHRONOS API.

```{contents}
:local:
:depth: 2
```

## LPV models

System dynamics are generally expressed as a nonlinear ODE:

$$
\dot{x}=f(x,u,d),
$$

where $x\in\mathbb{R}^{n_x}$ is the state, $u\in\mathbb{R}^{n_u}$ is the
control input, and $d\in\mathbb{R}^{n_d}$ is a known disturbance or other
known input. 

CHRONOS does not work directly with nonlinear systems; instead, it expects the dynamics to be expressed in the standard LPV form:

$$
\dot x=A(\rho)x+B(\rho)u+B_d(\rho)d.
$$

Here $A(\rho)$, $B(\rho)$, and $B_d(\rho)$ depend on a scheduling vector
$\rho$. Through *nonlinear embedding*, nonlinear terms in the original
model are expressed through these varying coefficients.

```{note}
CHRONOS uses a fixed schedule of the varying parameter $\rho$ during the solve stage. The varying parameters could be constant across the prediction horizon, for example using the current measured state, or vary from stage to stage using previous MPC predictions or an external forecast.
```

With this schedule fixed, the dynamics are input-output linear in $x$, $u$, and $d$.
This structure lets CHRONOS use linear algebra to solve the prediction
problem, even when the original model is nonlinear.

The same idea applies to a nonlinear signal,

$$
y=f_y(x,u,d),
$$

which could be represented in the LPV form

$$
y=C(\rho)x+D(\rho)u+D_d(\rho)d.
$$


The following section shows three ways to derive LPV formulations from nonlinear models:

*   {ref}`Direct factorization <lpv-direct-factorization>`
*   {ref}`Algebraic manipulation <lpv-algebraic-manipulation>`
*   {ref}`(Local) Taylor linearization <lpv-taylor-linearization>`

## Three ways to represent nonlinear equations



(lpv-direct-factorization)=
### 1. Direct factorization

The most straightforward method is to directly factorize our model with respect to the states $x$, control inputs $u$, and known inputs $d$ whenever the model allows it. Factorizing the nonlinear model gives:

$$
\dot{x}=A(\rho)x+B(\rho)u+B_d(\rho)d,
$$

where the nonlinearity is in the coefficients as
functions of the schedule $\rho$, rather than in products among the modeled
state and control variables after the schedule is fixed.

#### Scalar Example 1

Let's assume the following scalar model:

$$
\dot x = -x\sin(x) + u
$$

By factorizing the right-hand side with respect to $x$ and $u$, we obtain the following LPV model:

$$
\dot x = \left[-\sin(x) \right]x + \left[1\right]u,\qquad \rho=\sin(x),\quad A(\rho)=-\rho,\quad B=1.
$$

Note that if we expand the LPV representation by multiplying the resulting $A(\rho)$ and $B(\rho)$ matrices by the state and control input, respectively, we get back the original nonlinear model.



#### Scalar Example 2

Let's assume the following scalar model:

$$
\dot x = -x^2 + u
$$

Now we factorize the right-hand side, writing the squared term as $x\cdot x$:

$$
\dot x = \left[-x \right]x + \left[1\right]u,\qquad \rho=x,\quad A(\rho)=-\rho,\quad B=1.
$$

As we can see, this process is quite straightforward for scalar systems and for terms with a single dependency on the states or inputs. 

#### Coupled two-state system

Consider a system with two states and coupled nonlinearities, assuming
$x>0$ and $y>0$ throughout this example:

$$
\begin{aligned}
\dot x &= -xy+u\\
\dot y &= -x^2y
\end{aligned}
$$

We are now presented with multiple choices for factorizing with respect to the states. However, in these scenarios, we must be careful and always check whether the representation we pick satisfies the controllability condition.

For example, let's consider the following LPV representation:

$$
\begin{bmatrix}
\dot x \\
\dot y
\end{bmatrix} = 
\begin{bmatrix}
0 & -x \\
0 & -x^2
\end{bmatrix}\begin{bmatrix}
x \\
y
\end{bmatrix} + 
\begin{bmatrix}
1 \\
0
\end{bmatrix} u
$$

If we expand this system by multiplying the $A(\rho)$ and $B(\rho)$ matrices by the state and control input vectors, respectively, we get back the original nonlinear model. However, if we compute the controllability matrix $\left[ B\; AB\right]$ for this system, we find that it is rank deficient.

On the other hand, the following LPV representation of the system does not present such a problem:

$$
\begin{bmatrix}
\dot x \\
\dot y
\end{bmatrix} = 
\begin{bmatrix}
0 & -x \\
-xy & 0
\end{bmatrix}\begin{bmatrix}
x \\
y
\end{bmatrix} + 
\begin{bmatrix}
1 \\
0
\end{bmatrix} u
$$

```{note}
When translating the dynamics from a nonlinear model to an LPV one, always check the controllability conditions of the LPV model you obtained. You can do the check analytically or by replacing the terms of your $A(\rho)$ and $B(\rho)$ matrices with appropriate values.
```

#### Squared distance to a reference point

The LPV modeling technique also applies to system signals. Suppose an MPC controller uses the distance from a planar position $(x,y)$ to a fixed reference point $(x_0,y_0)$. To expose the states, we can square the distance signal and expand the resulting expression.

$$
\begin{aligned}
r &= \sqrt{(x-x_0)^2+(y-y_0)^2},\\
s=r^2 &= x^2-2x_0x+x_0^2+y^2-2y_0y+y_0^2.
\end{aligned}
$$

Now factor the terms involving $x$ and $y$ and use the scheduled position $\rho=[\hat x(k)\; \hat y(k)]^T$ in the coefficients. The scheduled coordinates $[\hat x(k)\; \hat y(k)]^T$ are position estimates for each stage $k$ along the prediction horizon, while $x$ and $y$ are the position variables being optimized at that stage. The resulting LPV signal is

$$
s_{\mathrm{LPV}}=
\begin{bmatrix}\hat x(k)-2x_0 & \hat y(k)-2y_0\end{bmatrix}
\begin{bmatrix}x\\y\end{bmatrix}
+0u+[x_0^2+y_0^2]d,\qquad d=1.
$$

We can use the previous MPC state predictions for these scheduled positions. This gives an accurate estimate of the squared distance to the reference point along the prediction horizon.

The same formulation can be applied to a changing reference point, such as a moving obstacle. The offset term can instead be represented by setting $D_d=1$ and $d=x_0(k)^2+y_0(k)^2$, where $x_0(k)$ and $y_0(k)$ are the forecast coordinates of the reference point along the prediction horizon. Both the scheduled position and the reference position can now vary across the horizon.

$$
s_{\mathrm{LPV}}=
\begin{bmatrix}\hat x(k)-2x_0(k) & \hat y(k)-2y_0(k)\end{bmatrix}
\begin{bmatrix}x\\y\end{bmatrix}
+0u+[1]d,\qquad d=x_0^2(k)+y_0^2(k).
$$

(lpv-algebraic-manipulation)=
### 2. Algebraic manipulation

The nonlinear model is not always directly factorizable with respect to the states and/or control inputs. However, an exact algebraic identity may expose the state or input factor needed for a matrix representation.

#### Scalar Example

Let's assume the following scalar model:

$$
\dot x = -\sin(x) + u
$$

Here we cannot factorize the model on the states at it is contained within the $\sin(x)$ term.Instead we can manipulate the model by multipliying and dividing the term $\sin(x)$ by the state value: 

$$
\sin(x)= \frac{\sin(x)\cdot x}{x}
$$

Under this manipulation we can now factorize the nonlinear model on the state, obtaining the following LPV model:

$$
\dot x = \left[-\frac{\sin(x)}{x} \right]x + \left[1\right]u,\qquad \rho=\frac{\sin(x)}{x},\quad A(\rho)=-\rho,\quad B=1.
$$

Note that if we expand the LPV representation by multiplying the resulting $A(\rho)$ and $B(\rho)$ matrices by the state and control input, respectively, we get back the original nonlinear model. However this model ill-defined for $x=0$, for this reason, this approach is only recommended in an scenario where $x>0$, and 0 is not part of the set.

#### Two-Tank Water-Level Example

Let's look now at a slightly more complex problem, the nonlinear water-level ODE from the [two-tank example](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Two%20Tank/readme.md):

$$
\begin{aligned}
\dot h_1 &= u/A_b-\sqrt{2gh_1}/A_b,\\
\dot h_2 &= \sqrt{2gh_1}/A_b-\sqrt{2gh_2}/A_b.
\end{aligned}
$$

The states $x = [h_1\; h_2]^T$ are not factorizable as they are inside square root terms. For positive tank heights, we can do the following algebraic manipulation:

$$
\sqrt{2gh_i}=\frac{\sqrt{2gh_i}}{h_i}h_i.
\qquad h_i>0.
$$

Using this modification, we can now rewrite the same equations as the following LPV model with $\rho=[h_1\;h_2]^T$:

$$
\dot{x}=
\begin{bmatrix}
-\dfrac{\sqrt{2gh_1}}{A_bh_1} & 0\\
\dfrac{\sqrt{2gh_1}}{A_bh_1} & -\dfrac{\sqrt{2gh_2}}{A_bh_2}
\end{bmatrix}x
+\begin{bmatrix}1/A_b\\0\end{bmatrix}u.
$$

Note that the coefficient expression is undefined at $h_i=0$, something to always consider when using this approach to obtain a factorizable model.

(lpv-taylor-linearization)=
### 3. Taylor linearization

As we have seen, direct factorization can sometimes lead to non-controllable LPV representations, and algebraic manipulation might not be suitable if $\rho=0$ is part of the set of varying parameters.

In these cases, we can use local linearization of only some nonlinear
terms, e.g. those that carry the missing coupling or contain nonfactorizable nonlinearities, while the remaining dynamics can keep their exact factorization. With this selective linearization, the physical state and control coordinates do not need to change.


#### Coupled two-state system

Consider again the two-state system, with the same assumptions $x>0$ and $y>0$:

$$
\begin{aligned}
\dot x &= -xy+u\\
\dot y &= -x^2y
\end{aligned}
$$

Rather than directly factorizing the first ODE, let's do a Taylor expansion of the $xy$ term, from which we obtain:

$$
\begin{aligned}
&\left\{
\begin{aligned}
\dot x &= -x^*y^*-y^*(x-x^*)-x^*(y-y^*)+u \\
\dot y &= -x^2y 
\end{aligned} 
\right. \\
&\left\{
\begin{aligned}
\dot x &= -y^*x-x^*y + u + x^*y^* \\
\dot y &= -x^2y 
\end{aligned}
\right.
\end{aligned}
$$

By directly factorizing with respect to the states and control inputs and treating the fixed term as a known input $d=x^*y^*$, we obtain the LPV model:

$$
\begin{bmatrix}
\dot x \\
\dot y
\end{bmatrix} = 
\begin{bmatrix}
-y^* & -x^* \\
-xy & 0
\end{bmatrix}\begin{bmatrix}
x \\
y
\end{bmatrix} + 
\begin{bmatrix}
1 \\
0
\end{bmatrix} u + 
\begin{bmatrix}
1 \\
0
\end{bmatrix} d, \qquad d = x^*y^*
$$

Note that, in practice, a scheduling parameter $x^*$ coming from the Taylor expansion will be treated in the same way as a state $x$ from the direct factorization of the nonlinear model.

#### Stirring-Tank Example

Let's consider the nonlinear dynamics of the end-product concentration and temperature used in the [stirring-tank example](https://github.com/arielmb94/CHRONOS-MPC/blob/docs/examples/Stirring%20Tank%20Reactor/readme.md).

$$
\begin{aligned}
\dot{c} &= (1-c)/\theta_f-kce^{-M/v},\\
\dot{v} &= (x_f-v)/\theta_f+kce^{-M/v}-\alpha u(v-x_c),
\end{aligned}
$$

where the coolant input $u$ acts on temperature $v$, and the concentration
$c$ depends on temperature through the reaction term $-kce^{-M/v}$. A direct
factorization gives

$$
\begin{bmatrix}\dot{c}\\\dot{v}\end{bmatrix}
=
\begin{bmatrix}
-1/\theta_f-ke^{-M/v} & 0\\
ke^{-M/v} & -1/\theta_f
\end{bmatrix}
\begin{bmatrix}c\\v\end{bmatrix}
+\begin{bmatrix}0\\-\alpha(v-x_c)\end{bmatrix}u
+\begin{bmatrix}1/\theta_f\\x_f/\theta_f\end{bmatrix}d,
\qquad d=1.
$$

This factorization is exact when the coefficient schedule uses the current
temperature. However, in this LPV model, there is no path from the control input to the concentration $c$ dynamics; thus, it is not controllable.

The nonlinear model does contain this coupling through $-kce^{-M/v}$. To
expose it, let's approximate only that term with a first-order Taylor expansion
around $v^o$:

$$
-kce^{-M/v}\approx-k e^{-M/v^o}c-q(v-v^o),\qquad q=\frac{k c M e^{-M/v^o}}{(v^o)^2}.
$$

The term $-q v$ adds a temperature coefficient to the concentration
equation; the offset $q v^o$ is carried through the known-input term. With
$x=[c\;v]^T$ and $d=[1\;v^o]^T$, the resulting LPV model is

$$
\dot{x}=
\begin{bmatrix}
-1/\theta_f-ke^{-M/v^o} & -q\\
ke^{-M/v} & -1/\theta_f
\end{bmatrix}x
+\begin{bmatrix}0\\-\alpha(v-x_c)\end{bmatrix}u
+\begin{bmatrix}
1/\theta_f & q\\
x_f/\theta_f & 0
\end{bmatrix}d.
$$

The temperature-to-concentration coefficient now connects the state affected
by the control to the concentration equation. Only the term needed to expose
this path is approximated; the other terms retain their factorized form.

#### Full Taylor Expansion

More generally, a first-order expansion of the full vector field at
$(x^o,u^o)$ gives

$$
f(x,u)\approx f(x^o,u^o)+A_o(x-x^o)+B_o(u-u^o)
=A_o x+B_o u+c_o,
$$

where $A_o=\partial f/\partial x$, $B_o=\partial f/\partial u$, and

$$
c_o=f(x^o,u^o)-A_o x^o-B_o u^o.
$$

At an equilibrium, $f(x^o,u^o)=0$, so
$c_o=-A_o x^o-B_o u^o$. This affine offset can be represented by
CHRONOS's known-input term, for example $d=1$ and $B_{d,o}=c_o$. This form
keeps the original state and control coordinates rather than introducing
deviation variables. Partial and full Taylor expansions are local
approximations whose accuracy is expected near the chosen operating point.

## Discretize the model

So far, we have been using continuous-time models in this guide. However, CHRONOS, like any other MPC solver, assumes discrete prediction models. Thus, the LPV representation we obtain for the dynamics needs to be discretized.

For a continuous LPV model

$$
\dot{x}=A_c(\rho)x+B_c(\rho)u+B_{d,c}(\rho)d,
$$

the discrete prediction model is

$$
x_{k+1}=A_kx_k+B_ku_k+B_{d,k}d_k,
\qquad k=0,\ldots,N-1.
$$

Under forward Euler with sample time $T_s$, the coefficients are

$$
A_k=I+T_sA_{c,k},\qquad
B_k=T_sB_{c,k},\qquad
B_{d,k}=T_sB_{d,c,k}.
$$

Other discretization methods can be used when appropriate to the plant and
sampling assumptions.

## Next steps

In this guide, you learned how to represent your nonlinear model in an LPV
form compatible with CHRONOS. Continue to
[Defining the MPC problem](defining-the-mpc-problem.md) to learn how to use
your LPV model to set up an MPC problem in CHRONOS.

# Triple Pendulum with Angle Control

## Kinematic Chain:

The [SVG graphic](docs/images/triple-pendulum.svg) graphic depicts a 2D kinematic diagram of an actuated triple pendulum hanging from a fixed horizontal support.

- Link 1: Length $l_1$, terminates at passive joint mass $m_1$. Its position is defined by the free absolute angle $\theta_1$ measured from the vertical downward axis.
- Link 2: Length $l_2$, terminates at actuated joint mass $m_2$ (depicted as a square servo housing). Its position is defined by the free relative angle $\theta_2$ measured from the extended axis of Link 1.
- Link 3: Length $l_3$, terminates at tip mass $m_3$. Its position is defined by the prescribed time-dependent relative angle $\alpha(t)$ driven by the servo, measured from the extended axis of Link 2.

## Visual Elements:
Dashed lines represent the zero-angle reference axes. Dashed arcs with arrowheads indicate the direction and measurement of angles $\theta_1$, $\theta_2$, and $\alpha(t)$.

# Coordinate System
The system possesses two free degrees of freedom. The required coordinate system for formulating pure Ordinary Differential Equations (ODEs) and avoiding Differential-Algebraic Equations (DAEs) uses minimal generalized coordinates based on the independent angles.
- $q_1(t)$: The absolute angle $\theta_1(t)$ of the first leg relative to the vertical.
- $q_2(t)$: The relative angle $\theta_2(t)$ of the second leg relative to the first leg.
The orientation of the third leg is not a free variable. Its absolute angle from the vertical is the sum of the system angles and the prescribed servo angle: $\theta_{3_{abs}}(t) = \theta_1(t) + \theta_2(t) + \alpha(t)$.

# Mathematical Approach
Lagrangian mechanics is the necessary mathematical framework. It bypasses the calculation of internal joint reaction forces, directly yielding the ODEs.
1. Kinematics: Define the Cartesian positions $(x_i, y_i)$ of each mass $m_i$ purely in terms of $q_1$, $q_2$, and the known function $\alpha(t)$. • $$x_1 = l_1 \sin(q_1)$$ • $$y_1 = -l_1 \cos(q_1)$$ • $$x_2 = x_1 + l_2 \sin(q_1 + q_2)$$ • $$y_2 = y_1 - l_2 \cos(q_1 + q_2)$$ • $$x_3 = x_2 + l_3 \sin(q_1 + q_2 + \alpha(t))$$ • $$y_3 = y_2 - l_3 \cos(q_1 + q_2 + \alpha(t))$$
2. Derivatives: Compute the velocities $(\dot{x}_i, \dot{y}_i)$. The time derivative of $x_3$ and $y_3$ will include $\dot{\alpha}(t)$.
3. Energies: Formulate the Kinetic Energy $T = \frac{1}{2}\sum m_i (\dot{x}_i^2 + \dot{y}_i^2)$ and Potential Energy $V = \sum m_i g y_i$.
4. Lagrangian: Define $L = T - V$.
5. Euler-Lagrange Equations: Apply the operator for the two free coordinates $i \in \{1, 2\}$:  $$\frac{d}{dt}\left(\frac{\partial L}{\partial \dot{q}_i}\right) - \frac{\partial L}{\partial q_i} = 0$$

# Julia Implementation Strategy
The analytic expansion of the Euler-Lagrange equations for this specific system produces hundreds of trigonometric terms, making manual derivation highly susceptible to algebraic errors.
To write ODE equations for DifferentialEquations.jl while bypassing ModelingToolkit.jl DAE processing:
1. Symbolic Derivation: Use Symbolics.jl to define $q_1(t)$, $q_2(t)$, and $\alpha(t)$ as symbolic variables.
2. Programmatic Construction: Write Julia code to compute the Cartesian equations, energies, and the Lagrangian symbolically.
3. Differentiation: Use Symbolics.derivative to execute the Euler-Lagrange operations.
4. Mass Matrix Isolation: Rearrange the resulting equations into the standard mass-matrix form:  $$M(q, t) \ddot{q} = f(q, \dot{q}, t)$$  Where $M$ is a $2 \times 2$ inertia matrix and $f$ contains the Coriolis, centrifugal, gravitational, and $\ddot{\alpha}(t)$ terms.
5. Code Generation: Use Symbolics.build_function to automatically generate optimized, allocation-free Julia code for $M$ and $f$.
6. Integration: Pass the generated functions to a DynamicalODEProblem or standard ODEProblem in DifferentialEquations.jl.

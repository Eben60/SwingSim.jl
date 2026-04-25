# Equation of Motion: Damped Pendulum with Variable Length (Cartesian Coordinates)

This document describes the Differential-Algebraic Equation (DAE) formulation of a damped pendulum with a time-varying length $L(t)$ in Cartesian coordinates $(x, y)$. This format is perfectly suited for acausal modeling with `ModelingToolkit.jl`.

## Mathematical Formulation

Instead of a single ODE for the angle $\theta$, the system is described by Newton's laws of motion in the $x$ and $y$ directions, coupled with an algebraic constraint for the length of the pendulum. 

Let the pivot be at the origin $(0,0)$, with the $y$-axis pointing upwards.

**Equations of Motion:**
$$ \ddot{x}(t) = -\lambda(t) x(t) - c \dot{x}(t) $$
$$ \ddot{y}(t) = -\lambda(t) y(t) - c \dot{y}(t) - g $$

**Algebraic Constraint (Holonomic):**
$$ x(t)^2 + y(t)^2 = L(t)^2 $$

Here, $\lambda(t)$ is an algebraic variable (Lagrange multiplier) representing the specific tension force scaled by the length ($\lambda = T / (m L)$). The mass $m$ is factored out (effectively $m=1$) since we are modeling kinematics and specific forces.

The linear damping coefficient $c$ corresponds to the angular damping factor used in the polar coordinates model:
$$ c = \frac{1}{Q}\sqrt{\frac{g}{L_0}} $$

## Variable Definitions

| Symbol | Description | Units (SI) |
| :--- | :--- | :--- |
| $x(t)$, $y(t)$ | Cartesian position of the pendulum mass | m |
| $\dot{x}(t)$, $\dot{y}(t)$ or $v_x, v_y$ | Velocity components | m/s |
| $\ddot{x}(t)$, $\ddot{y}(t)$ | Acceleration components | m/s² |
| $L(t)$ | Instantaneous length of the pendulum | m |
| $T(t)$ | Tension force in the pendulum string | N |
| $m$ | Mass of the pendulum bob (factored out as $m=1$) | kg |
| $\lambda(t)$ | Lagrange multiplier (scaled tension) | s⁻² |
| $c$ | Linear damping coefficient per unit mass | s⁻¹ |
| $g$ | Acceleration due to gravity ($\approx 9.80665$) | m/s² |
| $Q$ | Quality Factor (dimensionless damping parameter) | - |
| $L_0$ | Reference length at which $Q$ is specified | m |

## How the Tension $\lambda(t)$ is Solved

For a detailed explanation of how $\lambda(t)$ and $T(t)$ are determined algebraically through index reduction, see:
[How the Tension is Solved](TensionIndexReduction.md)

## Notes for `ModelingToolkit.jl` Implementation

1. **DAE Index Lowering:** The constraint $x^2 + y^2 = L(t)^2$ makes this an Index-3 DAE. When implementing this system in `ModelingToolkit.jl` (via `@mtkmodel` or `ODESystem`), you must use `structural_simplify` to automatically perform Pantelides index reduction.
2. **First-Order Form:** While `ModelingToolkit` can accept second-order derivatives natively (`D(D(x))`), it is often safer and more robust to explicitly define first-order states: $\dot{x} = v_x$ and $\dot{y} = v_y$, and apply the acceleration equations to $\dot{v}_x$ and $\dot{v}_y$.
3. **Equivalence:** This Cartesian formulation perfectly reproduces the inertial term ($\frac{2\dot{L}}{L}\dot{\theta}$) found in the polar equation automatically through the constraint differentiation performed by the ModelingToolkit compiler.

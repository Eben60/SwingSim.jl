# Equation of Motion: Damped Pendulum with Variable Length

The following second-order ordinary differential equation (ODE) describes the motion of a damped pendulum where the length is a function of time $L(t)$. The damping constant is expressed through the Quality Factor $Q$ defined at a reference length $L_0$.

## Mathematical Formulation

$$\ddot{\theta}(t) + \left( \frac{2\dot{L}(t)}{L(t)} + \frac{1}{Q}\sqrt{\frac{g}{L_0}} \right) \dot{\theta}(t) + \frac{g}{L(t)}\sin(\theta(t)) = 0$$

To implement this in a first-order ODE solver (such as `DifferentialEquations.jl`), the system is decomposed into:

$$\frac{d\theta}{dt} = \omega$$
$$\frac{d\omega}{dt} = -\left( \frac{2\dot{L}(t)}{L(t)} + \frac{1}{Q}\sqrt{\frac{g}{L_0}} \right) \omega - \frac{g}{L(t)}\sin(\theta)$$

## Variable Definitions

| Symbol | Description | Units (SI) |
| :--- | :--- | :--- |
| $\theta(t)$ | Angular displacement from vertical | rad |
| $\dot{\theta}(t)$ or $\omega$ | Angular velocity | rad/s |
| $\ddot{\theta}(t)$ | Angular acceleration | rad/s² |
| $L(t)$ | Instantaneous length of the pendulum | m |
| $\dot{L}(t)$ | First time-derivative of the length | m/s |
| $g$ | Acceleration due to gravity ($\approx 9.80665$) | m/s² |
| $Q$ | Quality Factor (dimensionless damping parameter) | - |
| $L_0$ | Reference length at which $Q$ is specified | m |

## Terms Analysis

* **Inertial Change Term** ($\frac{2\dot{L}}{L}\dot{\theta}$): Accounts for the conservation of angular momentum as the moment of inertia changes with $L(t)$.
* **Dissipative Term** ($\frac{1}{Q}\sqrt{\frac{g}{L_0}}\dot{\theta}$): Represents external damping (e.g., air resistance or pivot friction) based on the system's characteristics at reference length $L_0$.
* **Restoring Term** ($\frac{g}{L}\sin\theta$): The gravitational component driving the oscillation.
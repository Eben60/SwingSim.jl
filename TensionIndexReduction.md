# How the Tension $\lambda(t)$ is Solved

The mathematical system presented for the Cartesian pendulum has 3 unknown variables ($x, y, \lambda$) and 3 equations (2 ODEs, 1 algebraic constraint). This makes it a well-posed Index-3 DAE system.

You do not need to provide $\lambda(t)$ or $T(t)$ externally. The tension is an **algebraic dependent variable** that strictly arises to enforce the length constraint. Analytically, the system solves for $\lambda$ by taking the second time derivative of the length constraint $x^2 + y^2 = L(t)^2$:

1. **First derivative:** $x\dot{x} + y\dot{y} = L\dot{L}$
2. **Second derivative:** $\dot{x}^2 + x\ddot{x} + \dot{y}^2 + y\ddot{y} = \dot{L}^2 + L\ddot{L}$

Substituting the equations for $\ddot{x}$ and $\ddot{y}$ into the second derivative, applying the constraints, and solving for $\lambda$ yields an explicit formula for the tension purely in terms of the current kinematic state and the known function $L(t)$:
$$ \lambda(t) = \frac{\dot{x}^2 + \dot{y}^2 - \dot{L}^2 - L\ddot{L} - c L \dot{L} - g y}{L^2} $$

This process (differentiating constraints to find algebraic variables) is known as **index reduction**.

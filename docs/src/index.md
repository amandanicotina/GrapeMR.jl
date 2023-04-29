# GrapeMR.jl

Documentation for GrapeMR.jl

# Optimal Control
## Description of System Dynamics
The dynamics of a system can be mathematically described by the following system of equations:
```math
\begin{cases}
    \dot{\vec{x}}(t) = \vec{f}(x(t)) \\
    x(0) = x_0 (t > 0)
\end{cases}
```
Here, $\vec{f}$ represents all of the forces acting on the system. In a control dynamics, a function $u(t)$, called the control function, is introduced to allow for manipulation of the system. The system dynamics can now be written as:
```math
\begin{cases}
    \dot{\vec{x}}(t) = \vec{f}(x(t), u(t)) \\
    x(0) = x_0 (t > 0). 
\end{cases}
```
The goal in this type of system is to find the "best set" of controls. To achieve this, a reward function $C$ is introduced. If $u^*(t)$ is the optimal control, then for all $u(t)$, the reward function satisfies:
```math
C(u^*) ≥ C(u).
```

## Optimal Control in NMR
In Nuclear Magnetic Resonance (NMR), the system dynamics are governed by the Bloch Equations, which can be expressed as follows:
```math
\begin{bmatrix}
    \dot{M_x} \\ \dot{M_y} \\ \dot{M_z}
\end{bmatrix} = 
\begin{bmatrix}
    -1/T_2 & \Delta B_0 & -\omega_y \\
    -\Delta B_0 & -1/T_2 & ω_x \\
    \omega_y & -\omega_x & 1/T_1
\end{bmatrix}
\begin{bmatrix}
    M_x \\ M_y \\ M_z
\end{bmatrix} + 
\begin{bmatrix}
    0 \\ 0 \\ M_0/T_1
\end{bmatrix}
```
In this equation, $(\omega_x, \omega_y)$ represent the RF pulses used in an NMR experiment to control the dynamics of the magnetization, which is described by the vector $\vec{M} = (M_x, M_y, M_z)$. The magnetization dynamics is also governed by the intrinsic relaxation values, $T_1$ and $T_2$, of each sample.

## Cost Function
Defining a cost function is essencial in this processes since it is based on it that the system will be optimized. In the case of contrast saturation, given two samples $a$ and $b$, the goal is to supress one signal while maximizing the other. This can be achieved using the following cost function:
```math
C^{(a>b)} = ||\vec{M}^{(b)}(t_f)|| - (M_z^{(a)}(t_f))^2
```
This cost function yields the best value of $-1$ when the signal of sample $a$ is zero and the signal of sample $b$ is maximum.

## GRAPE Algorithm 
The basic algorithm for optimization is Gradient Ascent, which can be described using the following steps:
1. Choose an initial control field (can be an arbitrary smart guess);
2. Compute the cost function $C(u)$;
3. Compute the gradient $\frac{\delta C}{\delta u}$;
4. Update the current control field using the formula: $\vec{\omega}^{(n + 1)} = \vec{\omega}^{(n)} - \gamma \frac{\delta C}{\delta u}$;
5. Repeat 2-4 until convergence.

While this is the simplest way to apply the gradient ascent algorithm, a more precise approach involves using forward and backward propagation to update the control field. This approach will be discussed in the following section.

![GRAPE optimization](../images/grape_pulse.jpeg)



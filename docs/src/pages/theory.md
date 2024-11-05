
# Magnetic systems 
In a homogeneous magnetic field $\vec{B}_0 = |B_0|\mathbf{z}$, magnetic moments align with $\vec{B}_0$, resulting in a bulk magnetization $M_0$ along the z-axis. The presence of $\vec{B}_0$ induces a resonance condition:

```math
\omega_0 = -\gamma B_0
```

where $\omega_0$ [rad/s] is the Larmor frequency, $\gamma$ [rad/s/T] is the gyromagnetic ratio, and $B_0$ [T] is the homogeneous magnetic field. When an oscillating magnetic field $B_1(t)$, perpendicular to $B_0$, is applied at the resonance frequency, it rotates the magnetization. The rotation angle is given by:

```math
\alpha = 2 \pi \gamma B_1 \Delta t
```

where $\Delta t$ is the duration of the $B_1(t)$ field. Since the gyromagnetic ratio for the proton is on the order of MHz, the $B_1(t)$ field operates at radio frequencies. This field can be represented as a combination of its x and y components: 

```math
B_1(t) = B_1x + i B_1y.
```

Once $B_1(t)$ is turned off, the system relaxes to equilibrium, generating a signal known as Free Induction Decay (FID). The relaxation parameters, $T_1$ and $T_2$, are unique to each molecule, resulting in distinct decay patterns for different molecules under the same RF pulse.

# Optimal Control
## Dynamics *vs* Controlled Dynamics
The following set of equations mathematically represent the dynamics of a system:

```math
\begin{cases}
    \dot{\vec{x}}(t) = \vec{f}(x(t)) \\
 x(0) = x_0 (t > 0)
\end{cases}
```

Here, $\vec{f}$ represents all forces acting on the system. In a controlled system, a function $u(t)$, known as the control function, is introduced to allow for dynamics manipulation. The system dynamics now become:

```math
\begin{cases}
    \dot{\vec{x}}(t) = \vec{f}(x(t), u(t)) \\
 x(0) = x_0 (t > 0). 
\end{cases}
```

Such a system aims to find the optimal set of controls. To achieve this, a reward function $C$ is defined. If $u^*(t)$ is the optimal control, then for any other control $u(t)$, the reward function satisfies:

```math
C(u^*) ≥ C(u).
```

## Optimal Control in NMR
In Nuclear Magnetic Resonance (NMR), controlled dynamics are achieved using RF pulses as control inputs. The magnetization dynamics are governed by the Bloch Equations, expressed as follows:

```math
    \frac{d}{dt} \begin{pmatrix} 1 \\ M_x \\ M_y \\ M_z \end{pmatrix} = 
        \begin{pmatrix}
            1 & 0 & 0 & 0 \\
            0 & -1/T_2 & \Delta B_0 & -\mathbf{u}_y (t) \\
            0 & -\Delta B_0 & -1/T_2 & \mathbf{u}_x (t) \\
            1/T_1 & \mathbf{u}_y (t) & -\mathbf{u}_x (t) & -1/T_1
        \end{pmatrix}
        \begin{pmatrix} 1 \\ M_x \\ M_y \\ M_z \end{pmatrix}.
```

where the magnetization vector $\vec{M} = (M_x, M_y, M_z)$ describes the system's state.

## Cost Function
Defining a cost function is crucial for the optimization. For contrast saturation, given two samples $a$ and $b$, the goal is to suppress one signal while maximizing the other. This objective can be achieved with the following cost function:

```math
C^{(a>b)} = ||\vec{M}^{(b)}(t_f)|| - (M_z^{(a)}(t_f))^2
```

This cost function achieves an optimal value of $-1$ when the signal from sample $a$ is zero and the signal from sample $b$ is maximized to 1.

# GRAPE Algorithm 
The algorithm step-by-step proceeds as follows:

1. Choose an initial control field $\mathbf{u}(t_0)$ (could be a smart guess or arbitrary).
2. Calculate forward propagation of the magnetization vector $\mathbf{M}(t_0) \rightarrow \mathbf{M}(t_n)$.
3. Calculate cost function value $C(\mathbf{M}(t_n))$.
4. Calculate backward propagation of the adjoint state $\boldsymbol{\chi}(t_f) \rightarrow \boldsymbol{\chi}(t_n)$.
5. Update each of the current control fields as follows: $\mathbf{u}^{n+1}_{i} = \mathbf{u}^{n}_{i} - \epsilon \frac{\partial C}{\partial \mathbf{u}^{n}_{i}}$, for $\epsilon \geq 0$ and $i = x, y$.
6. Repeat 2-5 until convergence.

## Forward Propagation
During forward propagation, the system's state evolves under the influence of the current control field. The Bloch Equations are solved for each time step. This step yields the system state at the final time $t_f$, which is then used to evaluate the cost function $C$.

## Backward Propagation
After the cost function is computed, its gradient is calculated and used as the initial state, known as the adjoint state. This state is then propagated "backward" in time.

## Gradient Update
The control field is updated using the results computed from the backward and forward propagation, according to the rule:
```math
    \frac{\partial C}{\partial \mathbf{u}^{n}} = \chi_{n+1}^{\intercal} I \mathbf{M}_n \Delta t 
```
Where $I = I_x, I_y$ is the matrix representation for each component of the control field. The update process repeats until the cost function reaches a minimum or predefined number of iterations.

## Step size
The step size $\epsilon$ controls how much the control field is adjusted in each iteration. A larger $\epsilon$ can speed up convergence but risks overshooting the optimal solution. A smaller $\epsilon$ results in slower convergence. Adaptive step size approaches adjust $\epsilon$ dynamically, effectively balancing speed and accuracy. These values are typically selected through initial testing or automatic search.
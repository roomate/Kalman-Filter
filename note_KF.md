# Kalman filter Note

Kalman filters are a class of algorithm for estimation under uncertainty, or propagation belief. Like often in control theory, one has a state evolution equation, and an observation equation. There is two declinations, a continuous and discrete models. I'll present both. The state and observation space are respectively noted $\mathcal{X}$ and $\mathcal{Y}$.

- The discrete model: $\mathcal{X} \equiv \mathbb{R}^n$ and $\mathcal{Y} \equiv \mathbb{R}^m$.
```math
x_{k+1} &= F_k x_k + s_k + \xi_k, \, k \geq 1, \\
y_k &= H_k x_k + \eta_k.
```

$F\_k$ and $H\_k$ are respectively evolution and observation operators. $s\_k$ is a deterministic control parameter and $\xi\_k \sim N(0, Q\_k)$, $\eta\_k \sim N(0, R\_k)$ are independant random variables. Initial state is $x\_0 = m\_0 + \xi\_0$, with $m\_0 \in \mathcal{X}$ known, and $\hat x_{0|0} = m\_0$. Also, $\forall k \in \mathbb{N}, F\_k \in \mathbb{R}^{n \times n}, H\_k \in \mathbb{R}^{m \times n}, Q\_k \in \mathbb{S}_+^{n \times n} \, \text{and} \, R_k \in \mathbb{S}_+^{m \times m}$. Kalman filter aims to propose an estimator $\hat x_{k+1}$ of $x_{k+1}$ knowing the observations $y\_1, \cdots, y\_k, y_{k+1}$ of the $\textbf{true}$ state via an iterative procedure in two steps:

At step $k$:
- 'Evolution/Prediction step'. Knowing an estimator $\hat{x}_{k|k}$ of $x\_k$, propagate it thanks to the state equation; $\hat{x}_{k+1|k} = F\_k \hat x\_k + s\_k $, and make an observation $y_{k+1}$.
- 'Update step'. Update the estimator thanks the measurement of $y_{k+1}$ via conditioning, that is, compute the law of $\hat{x}_{k+1}|_{y\_1, \cdots, y_{k+1}}$.

While the evolution step is quite trivial if all parameters are determined, assumption of gaussianity on the noises make the prediction step analytically tractable. The notation $|k$ signifies knowing $y\_1, \cdots, y\_k$. It can be easily showned that the couple $(x\_k, y\_k)$ follows a Gaussian law.

I skip the computations, and only outline the results. I recall the method though. The initial condition is $\hat x_{0|0} = m\_0$ and $C_{0|0} = Q\_0$.

Assuming established that $x_{k-1|k-1} \sim N(\hat x_{k-1|k-1}, C_{k-1|k-1})$, the prediction step implies $x_{k|k-1} \sim N(\underbrace{F\_k \hat x_{k-1|k-1} + s\_k}_{\hat x_{k | k -1}}, \underbrace{F\_k C_{k-1|k-1} F\_k^T + Q\_k}_{C_{k|k-1}})$. Also, since $y\_k \sim N(H\_k\hat x_{k-1|k-1}, F\_k C_{k-1|k-1} H\_k^T + R\_k)$, one can show that the gaussian conditioning over $y\_1, \cdots y\_k$ reads:
$$
x_{k|\{k - 1\} \cup \{y\_k\}} = x_{k|k} \sim N(\underbrace{\hat x_{k | k -1} + C_{k|k - 1}H\_k^\top S\_k^{-1}(y\_k - H_k \hat x_{k|k - 1})}_{\hat x_{k|k}}, \underbrace{C_{k|k - 1} - C_{k|k - 1}H\_k^TS\_k^{-1}H\_kC_{k|k-1}}_{C_{k|k}}),
$$
where $S\_k = H\_k C_{k|k - 1}H\_k^T + R\_k$, and that's it, everything is set-up. At the end, the parameters of the model are $(s\_k)_{k \in \mathbb{N}}$, $(F_k)_{k \in \mathbb{N}}$, $(H\_k)_{k \in \mathbb{N}}$, $(Q_k)_{k \in \mathbb{N}}$, $(R\_k)_{k \in \mathbb{N}}$, and $m_0$. The control $(s_k)_{k \in \mathbb{N}}$ could for example be used to minimize a certain criteria in an optimal control context.

If $F\_k$, $H\_k$, $Q\_k$ and $R\_k$ do not actually depend on $k$, the system is said to be autonomous.

- The continuous model
```math
\dot x(t) &= Ax(t) + Bu(t) + G\xi(t) \\
y(t) &= C x(t) + E\mu(t)
```

Reference:

[1] MAP568, Gestion des incertitudes et analyse de risque, Josselin Garnier, Ecole Polytechnique.

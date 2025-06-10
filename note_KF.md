# Kalman filter Note

Kalman filters are a class of algorithm for estimation under uncertainty, or propagation belief. Like often in control theory, one has a state evolution equation, and an observation equation. There is two declinations, a continuous and discrete models. I'll present both. The state and observation space are respectively noted $\mathcal{X}$ and $\mathcal{Y}$.

- The discrete model: $\mathcal{X} \equiv \mathbb{R}^n$ and $\mathcal{Y} \equiv \mathbb{R}^m$.
\begin{align*}
x_{k+1} &= F_k x_k + s_k + \xi_k, \, k \geq 1, \\
y_k &= H_k x_k + \eta_k.
\end{align*}

$F_k$ and $H_k$ are respectively evolution and observation operators. $s_k$ is a deterministic control parameter and $\xi_k \sim N(0, Q_k)$, $\eta_k \sim N(0, R_k)$ are independant random variables. Initial state is $x_0 = m_0 + \xi_0$, with $m_0 \in \mathcal{X}$ known, and $\hat x_{0|0} = m_0$. Also, $\forall k \in \mathbb{N}, F_k \in \mathbb{R}^{n \times n}, H_k \in \mathbb{R}^{m \times n}, Q_k \in \mathbb{S}_+^{n \times n} \, \text{and} \, R_k \in \mathbb{S}_+^{m \times m}$. Kalman filter aims to propose an estimator $\hat x_{k+1}$ of $x_{k+1}$ knowing the observations $y_1, \cdots, y_k, y_{k+1}$ of the $\textbf{true}$ state via an iterative procedure in two steps:

At step $k$:
- 'Evolution/Prediction step'. Knowing an estimator $\hat x_{k|k}$ of $x_k$, propagate it thanks to the state equation; $\hat x_{k+1|k} = F_k \hat x_k + s_k $, and make an observation $y_{k+1}$.
- 'Update step'. Update the estimator thanks the measurement of $y_{k+1}$ via conditioning, that is, compute the law of $\hat x_{k+1}|_{y_1, \cdots, y_{k+1}}$.

While the evolution step is quite trivial if all parameters are determined, assumption of gaussianity on the noises make the prediction step analytically tractable. The notation $|k$ signifies knowing $y_1, \cdots, y_k$. It can be easily showned that the couple $(x_k, y_k)$ follows a Gaussian law.

I skip the computations, and only outline the results. I recall the method though. The initial condition is $\hat x_{0|0} = m_0$ and $C_{0|0} = Q_0$.

Assuming established that $x_{k-1|k-1} \sim N(\hat x_{k-1|k-1}, C_{k-1|k-1})$, the prediction step implies $x_{k|k-1} \sim N(\underbrace{F_k \hat x_{k-1|k-1} + s_k}_{\hat x_{k | k -1}}, \underbrace{F_k C_{k-1|k-1} F_k^T + Q_k}_{C_{k|k-1}})$. Also, since $y_k \sim N(H_k\hat x_{k-1|k-1}, F_k C_{k-1|k-1} H_k^T + R_k)$, one can show that the gaussian conditioning over $y_1, \cdots y_k$ reads:
$$
x_{k|\{k - 1\} \cup \{y_k\}} = x_{k|k} \sim N(\underbrace{\hat x_{k | k -1} + C_{k|k - 1}H_k^\top S_k^{-1}(y_k - H_k \hat x_{k|k - 1})}_{\hat x_{k|k}}, \underbrace{C_{k|k - 1} - C_{k|k - 1}H_k^TS_k^{-1}H_kC_{k|k-1}}_{C_{k|k}}),
$$
where $S_k = H_k C_{k|k - 1}H_k^T + R_k$, and that's it, everything is set-up. At the end, the parameters of the model are $(s_k)_{k \in \mathbb{N}}$, $(F_k)_{k \in \mathbb{N}}$, $(H_k)_{k \in \mathbb{N}}$, $(Q_k)_{k \in \mathbb{N}}$, $(R_k)_{k \in \mathbb{N}}$, and $m_0$. The control $(s_k)_{k \in \mathbb{N}}$ could for example be used to minimize a certain criteria in an optimal control context.

If $F_k$, $H_k$, $Q_k$ and $R_k$ do not actually depend on $k$, the system is said to be autonomous.

- The continuous model
\begin{align*}
\dot x(t) &= Ax(t) + Bu(t) + G\xi(t) \\
y(t) &= C x(t) + E\mu(t)
\end{align*}

Reference:

[1] MAP568, Gestion des incertitudes et analyse de risque, Josselin Garnier, Ecole Polytechnique.
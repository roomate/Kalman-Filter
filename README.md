# Kalman Filter

## Description
The kalman filter is a famous estimation technique to build estimators of the state vector $x$. The underlying dynamics and observations schemes are known, but stochastic; perturbed by a gaussian noise with known covariance matrices. Give a look at the [note](note_KF.md) to make yourself a clear idea.

## Installation

Type the following commands in your terminal to install the project locally:
`code`
mkdir KF\_dir; cd KF\_dir; git clone git@github.com:roomate/Kalman-Filter.git;
`code`

## Running the code
You can tune the different configurations parameters via your command line interface with the flags:


## Results
The parameters are 
\begin{align*}
F &= \bigg[
\cos(\pi/3) & -\sin(\pi/3) \\
\sin(\pi/3) & \cos(\pi/3)
\bigg]
H &= \bigg[
1 & 1 \\
1 & 0 \\
\bigg]
Q = \bigg[
1 & 0.5 \\
0.5 & .4 \\
\bigg]
R = \bigg[
.2 & 0 \\
0 .4 \\
\bigg]
\end{align*}
and $m\_0 = [1, 0]$.

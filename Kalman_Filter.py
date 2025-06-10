from scipy.stats import norm as normal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.widgets import Slider
import argparse

class Kalman_Filter(object):
    def __init__(self, n, mode="autonomous", F: np.ndarray = None, H: np.ndarray = None, Q: np.ndarray = None, R: np.ndarray = None):
        """
        n (int). Size of the state vector.
        mode (str). 2 modes possible: 'autonomous' and 'non-autonomous'. In the first case, define F, H, Q and R here, since it is the same at each step. In the second case, they 
        should be defined at each step in the "Initialize' and 'Prediction' method.
        """
        self.n = n
        self.mode = mode
        if mode == 'autonomous':
            #Evolution operator
            F = np.squeeze(F)
            if len(F.shape) == 1:
                self.F = np.diag(F)
            else:
                self.F = F
            #Observation operator
            self.H = np.squeeze(H)
            #Covariance matrix of noise for observation x
            if len(Q.shape) == 1:
                self.Q = np.diag(Q)
            else:
                self.Q = Q
            self.Q_sq_r = np.linalg.cholesky(Q)
            #Covariance matrix of the noise for observation y
            if len(R.shape) == 1:
                self.R = np.diag(R)
            else:
                self.R = R
            self.R_sq_r = np.linalg.cholesky(R)
        elif mode == 'non-autonomous':
            self.F = None
            self.H = None
            self.Q = None
            self.R = None

        self.x = np.array([[]]) #Estimator of the state
        self.x_post = np.array([[0]*n]) #posterior state
        self.x_true = np.array([[]]) #True state

        self.x_error = np.array([])

        self.y = np.array([[]])

        self.m0 = None
        self.C = None #C_k|k
        self.C_post = None #C_k|k-1
        self.s = np.zeros(n)

        self.min = None
        self.true_min = None

        self.x_max = None
        self.x_true_max = None


        self.error_max = 0


    def initialize(self, m0: np.ndarray, H: np.ndarray|None = None, Q: np.ndarray|None = None, R: np.ndarray|None = None):
        assert len(m0) == self.n, "m0 should be of size n."
        self.m0 = m0

        assert not len(np.squeeze(self.x)), "x should be empty"
        assert not len(np.squeeze(self.y)), "y should be empty"


        if self.mode == 'autonomous':
            H = self.H
            Q_sq_r = self.Q_sq_r
            R_sq_r = self.R_sq_r


        elif self.mode == "non-autonomous":
            assert H is not None, "H should not be None if the system is non-autonomous"
            assert Q is not None, "Q should not be None if the system is non-autonomous"
            assert R is not None, "R should not be non if the system is non-autonomous"
            R_sq_r = np.linalg.cholesky(R)
            Q_sq_r = np.linalg.cholesky(Q)

        self.x = np.concatenate((self.x, [m0]), axis=-1)
        self.min, self.max = m0, m0

        self.x_true = np.concatenate((self.x_true, [m0 + Q_sq_r @ normal.rvs(size=self.n)]), axis=-1)
        self.true_min, self.true_max = self.x_true[-1], self.x_true[-1]

        x = self.x_true[-1]
        self.y = np.concatenate((self.y, [x @ H.T + R_sq_r @ normal.rvs(size=self.n)]), axis = -1)
        self.C = self.Q


    def Predict(self, F: np.ndarray|None = None, H: np.ndarray|None = None, s: np.ndarray|float = 0, Q: np.ndarray|None = None, R: np.ndarray|None = None):
        """
        Sample from x_k|k-1, and take the observation from it.
        """
        if isinstance(s, float):
            s = np.ones(self.n)*s

        if self.mode == 'autonomous':

            F = self.F
            Q = self.Q
            H = self.H
            R = self.R
            Q_sq_r = self.Q_sq_r
            R_sq_r = self.R_sq_r

        elif self.mode == "non-autonomous":
            assert F is not None, "F should not be None if the system is non-autonomous"
            assert H is not None, "H should not be None if the system is non-autonomous"
            assert Q is not None, "Q should not be None if the system is non-autonomous"
            assert R is not None, "R should not be None if the system is non-autonomous"

            Q_sq_r = np.linalg.cholesky(Q)
            R_sq_r = np.linalg.cholesky(R)

        x = self.x[-1]
        self.x_post = np.concatenate((self.x_post, [x @ F.T + s]), axis=0)

        x = self.x_true[-1]
        self.x_true = np.concatenate((self.x_true, [x @ F.T + s + Q_sq_r @ normal.rvs(size=self.n)]), axis=0)

        x = self.x_true[-1]

        for i in range(self.n):
            if x[i] < self.true_min[i]:
                self.true_min[i] = x[i]
            elif x[i] > self.true_max[i]:
                self.true_max[i] = x[i]

        self.y = np.concatenate((self.y, [x @ H.T + R_sq_r @ normal.rvs(size=self.n)]), axis = 0)

        self.C_post = F @ self.C @ F.T + Q

    def Update(self, F: np.ndarray|None = None, H: np.ndarray|None = None, Q: np.ndarray|None = None, R: np.ndarray|None = None):
        """
        Compute x_k|k.
        """
        if self.mode=="autonomous":
            F = self.F
            Q = self.Q
            H = self.H
            R = self.R

        elif self.mode=="Non-autonomous":
            assert F is not None, "F should not be None if the system is non-autonomous"
            assert H is not None, "H should not be None if the system is non-autonomous"
            assert Q is not None, "Q should not be None if the system is non-autonomous"
            assert R is not None, "R should not be None if the system is non-autonomous"

        S = self.H @ self.C_post @ H.T + R
        y = self.y[-1]

        #Update x
        x_post = self.x_post[-1]
        tmp = np.linalg.solve(S, y - x_post @ H.T)
        self.x = np.concatenate((self.x, [x_post + self.C_post @ H.T @ tmp]), axis=0)

        x = self.x[-1]
        for i in range(self.n):
            if x[i] < self.min[i]:
                self.min[i] = x[i]
            elif x[i] > self.max[i]:
                self.max[i] = x[i]

        #Update c
        tmp = np.linalg.solve(S, H @ self.C_post)
        self.C = self.C_post - self.C_post @ H.T @ tmp

        #Compute error
        self.x_error = np.append(self.x_error, np.linalg.norm(self.x[-1] - self.x_true[-1])/np.linalg.norm(self.x_true[-1]))
        if self.x_error[-1] > self.error_max:
            self.error_max = self.x_error[-1]

    def run_autonomous(self, n_iter: int, s):
        if s is None:
            s = self.s
        for _ in range(n_iter):
            self.Predict(s=s)
            self.Update()

    def run_non_autonomous(self, n_iter: int, F: np.ndarray|None, H: np.ndarray|None, s: np.ndarray|float = 0, Q: np.ndarray|None = None, R: np.ndarray|None = None):
        if s in None:
            s = self.s
        for _ in range(n_iter):
            self.Predict(F, H, s, Q, R)
            self.Update(F, H, Q, R)

    def run(self, n_iter: int, F: np.ndarray = None, H: np.ndarray|None = None, s: np.ndarray|float = None, Q: np.ndarray|None = None, R: np.ndarray|None = None):
        if self.mode == "autonomous":
            self.run_autonomous(n_iter, s)
        elif self.mode == "Non-autonomous":
            self.run_non_autonomous(n_iter, F, H, s, Q, R)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument("--mode", default="autonomous", action='store', help="Is the system autonomous ?", choices=["autonomous", "Non-autonomous"])
    parser.add_argument("-s_x", default=0, action='store', help="Action along x-axis")
    parser.add_argument("-s_y", default=0, action='store', help="Action along y-axis")
    parser.add_argument("-n", default=2, action='store', help="dimension of state space", choices = ("2", "3"))
    args = parser.parse_args()

    n = int(args.n)

    if args.mode == "Non-autonomous":
        raise ValueError("Non-autonomous mode not yet available.")

    s = 0
    n_iter = 0
    theta = np.pi/3

    if n == 2:
        m0 = np.array([1, 0])

        F = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        H = np.array([[1, 1], [1, 0]])
        Q = np.array([[1, 0.5], [0.5, .4]])
        R = np.array([[.2, 0], [0, .4]])
    elif n == 3:
        m0 = np.array([1, 0, .5])

        phi = np.pi/6
        F = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]]) @ np.array([[np.cos(phi), 0, -np.sin(phi)], [0, 1, 0], [np.sin(phi), 0, np.cos(phi)]])
        H = np.array([[1, 1, 0], [1, 0, 1], [0, 0, 1]])
        Q = np.array([[1, 0.5, .2], [0.5, .4, 1], [.2, .4, 1]])
        R = np.array([[.4, 0, 0], [0, .6, .6], [0, .6, .7]])

    KF = Kalman_Filter(n, mode=args.mode, F=F, H=H, Q=Q, R=R)
    KF.initialize(m0)

    fig = plt.figure(figsize=(10, 5))
    if n == 3:
        ax1 = fig.add_subplot(1, 2, 1, projection="3d")
        line, = ax1.plot([], [], [], '.', c = 'blue', label="x", linestyle = ':', linewidth = .4)
        line_true, = ax1.plot([], [], [], '.', c = 'red', label="x_True", linestyle = ':', linewidth = .4)

    elif n == 2:
        ax1 = fig.add_subplot(1, 2, 1)
        line, = ax1.plot([], [], '.', c = 'blue', label="x", linestyle = ':', linewidth = .4)
        line_true, = ax1.plot([], [], '.', c = 'red', label="x_True", linestyle = ':', linewidth = .4)

    ax1.set_title("State of the system & its estimation")
    ax1.legend()

    ax2 = fig.add_subplot(1, 2, 2)
    line_error, = ax2.plot([], [], c = 'black', label="relative error")
    ax2.set_title("||x - x_True||/||x_True||")
    ax2.set_xlabel("N° of iterations")
    ax2.legend()

    def change_bound(x: float, max: bool):
        """
        Change the bound of axis depending on its sign.
        If 'max' is set to True, then it is an upper bound. Otherwise, it is a lower bound.
        """
        c = 1.2
        if max:
            if x > 0:
                return x*1.2
            else:
                return x/1.2
        else:
            if x > 0:
                return x/1.2
            else:
                return x*1.2

    def init_func():
        global n_iter

        ax1.set_xlim(-10, 10)
        ax1.set_ylim(-10, 10)
        if KF.n == 2:
            line.set_data(*[KF.x[:,i] for i in range(KF.x.shape[1])])
            line_true.set_data(*[KF.x_true[:,i] for i in range(KF.x_true.shape[1])])
        elif KF.n == 3:
            ax1.set_zlim(-10, 10)
            line.set_data_3d(*[KF.x[:,i] for i in range(KF.x.shape[1])])
            line_true.set_data_3d(*[KF.x_true[:,i] for i in range(KF.x_true.shape[1])])

        line_error.set_data(range(len(KF.x_error)), KF.x_error)

        ax2.set_ylim(0, .6)
        ax2.set_xlim(0, 100)

        n_iter += 1

    def update_fig(*args):
        global n_iter
        KF.run(1)
        
        ax1.set_xlim(np.min((-10, change_bound(KF.true_min[0], 0), change_bound(KF.min[0], 0))), np.max((10, change_bound(KF.true_max[0], 1), change_bound(KF.max[0], 1))))
        ax1.set_ylim(np.min((-10, change_bound(KF.true_min[1], 0), change_bound(KF.min[1], 0))), np.max((10, change_bound(KF.true_max[1], 1), change_bound(KF.max[1], 1))))
        if KF.n == 3:
            ax1.set_zlim(np.min((-10, KF.true_min[2] - 1, KF.min[2] - 1)), np.max((10, KF.true_max[2] + 1, KF.max[2] + 1)))
            line.set_data_3d(*[KF.x[:,i] for i in range(KF.x.shape[1])])
            line_true.set_data_3d(*[KF.x_true[:,i] for i in range(KF.x_true.shape[1])])
        elif KF.n == 2:
            line.set_data(*[KF.x[:,i] for i in range(KF.x.shape[1])])
            line_true.set_data(*[KF.x_true[:,i] for i in range(KF.x_true.shape[1])])

        if KF.x_error.shape[0] > 100:
            KF.x_error = np.delete(KF.x_error, 0)
        line_error.set_data(range(np.maximum(0, n_iter - 100), np.maximum(len(KF.x_error), n_iter)), KF.x_error)

        ax2.set_ylim(0, np.maximum(.6, KF.error_max*1.1))
        ax2.set_xlim(np.maximum(0, n_iter - 100), np.maximum(100, n_iter))

        n_iter += 1

    #Add axes on the figure
    axamp = fig.add_axes(rect=[0.05, 0.25, 0.0225, 0.33])

    #Declare the Slider
    Noise_x = Slider(
        ax=axamp,
        label="Noise in state evolution",
        valmin=0,
        valmax=10,
        valinit=1,
        orientation="vertical"
    )

    #Callback
    def update_Noise_x(val):
        """
        Control the level of noise in the state evolution
        """
        KF.Q *= Noise_x.val
        KF.Q_sq_r *= np.sqrt(Noise_x.val)

    #Declare the callack when the slider is triggered
    Noise_x.on_changed(update_Noise_x)

    #Add axes on the figure
    axmag = fig.add_axes(rect=(0.18, 0.05, 0.25, 0.03))

    #Declare the Slider
    Noise_y = Slider(
        ax=axmag,
        label="Noise in measurement",
        valmin=0,
        valmax=10,
        valinit=1,
        orientation="horizontal"
    )

    #Callback
    def update_Noise_y(val):
        """
        Control the level of noise in the measurement.
        """
        KF.R *= Noise_y.val
        KF.R_sq_r *= np.sqrt(Noise_y.val)

    #Declare the callack when the slider is triggered
    Noise_y.on_changed(update_Noise_y)

    anim = FuncAnimation(fig, func=update_fig, init_func=init_func, interval = 300, blit=False, frames = 100)

    #Display the plot
    plt.show()
from scipy.stats import norm as normal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider

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

        self.x_min = None
        self.x_true_min = None

        self.x_max = None
        self.x_true_max = None

        self.y_min = None
        self.y_true_min = None

        self.y_max = None
        self.y_true_max = None

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
        self.x_min, self.y_min = m0
        self.x_max, self.y_max = m0

        self.x_true = np.concatenate((self.x_true, [m0 + Q_sq_r @ normal.rvs(size=self.n)]), axis=-1)
        self.x_true_min, self.y_true_min = self.x_true[-1]
        self.x_true_max, self.y_true_max = self.x_true[-1]

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
        if x[0] < self.x_true_min:
            self.x_true_min = x[0]
        elif x[0] > self.x_true_max:
            self.x_true_max = x[0]

        if x[1] < self.y_true_min:
            self.y_true_min = x[1]
        elif x[1] > self.y_true_max:
            self.y_true_max = x[1]


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
        if x[0] < self.x_min:
            self.x_min = x[0]
        elif x[0] > self.x_max:
            self.x_max = x[0]

        if x[1] < self.y_min:
            self.y_min = x[1]
        elif x[1] > self.y_max:
            self.y_max = x[1]

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
    m0 = np.array([1, 0])
    s = 0
    n_iter = 0
    theta = np.pi/3

    F = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    H = np.array([[1, 1], [1, 0]])
    Q = np.array([[1, 0.5], [0.5, .4]])
    R = np.array([[.2, 0], [0, .4]])

    KF = Kalman_Filter(2, mode="autonomous", F=F, H=H, Q=Q, R=R)
    KF.initialize(m0)

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10,5))
    line, = ax[0].plot([], [], '.', c = 'blue', label="x", linestyle = ':', linewidth = .4)
    line_true, = ax[0].plot([], [], '.', c = 'red', label="x_True", linestyle = ':', linewidth = .4)
    ax[0].set_title("State of the system & its estimation")
    ax[0].legend()

    line_error, = ax[1].plot([], [], c = 'black', label="relative error")
    ax[1].set_title("||x - x_True||/||x_True||")
    ax[1].set_xlabel("N° of iterations")
    ax[1].legend()

    def init_func():
        global n_iter

        line.set_data(KF.x[:, 0], KF.x[:, 1])
        line_true.set_data(KF.x_true[:, 0], KF.x_true[:, 1])
        line_error.set_data(range(len(KF.x_error)), KF.x_error)
        ax[0].set_xlim(-10, 10)
        ax[0].set_ylim(-10, 10)

        ax[1].set_ylim(0, .6)
        ax[1].set_xlim(0, 100)

        n_iter += 1

    def update_fig(*args):
        global n_iter
        KF.run(1)
        line.set_data(KF.x[:, 0], KF.x[:, 1])
        line_true.set_data(KF.x_true[:, 0], KF.x_true[:, 1])
        ax[0].set_xlim(np.min((-10, KF.x_true_min - 1, KF.x_min - 1)), np.max((10, KF.x_true_max + 1, KF.x_max + 1)))
        ax[0].set_ylim(np.min((-10, KF.x_true_min - 1, KF.x_min - 1)), np.max((10, KF.x_true_max + 1, KF.x_max + 1)))

        if KF.x_error.shape[0] > 100:
            KF.x_error = np.delete(KF.x_error, 0)
        line_error.set_data(range(np.maximum(0, n_iter - 100), np.maximum(len(KF.x_error), n_iter)), KF.x_error)

        ax[1].set_ylim(0, np.maximum(.6, KF.error_max + .1))
        ax[1].set_xlim(np.maximum(0, n_iter - 100), np.maximum(100, n_iter))

        n_iter += 1
    

    #Add axes on the figure
    axamp = fig.add_axes(rect=[0.05, 0.25, 0.0225, 0.33])

    #Declare the Slider
    Sx = Slider(
        ax=axamp,
        label="Control along x",
        valmin=-2,
        valmax=2,
        valinit=0,
        orientation="vertical"
    )

    #Callback
    def update_Sx(val):
        KF.s[0] = Sx.val

    #Declare the callack when the slider is triggered
    Sx.on_changed(update_Sx)

    #Add axes on the figure
    ayamp = fig.add_axes(rect=[0.0, 0.25, 0.0225, 0.33])

    #Declare the Slider
    Sy = Slider(
        ax=ayamp,
        label="Control along y",
        valmin=-2,
        valmax=2,
        valinit=0,
        orientation="vertical"
    )

    #Callback
    def update_Sy(val):
        KF.s[1] = Sy.val

    #Declare the callack when the slider is triggered
    Sy.on_changed(update_Sy)

    anim = FuncAnimation(fig, func=update_fig, init_func=init_func, interval = 300, blit=False)

    #Display the plot
    plt.show()
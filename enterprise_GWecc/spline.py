import numpy as np

def spline_coeffs(t12, s12, h12):
    t1, t2 = t12
    s1, s2 = s12
    h1, h2 = h12
    
    T = np.array([[  t1**3,   t1**2, t1, 1 ],
                  [3*t1**2, 2*t1,    1,  0 ],
                  [  t2**3,   t2**2, t2, 1 ],
                  [3*t2**2, 2*t2,    1,  0 ]])
    S = np.array([s1, h1, s2, h2])
    return np.linalg.solve(T, S)

class PTASignalSpline:
    def __init__(self, ptasignal_func, waveform_func, toas, Nspline):
        self.toas = toas
        self.tmax = np.max(toas)
        self.tmin = np.min(toas)
        self.tspan = self.tmax-self.tmin
        
        self.ptasignal_func = ptasignal_func
        self.waveform_func = waveform_func
        
        self.Nspline = Nspline
        self.ts = np.linspace(self.tmin, self.tmax, Nspline+1)
        self.ss = self.ptasignal_func(self.ts)
        self.hs = self.waveform_func(self.ts)
        
        self.coeffss = [spline_coeffs((self.ts[i], self.ts[i+1]), 
                                      (self.ss[i], self.ss[i+1]), 
                                      (self.hs[i], self.hs[i+1]))
                        for i in range(Nspline)]
        
        self.eval = np.vectorize(self.evaluate)
        
    def evaluate(self, t):
        i = int((t-self.tmin)/self.tspan*self.Nspline)
        if i==self.Nspline:
            i = self.Nspline-1
        
        A = self.coeffss[i]
        T = np.array([t**3, t**2, t, 1])
        return np.dot(A, T)
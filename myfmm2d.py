import numpy as np 
from sklearn.manifold import TSNE 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec 
plt.rcParams['figure.figsize'] = (16, 8)
plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.sans-serif'] = "Times New Roman"
import skfmm

class MyFmmDe():
    def __init__(self, velo, dx=[1, 1]):
        self.dx = dx 
        self.status = np.zeros_like(velo) # 0 alive:0, trail:1, dead:2 
        self.velo = velo 
        self.T = np.ones_like(velo) * 1e6
        self.ny, self.nx = velo.shape
        self.status[:, 0] = 2 
        self.status[0, :] = 2 
        self.status[:, -1] = 2 
        self.status[-1, :] = 2 
        self.status[:, 1] = 2 
        self.status[1, :] = 2 
        self.status[:, -2] = 2 
        self.status[-2, :] = 2 
    def get_velo(self, t, v):
        t1, t2 = t 
        #"""
        d1, d2 = self.dx 
        dd1, dd2 = d1**2, d2**2
        d12 = d1*d2
        v = 1/v
        
        sqrt = (
            dd1 * v**2 + dd2 * v**2 - (t1-t2)**2
        )
        
        t = 1/(dd1+dd2) * (dd1 * t2 + dd2 * t1 + d12 * np.sqrt(np.abs(sqrt))) 
        if t > np.max([t1, t2]) and sqrt > 0:
            t = 1/(dd1+dd2) * (dd1 * t2 + dd2 * t1 + d12 * np.sqrt(sqrt)) 
            print("MAX", t, t1, t2)
        else:
            t4 = t1 + d1 * v 
            t5 = t2 + d2 * v
            t = np.min([t5, t4])
            print("MIN", t, t1, t2)
        return t 
    def get_neighbor(self, ix, iy):
        neib = np.array([
            [0, 1], 
            [1, 0], 
            [0, -1], 
            [-1, 0]
        ])
        sel = np.array([self.nx-1, self.ny-1, 0, 0]) 
        #print(sel_idx)
        array = np.array([iy, ix]) 
        neib = neib + array 
        #neib = neib[sel_idx]
        return neib 
    def source(self, sources=[[0, 0, 0]]):
        d1, d2 = self.dx
        for source in sources:
            iy, ix, t = source 
            self.status[iy, ix] = 2 
            self.T[iy, ix] = 0.
            neighbor = self.get_neighbor(ix, iy) 
            for neib in neighbor:
                i1, i2 = neib 
                if self.status[i1, i2] == 2:
                    continue 
                #r, theta = self.coord[i1, i2, i3]
                t1 = np.min([self.T[i1, i2-1], self.T[i1, i2+1]]) 
                t2 = np.min([self.T[i1-1, i2], self.T[i1+1, i2]]) 
                v = self.velo[i1, i2]
                t = self.get_velo([t1, t2], v)
                #print(i1, i2, i3, t, r, theta, r*theta)
                self.T[i1, i2] = np.min([t, self.T[i1, i2]]) 
                self.status[i1, i2] = 1  
    def fast_marching(self):
        d1, d2 = self.dx 
        is_finished = True
        while is_finished:
            #plt.clf()
            #gs = GridSpec(1, 2) 
            #fig = plt.figure(1) 
            #ax = fig.add_subplot(gs[0, 0])
            #ax.matshow(self.status[:, :]) 
            #ax = fig.add_subplot(gs[0, 1])
            #ax.matshow(np.clip(self.T[2:-2, 2:-2], 0, 30)) 
            #print(self.T[2:5, 2:5])
            #plt.show()
            id1, id2 = np.where(self.status==1)
            if len(id1)==0:
                break 
            t = self.T[id1, id2] 
            maxid = np.argmin(t) 
            iy, ix = id1[maxid], id2[maxid]
            self.status[iy, ix] = 2 
            neighbor = self.get_neighbor(ix, iy) 
            for neib in neighbor:
                i1, i2 = neib 
                if self.status[i1, i2] == 2:
                    continue 
                #r, theta = self.coord[i1, i2, i3]
                t1 = np.min([self.T[i1, i2-1], self.T[i1, i2+1]]) 
                t2 = np.min([self.T[i1-1, i2], self.T[i1+1, i2]]) 
                v = self.velo[i1, i2]
                t = self.get_velo([t1, t2], v)
                #print(i1, i2, i3, t, r, theta, r*theta)
                self.T[i1, i2] = np.min([t, self.T[i1, i2]]) 
                self.status[i1, i2] = 1  


d1 = np.linspace(1, 6, 10)
d2 = np.linspace(1, 6, 10) 
velo = np.ones([10, 10])
c1, c2 = np.meshgrid(d1, d2)
method = MyFmmDe(velo, [0.2, 0.1])
method.source([[3, 3, 0]]) 
method.fast_marching()

X = method.T[2:-2, 2:-2]

phi = np.ones([6, 6]) 
velo = np.ones([6, 6]) 
phi[1, 1] = 0
d1, d2 = np.meshgrid(np.linspace(0, 1, 6), np.linspace(0, 1, 6))
Y = skfmm.travel_time(phi, velo, dx=[0.1, 0.2], order=1)

gs = GridSpec(2, 2) 
fig = plt.figure(2) 
ax = fig.add_subplot(gs[0, 0])
ax.contour(c1[2:-2, 2:-2], c2[2:-2, 2:-2], X)
ax = fig.add_subplot(gs[0, 1])
ax.matshow(X)
ax = fig.add_subplot(gs[1, 0])
ax.contour(d1, d2, Y)
print(X[:6, :6])
print(Y[:6, :6])
ax = fig.add_subplot(gs[1, 1])
ax.matshow(Y)
plt.show()

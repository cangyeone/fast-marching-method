import numpy as np 
from sklearn.manifold import TSNE 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.gridspec import GridSpec 
plt.rcParams['figure.figsize'] = (16, 8)
plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.sans-serif'] = "Times New Roman"
import skfmm 
class MyFmm():
    def __init__(self, velo, coord, dx=[1, 1, 1]):
        self.dx = dx 
        self.status = np.zeros_like(velo) # 0 alive:0, trail:1, dead:2 
        self.velo = velo 
        self.T = np.ones_like(velo) * 1e6
        self.coord = coord
        self.nz, self.ny, self.nx = velo.shape
        self.status[:, :, 0] = 2 
        self.status[:, 0, :] = 2 
        self.status[0, :, :] = 2
        self.status[:, :, -1] = 2 
        self.status[:, -1, :] = 2 
        self.status[-1, :, :] = 2
        self.status[:, :, 1] = 2 
        self.status[:, 1, :] = 2 
        self.status[1, :, :] = 2
        self.status[:, :, -1] = 2 
        self.status[:, -2, :] = 2 
        self.status[-2, :, :] = 2
    def get_velo(self, t, v, r, theta):
        o = theta 
        t1, t2, t3 = t 
        d1, d2, d3 = self.dx 
        r2 = r ** 2
        
        cosd3d1 = (np.cos(o) * d3 * d1) ** 2 
        cosd3d2 = (np.cos(o) * d3 * d2) ** 2 
        cosd3 = np.cos(o) * d3 
        d1d2 = (d1 * d2) ** 2 
        sqrt = (
            (cosd3d1 + cosd3d2 * r2 + d1d2) * r2 * v \
                - cosd3 ** 2 * r2 * (t1 - t2) ** 2\
                    - d1 ** 2 * (t2 - t3) ** 2 \
                       - d2 ** 2 * r2 * (t1 - t3) ** 2
        )
        if sqrt > 0:
            t = 1/(cosd3d1 + cosd3d2 * r2 + d1d2) * (cosd3d1 * t2 + cosd3d2 * r2 * t1 + d1d2 * t3 + cosd3d1 * d2 * np.sqrt(sqrt)) 
            #print("Great", t)
        else:
            t = np.min([t1+d1/v, t2+r*d2/v, t3+r*np.cos(theta)*d3/v]) 
            #print("Little", t, [t1, t2, t3], [t1+d1/v, t2+r*d2/v, t3+r*np.cos(theta)*d3/v])
        return t 
    def get_neighbor(self, ix, iy, iz):
        neib = np.array([
            [0, 0, 1], 
            [0, 1, 0], 
            [1, 0, 0], 
            [0, 0, -1], 
            [0, -1, 0],
            [-1, 0, 0]
        ])
        sel = np.array([self.nx-1, self.ny-1, self.nz-1, 0, 0, 0]) 
        sel_idx = (sel!=ix)*(sel!=iy)*(sel!=iz) 
        #print(sel_idx)
        array = np.array([iz, iy, ix]) 
        neib = neib + array 
        #neib = neib[sel_idx]
        return neib 
    def source(self, sources=[[0, 0, 0, 0]]):
        d1, d2, d3 = self.dx
        for source in sources:
            iz, iy, ix, t = source 
            self.status[iz, iy, ix] = 2 
            self.T[iz, iy, ix] = 0.
            neighbor = self.get_neighbor(ix, iy, iz) 
            for neib in neighbor:
                i1, i2, i3 = neib 
                if self.status[i1, i2, i3] == 2:
                    continue 
                r, theta = self.coord[i1, i2, i3]
                t1 = np.min([self.T[i1, i2, i3-1], self.T[i1, i2, i3+1]]) 
                t2 = np.min([self.T[i1, i2-1, i3], self.T[i1, i2+1, i3]]) 
                t3 = np.min([self.T[i1-1, i2, i3], self.T[i1+1, i2, i3]]) 
                v = self.velo[i1, i2, i3]
                t = self.get_velo([t1, t2, t3], v, r, theta)
                print(i1, i2, i3, t, r, theta, r*theta)
                self.T[i1, i2, i3] = np.min([t, self.T[i1, i2, i3]]) 
                self.status[i1, i2, i3] = 1  
    def fast_marching(self):
        d1, d2, d3 = self.dx 
        is_finished = True
        while is_finished:
            #plt.clf()
            #gs = GridSpec(1, 2) 
            #fig = plt.figure(1) 
            #ax = fig.add_subplot(gs[0, 0])
            #ax.matshow(self.status[:, :, 3]) 
            #ax = fig.add_subplot(gs[0, 1])
            #ax.matshow(np.clip(self.T[2:-2, 2:-2, 3], 0, 30)) 
            #plt.show()
            id1, id2, id3 = np.where(self.status==1)
            if len(id1)==0:
                break 
            t = self.T[id1, id2, id3] 
            maxid = np.argmin(t) 
            iz, iy, ix = id1[maxid], id2[maxid], id3[maxid]
            self.status[iz, iy, ix] = 2 
            neighbor = self.get_neighbor(ix, iy, iz) 
            for neib in neighbor:
                i1, i2, i3 = neib 
                if self.status[i1, i2, i3] == 2:
                    continue 
                r, theta = self.coord[i1, i2, i3]
                t1 = np.min([self.T[i1, i2, i3-1], self.T[i1, i2, i3+1]]) 
                t2 = np.min([self.T[i1, i2-1, i3], self.T[i1, i2+1, i3]]) 
                t3 = np.min([self.T[i1-1, i2, i3], self.T[i1+1, i2, i3]]) 
                v = self.velo[i1, i2, i3]
                t = self.get_velo([t1, t2, t3], v, r, theta)
                self.T[i1, i2, i3] = np.min([t, self.T[i1, i2, i3]]) 
                self.status[i1, i2, i3] = 1 


class MyFmmDe():
    def __init__(self, velo, coord, dx=[1, 1, 1]):
        self.dx = dx 
        self.status = np.zeros_like(velo) # 0 alive:0, trail:1, dead:2 
        self.velo = velo 
        self.T = np.ones_like(velo) * 1e6
        self.nz, self.ny, self.nx = velo.shape
        self.status[:, :, 0] = 2 
        self.status[:, 0, :] = 2 
        self.status[0, :, :] = 2
        self.status[:, :, -1] = 2 
        self.status[:, -1, :] = 2 
        self.status[-1, :, :] = 2
        self.status[:, :, 1] = 2 
        self.status[:, 1, :] = 2 
        self.status[1, :, :] = 2
        self.status[:, :, -2] = 2 
        self.status[:, -2, :] = 2 
        self.status[-2, :, :] = 2
    def get_velo(self, t, v):
        t1, t2, t3 = t 
        d1, d2, d3 = self.dx 
        dd1, dd2, dd3 = d1**2, d2**2, d3**3 
        d12, d23, d13 = dd1*dd2, dd2*dd3, dd3*dd1
        v = 1/v
        sqrt = (
            d12 * v + d13 * v + d23 * v\
                - dd1 * (t2-t3)**2 \
                - dd2 * (t1-t3)**2 \
                - dd3 * (t1-t2)**2 
                    
        )
        t = 1/(d12 + d13 + d23) * (d12 * t3 + d13 * t2 + d23 * t1 + d1 * d2 * d3 * np.sqrt(np.abs(sqrt))) 
        if t > np.max(t1, t2, t3) and sqrt > 0:
            t = t
        else:
            sqrt1 = (
                        dd1 * v**2 + dd2 * v**2 - (t1-t2)**2
                    )
            sqrt2 = (
                        dd2 * v**2 + dd3 * v**2 - (t2-t3)**2
                    )
            sqrt3 = (
                        dd3 * v**2 + dd1 * v**2 - (t3-t1)**2
                    )
        return t 
    def get_neighbor(self, ix, iy, iz):
        neib = np.array([
            [0, 0, 1], 
            [0, 1, 0], 
            [1, 0, 0], 
            [0, 0, -1], 
            [0, -1, 0],
            [-1, 0, 0]
        ])
        sel = np.array([self.nx-1, self.ny-1, self.nz-1, 0, 0, 0]) 
        sel_idx = (sel!=ix)*(sel!=iy)*(sel!=iz) 
        #print(sel_idx)
        array = np.array([iz, iy, ix]) 
        neib = neib + array 
        #neib = neib[sel_idx]
        return neib 
    def source(self, sources=[[0, 0, 0, 0]]):
        d1, d2, d3 = self.dx
        for source in sources:
            iz, iy, ix, t = source 
            self.status[iz, iy, ix] = 2 
            self.T[iz, iy, ix] = 0.
            neighbor = self.get_neighbor(ix, iy, iz) 
            for neib in neighbor:
                i1, i2, i3 = neib 
                if self.status[i1, i2, i3] == 2:
                    continue 
                #r, theta = self.coord[i1, i2, i3]
                t1 = np.min([self.T[i1, i2, i3-1], self.T[i1, i2, i3+1]]) 
                t2 = np.min([self.T[i1, i2-1, i3], self.T[i1, i2+1, i3]]) 
                t3 = np.min([self.T[i1-1, i2, i3], self.T[i1+1, i2, i3]]) 
                v = self.velo[i1, i2, i3]
                t = self.get_velo([t1, t2, t3], v)
                #print(i1, i2, i3, t, r, theta, r*theta)
                self.T[i1, i2, i3] = np.min([t, self.T[i1, i2, i3]]) 
                self.status[i1, i2, i3] = 1  
    def fast_marching(self):
        d1, d2, d3 = self.dx 
        is_finished = True
        while is_finished:
            #plt.clf()
            #gs = GridSpec(1, 2) 
            #fig = plt.figure(1) 
            #ax = fig.add_subplot(gs[0, 0])
            #ax.matshow(self.status[:, :, 3]) 
            #ax = fig.add_subplot(gs[0, 1])
            #ax.matshow(np.clip(self.T[2:-2, 2:-2, 3], 0, 30)) 
            #plt.show()
            id1, id2, id3 = np.where(self.status==1)
            if len(id1)==0:
                break 
            t = self.T[id1, id2, id3] 
            maxid = np.argmin(t) 
            iz, iy, ix = id1[maxid], id2[maxid], id3[maxid]
            self.status[iz, iy, ix] = 2 
            neighbor = self.get_neighbor(ix, iy, iz) 
            for neib in neighbor:
                i1, i2, i3 = neib 
                if self.status[i1, i2, i3] == 2:
                    continue 
                #r, theta = self.coord[i1, i2, i3]
                t1 = np.min([self.T[i1, i2, i3-1], self.T[i1, i2, i3+1]]) 
                t2 = np.min([self.T[i1, i2-1, i3], self.T[i1, i2+1, i3]]) 
                t3 = np.min([self.T[i1-1, i2, i3], self.T[i1+1, i2, i3]]) 
                v = self.velo[i1, i2, i3]
                t = self.get_velo([t1, t2, t3], v)
                self.T[i1, i2, i3] = np.min([t, self.T[i1, i2, i3]]) 
                self.status[i1, i2, i3] = 1 

d1 = np.linspace(3, 6, 30)
d2 = np.linspace(0, np.pi/4, 30) 
d3 = np.linspace(0, np.pi/4, 30) 
d1 = np.linspace(1, 6, 10)
d2 = np.linspace(1, 6, 10) 
d3 = np.linspace(1, 6, 10) 
velo = np.ones([10, 10, 10])
c1, c2, c3 = np.meshgrid(d1, d2, d3)
coord = np.concatenate([itr[..., np.newaxis] for itr in [c1, c2]], axis=3)
method = MyFmmDe(velo, coord)
method.source([[3, 3, 3, 0]]) 
method.fast_marching()

X = method.T[2:-2, 2:-2, 2:-2]

phi = np.ones([6, 6, 6]) 
velo = np.ones([6, 6, 6]) 
phi[1, 1, 1] = 0
d1, d2 = np.meshgrid(np.linspace(0, 1, 6), np.linspace(0, 1, 6))
Y = skfmm.travel_time(phi, velo, dx=1.0, order=1)
print(X.shape, Y.shape)
gs = GridSpec(2, 2) 
fig = plt.figure(2) 
ax = fig.add_subplot(gs[0, 0])
ax.contour(d1, d2, X[:, 1, :])
ax = fig.add_subplot(gs[0, 1])
ax.matshow(X[:, 1, :])
ax = fig.add_subplot(gs[1, 0])
ax.contour(d1, d2, Y[:, 1, :])
print(X[:5, 1, :5])
print(Y[:5, 1, :5])
ax = fig.add_subplot(gs[1, 1])
ax.matshow(Y[:, 1, :])
plt.show()


gs = GridSpec(1, 2) 
fig = plt.figure(2) 
ax = fig.add_subplot(gs[0, 0])
ax.contour(c1[2:-2, 2:-2, 3], c2[2:-2, 2:-2, 3], X[2:-2, 2:-2, 3])
print(X[2:-2, 2:-2, 3])
ax = fig.add_subplot(gs[0, 1])
ax.matshow(X[2:-2, 3, 2:-2])
plt.show()

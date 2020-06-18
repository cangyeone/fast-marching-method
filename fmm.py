import os 
import struct 
import numpy as np 
print("=="*90)

num = 100
velo = np.ones([num * num]) * 1.0
n1 = n2 = n3 = num 
length = num * num 
files = open("dt", "wb") 
header = struct.pack("4i6f", 0, num, num, 0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0) #
files.write(header)
files.write(velo.astype(np.float32).tostring())
files.close()

files = open("ds", "wb") 
header = struct.pack("i", 2) 
files.write(header)
header = struct.pack("3if", 1, 1, 1, 0.0)
files.write(header)
header = struct.pack("3if", 50, 50, 50, 0.0) 
files.write(header)
files.close()

os.system("./test")

infiles = open("d-io", "rb") 
buff = infiles.read(4*length) 

mat = np.frombuffer(buff, dtype=np.float32)
mat = np.reshape(mat, [n1, n2])

import matplotlib.pyplot as plt  
plt.switch_backend("agg")
plt.matshow(mat)
plt.savefig("time3.png")
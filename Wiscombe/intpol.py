import numpy as np

data = np.genfromtxt("Johnson.dat")

wl = data[:,0]; n= data[:,1]; k =data[:,2]

wli = np.linspace(.4,1.,6001)
ni = np.interp(wli,wl,n)
ki = np.interp(wli,wl,k)
datai =np.zeros((len(wli),3))
datai[:,0] = wli; datai[:,1] = ni; datai[:,2] = ki;
np.savetxt("Johnson_nk.csv",datai)
print(wli)
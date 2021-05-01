# -*- coding:utf-8 -*-
"""
Author @zzeshen
Date: 2021.04.30
"""
import matplotlib.pyplot as plt
import numpy as np

# 标题字体
font1={
    'family':'Times New Roman',
    'weight':'normal',
    'size':14
}

# label字体
font2={
    'family':'Times New Roman',
    'weight':'normal',
    'size':12
}
I0=9.9
tau1=10e-3
tau2=350e-3
n=7e-3
c=3e8/1e3
v=1.5e8/1e3
y1=0
dz=1
d=500
mu=4*np.pi*1e-7
epsilon=8.8541e-10
hmax=7500 # 最大高度
exat=10000 # 精细度
def current(I0,tau1,tau2,n,t):
    temp=(t/tau1)**n/(1+(t/tau2))**n*np.exp(-t/tau2)
    eta=np.exp(-tau1/tau2*(n*tau2/tau1)**(1/n))
    i=temp*I0/eta
    return i
tsequence=np.linspace(0,3,exat+1)
i=current(I0,tau1,tau2,n,tsequence)

plt.figure(dpi=300)
plt.plot(tsequence,i,color="#1890ff")
plt.xlabel("t (ms)",font2)
plt.ylabel("Lightning Current (Ampere)",font2)
plt.title("Base Current",font1)
plt.show()
def eInduct(v,c,y1,dz,d,epsilon,hmax,exat,t):
    for z in range(0, hmax, exat):
        r = np.sqrt(z ^ 2+d ^ 2)
        sn = d/r
        t1 = t-r/c-z/v
        if t1 < 0:
            break
        else:
            y1 = y1+(2-3*sn**2)/(c*r**2)*current(I0,tau1,tau2,n,t1)*dz
    e=y1/(2*np.pi*epsilon)
    return e
def hInduct(v,c,y1,dz,d,mu,hmax,exat,t):
    for z in range(0,hmax,exat):
        r=np.sqrt(z**2+d**2)
        t1=t-r/c-z/v
        if t1<0:
            break
        else:
            y1=y1+d/r**3*current(I0,tau1,tau2,n,t)*dz
    h=y1/(2*np.pi)*mu
    return h
# 初始化数组
e=np.zeros((exat+1))
h=np.zeros((exat+1))

i=0

for t in np.linspace(0,3,exat+1):
    e[i]=eInduct(v,c,y1,dz,d,epsilon,hmax,exat,t)
    h[i]=hInduct(v,c,y1,dz,d,mu,hmax,exat,t)
    i=i+1
fig=plt.figure(dpi=300)
ay1=fig.add_subplot(111)
ay1.plot(tsequence, e, color="#1890ff", linewidth=1.2)
ay1.set_ylabel("Electirc Field (V/m)",fontdict=font1)
ay1.set_xlabel("t (ms)")
ay2=plt.twinx(ay1)
ay2.plot(tsequence, h, color="#f5222d", linewidth=1.2)
ay2.set_ylabel("Magnetic Field (A/m)")
plt.xlabel("t (ms)",fontdict=font1)
plt.legend("","Magnetic Field")
plt.title("Lightning EMF Simulation",fontdict=font2)
plt.show()
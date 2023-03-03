# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 03:52:19 2020

@author: lian_liu
"""

import numpy as np
import matplotlib.pyplot as plt

# define T_lamda
def Compute_T_lamda(rho,h,lamda):
    DC_lamda = np.ones([lamda.shape[0], lamda.shape[1]])*rho[-1]
    T_lamda = np.ones([lamda.shape[0], lamda.shape[1]])
    for j in range(np.size(rho)-2, -1, -1):
        A=1-np.exp(-2*lamda*h[j])
        B=1+np.exp(-2*lamda*h[j])
        V=rho[j]*A/B
#        print(V.size)
        DC_lamda=(V+DC_lamda)/(1+V*DC_lamda/(rho[j]**2))
#        T=rho[j]*[rho(i)*A+T*B)/(rho(i)*B+T*A)

    T_lamda = DC_lamda * lamda
#    print(np.size(T_lamda))
    return T_lamda

# define lamda
def Compute_lamda(r,a,C_num,s):
    C=np.arange(0,C_num)
    lamda=(1/r)*(10.0**(a+C*s))
    return lamda

#1 arrange receivers, rho, h
r_num = 100
r = 10**np.linspace(0,4,r_num).reshape(-1,1)
#rho = np.array([10,100,200,100,50,100])*1000
rho = np.array([10**1.5,10**0.5,10**2.5])
#rho=np.array([10,10,10,10,10,10])*10000
# nlayer = 10
h = np.array([10,30])
# h = 100.0 * (1.05**np.arange(0, nlayer-1, 1))

# 设置滤波系数、偏移量
C_Gup=np.array([3.17926147465e-06,-9.73811660718e-06,1.64866227408e-05,-1.81501261160e-05,1.87556556369e-05,
                    -1.46550406038e-05,1.53799733803e-05,-6.95628273934e-06,1.41881555665e-05,3.41445665537e-06,
                    2.13941715512e-05,2.34962369042e-05,4.84340283290e-05,7.33732978590e-05,1.27703784430e-04,
                    2.08120025730e-04,3.49803898913e-04,5.79107814687e-04,9.65887918451e-04,1.60401273703e-03,
                    2.66903777685e-03,4.43111590040e-03,7.35631696247e-03,1.21782796293e-02,2.01097829218e-02,
                    3.30096953061e-02,5.37143591532e-02,8.60516613299e-02,1.34267607144e-01,2.00125033067e-01,
                    2.74027505792e-01,3.18168749246e-01,2.41655667461e-01,-5.40549161658e-02,-4.46912952135e-01,
                    -1.92231885629e-01,5.52376753950e-01,-3.57429049025e-01,1.41510519002e-01,-4.61421935309e-02,
                    1.48273761923e-02,-5.07479209193e-03,1.83829713749e-03,-6.67742804324e-04,2.21277518118e-04,
                    -5.66248732755e-05,7.88229202853e-06])
# C_Gup=C_Gup.reshape(1,-1)
a_Gup=-3.05078187595e+00
s_Gup=1.10599010095e-01

App_resis_Gup = np.zeros(r_num)
lamda_Gup = np.zeros([r_num, C_Gup.shape[0]])
T_lamda_Gup = np.ones([r_num, C_Gup.shape[0]])

lamda_Gup = Compute_lamda(r,a_Gup,C_Gup.shape[0],s_Gup)
T_lamda_Gup = Compute_T_lamda(rho,h,lamda_Gup)
App_resis_Gup = np.sum(T_lamda_Gup * C_Gup, 1).reshape(-1,1) * r

plt.loglog(r,App_resis_Gup, label='Guptasarma_apparent_resisitivity')
plt.legend(loc='lower right')#lower left upper right
plt.title('DC forward Schlumberger')
plt.ylabel('Value (ohm.m)')
plt.xlabel('r=AB/2 (m)')
plt.show()
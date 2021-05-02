import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.stats import norm
from numpy import *
from random import *
from math import *

a=2
k=3
alfa=0.025
gamma=0.9
n=400
count = 1000

mx=k/a
dx=k/a**2

print("mx "+str(mx)+"\ndx "+str(dx))
xi = [0 for i in range(n)]


#input = open(r"C:\Users\Ksenia\terver.txt", "w")
"""
for i in range(n):
    xi[i]=0
    for j in range(k):
        t = random()
        xi[i]=xi[i]-1/a*log(t, math.e)
#    input.write(str(xi[i])+"\n")
input.close()
"""
input = open(r"C:\Users\Ksenia\terver.txt", "r")
text = open(r"C:\Users\Ksenia\output.txt","w")

for i in range(n):
    xi[i] = double(input.readline())
input.close()

#количество интервалов 4
N = floor(1+ 3.32*log(n, 10)) + 1
print("N> "+str(N)+"\n")
#Подсчитываем частоты Vk попадания выборочных значений в k-ый интервал kJ
u=(max(xi)-min(xi))/N
text.write("u> "+str(u)+"\n")
for i in range (N+1):
    text.write(str(min(xi)+u*i)+"\n")
text.write("> "+str(max(xi))+"\n")


vk = [0]*N
for i in range(N):
    for x in xi:
        if(min(xi)+ u*i <= x <= min(xi)+u*(i+1)):
            vk[i]+=1
text.write("chastota> "+str(vk)+"\n")

pk=[0]*N
for i in range(N):
    pk[i]=vk[i]/n
text.write("otnosit chastota> "+str(pk) +"\n")

hk=[0]*N
for i in range(N):
    hk[i]=vk[i]/(u*n)
text.write("visota> "+str(hk) +"\n")


f=[0]*N
for i in range (N):
    mid = u/2*(i+1)
    f[i] = (a*(a*mid)**(k-1)*e**(-a*mid))/(math.factorial(k-1))
text.write("teor znachenie f> "+str(f) + "\n")

def func(xi):
    #global xi
    global n
    global a
    global k
    global count
    f=[0]*count
    for i in range (count):
        f[i] = (a*(a*xi[i])**(k-1)*e**(-a*xi[i]))/(math.factorial(k-1))
    return f

x = np.linspace(min(xi), max(xi), count)

plt.hist(xi, bins = N, density=True)
plt.plot(x, func(x))
plt.show()

#2.2
mx1 = 0
for i in range(n):#выборочное среднее
    mx1 += xi[i] 
mx1/=n
text.write("mx1= " + str(mx1)+"\n")


dx1 = 0
for i in range(n):#выборочная дисперсия смещенная
    dx1+=(xi[i]-mx1)**2
dx1/=n
text.write("dx1= " +str(dx1)+"\n")
#2.3

#2.4

#2.5
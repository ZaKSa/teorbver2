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


xi = [0 for i in range(n)]

for i in range(n):
    xi[i]=0
    for j in range(k):
        t = random()
        xi[i]=xi[i]-1/a*log(t, math.e)
    print(str(xi[i])+"\n")
#input = open(r"C:\Users\Ksenia\terver.txt", "r")
#for i in range(n):
#    xi[i] = double(input.readline())
#input.close()

#text = open(r"C:\Users\Ksenia\output.txt","w")

"""
This is yet another multiline comment
For long comments the quotation marks are much easier than # comments.
"""

#количество интервалов 4
#N = floor(1+ 3.32*log(n, 10)) + 1
#text.write("N> "+str(N)+"\n")
#Подсчитываем частоты Vk попадания выборочных значений в k-ый интервал kJ
"""u=(max(xi)+1-min(xi))/N
text.write("u> "+str(u)+"\n")
for i in range (N):
    text.write("> "+str(min(xi)+u*i)+"\n")
text.write(str(max(xi))+"\n")
"""

"""vk = [0]*N
for k in xi:
    l = (k-min(xi))//u
    vk[int(l)] += 1
"""
#text.write("chastota> "+str(vk)+"\n")
"""
pk=[0]*N
for i in range(N):
    pk[i]=vk[i]/n
text.write("otnosit chastota> "+str(pk) +"\n")


hk=[0]*N
for i in range(N):
    hk[i]=vk[i]/(u*n)
text.write("visota> "+str(hk) +"\n")
"""

"""
for i in range (N):
    mid = u/2*(i+1)
    f = math.exp((mid-a)**2/(2*gamma2))/(gamma*sqrt(2*math.pi))
    text.write("teor znachenie f> "+str(f) + "\n")
"""

"""
x = np.linspace(a-8, a+8, 1000)
pdf = norm.pdf(x, loc = a, scale=gamma)

plt.hist(xi, bins = N, density=True)
plt.plot(x, pdf)
plt.show()"""
#plt.hist(xi, bins=np.arange(min(xi), max(xi)+1))
#plt.show()


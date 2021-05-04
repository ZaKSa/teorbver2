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

print("mx "+str(mx)+"\ndx "+str(dx))#истинные
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

#количество интервалов 10
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

sum0=0
for i in range(N):
    sum0+=vk[i]
text.write("summ> "+str(sum0) +"\n")

pk=[0]*N
for i in range(N):
    pk[i]=vk[i]/n
text.write("otnosit chastota> "+str(pk) +"\n")

sum=0
for i in range(N):
    sum+=pk[i]
text.write("summ> "+str(sum) +"\n")

hk=[0]*N
for i in range(N):
    hk[i]=vk[i]/(u*n)
text.write("visota> "+str(hk) +"\n")


f=[0]*N
for i in range (N):
    mid = ((min(xi)+u*i)+(min(xi)+u*(i+1)))/2
    f[i] = (a*(a*mid)**(k-1)*e**(-a*mid))/(math.factorial(k-1))
text.write("teor znachenie f> "+str(f) + "\n")

def func(xi):
    global a
    global k
    global count
    f=[0]*count
    for i in range (count):
        f[i] = (a*(a*xi[i])**(k-1)*e**(-a*xi[i]))/(math.factorial(k-1))
    return f

#x = np.linspace(min(xi), max(xi), count)

#plt.hist(xi, bins = N, density=True)
#plt.plot(x, func(x))
#plt.show()

#2.2
mx1 = 0
for i in range(n):#выборочное среднее
    mx1 += xi[i] 
mx1/=n
text.write("mx1= " + str(mx1)+"\n")


dx1 = 0
for i in range(n):#выборочная дисперсия
    dx1+=(xi[i]-mx1)**2
dx1/=n
text.write("dx1= " +str(dx1)+"\n")


####
dx_no = 0
for i in range(n):#выборочная дисперсия несмещенная
    dx_no+=(xi[i]-mx1)**2
dx_no/=(n-1)
text.write("dx1 nesmeshennaya= " +str(dx_no)+"\n")


x2 = 0
for i in range(N):#
    mid = ((min(xi)+u*i)+min(xi)+u*(i+1))/2
    x2 += vk[i]*mid
x2/=n
text.write("mx2= " +str(x2)+"\n")


dx2 = 0
for i in range(N):#
    mid = ((min(xi)+u*i)+min(xi)+u*(i+1))/2
    dx2+=(mid-x2)**2*vk[i]
dx2/=n
text.write("dx2= " + str(dx2)+"\n")


#2.3 метод моментов
k_m=(mx1)**2/dx1
text.write("metod momentov> k= "+ str(k_m)+"\n")

a_m=mx1/dx1
text.write("metod momentov> a= "+ str(a_m)+"\n")
#2.4 Построить довер интервалы

#(1+y)/2=0.95
c=1.645

#приближенный дов интервал для Д
t1_5=n*dx1/(n-1+c*sqrt(2*(n-1)))
t2_5=n*dx1/(n-1-c*sqrt(2*(n-1)))
text.write("приближенный дов интервал для Д\n"+"t1_5= "+str(t1_5)+"\n")
text.write("t2_5= "+str(t2_5)+"\n")
#приближенный дов для М
t1_6=mx1 - c*dx1/sqrt(n)
t2_6=mx1 + c*dx1/sqrt(n)
text.write("приближенный дов для М\n"+"t1_6= "+str(t1_6)+"\n")
text.write("t2_6= "+str(t2_6)+"\n")
#приближенный дов для Д
M4=0
for i in range(n):
    M4+=(xi[i]-mx1)**4
M4/=n
t1_7=dx1 -c*sqrt(M4-dx1**2)/sqrt(n)
t2_7=dx1 +c*sqrt(M4-dx1**2)/sqrt(n)
text.write("приближенный дов для Д\n"+"t1_7= "+str(t1_7)+"\n")
text.write("t2_7= "+str(t2_7)+"\n")

#2.5

#1.5
#привести значения статистики критерия и порога
#порог
#1-alfa=0.975, (N-1)-1=8
hi2=17.535
#1-alfa=0.975, (N-1)-1-r=9-1-2=6
hi2_un=14.449
#объединеним последний и предпоследний интервалы, первй и второй
N_=N-1
p_un = [0]*(N_)

interv_=[0]*(N_+1)
text.write("posle > ")
interv_[0]=min(xi)
text.write(str(min(xi))+"\n")
for i in range (1,N_):
    interv_[i]=min(xi)+u*(i)
    text.write(str(min(xi)+u*(i))+"\n")
interv_[N_]=max(xi)
text.write(str(max(xi))+"\n"+"\n")
for i in range (N_+1):
    text.write(str(interv_[i])+"\n")
text.write("\n")

vk_ = [0]*N_
text.write("posle chastota> ")
for i in range(N_-1):
    vk_[i]=vk[i]
vk_[-1]=vk[-1]+vk[-2]
for i in range(N_):
    text.write(str(vk_[i])+"\n")


pk_=[0]*N_
text.write("posle otnosit chastota> ")
for i in range(N_-1):
    pk_[i]=pk[i]
pk_[-1]=pk[-1]+pk[-2]
for i in range(N_):
    text.write(str(pk_[i]) +"\n")

sum1=0
for i in range(N_):
    sum1+=vk_[i]
text.write("summ> "+str(sum1) +"\n")

sum2=0
for i in range(N_):
    sum2+=pk_[i]
text.write("summ> "+str(sum2) +"\n")

text.write("posle teor znachenie f> ")
text.write(str(norm.pdf(((min(xi))+min(xi)+u*2)/2, loc = a, scale= gamma)) + "\n")
for i in range (1, N_-1):
    mid_ = ((min(xi)+u*(i+1))+min(xi)+u*(i+2))/2
    #f_ = math.exp((mid_-a)**2/(2*gamma2))/(gamma*sqrt(2*math.pi))
    text.write(str(norm.pdf(mid_, loc = a, scale= gamma)) + "\n")
text.write(str(norm.pdf(min(xi)+u*9, loc = a, scale= gamma)) + "\n")

######################################################################
#статистика критерия
#(a*(a*u_)**(k-1)*
def F_un(u_):
    F_u=0
    global k
    global a
    for i in range(k-1):
        F_u+=(a*u_)**i/(math.factorial(i))
    F_u=1-(e**(-a*u_)*F_u)
    
    return F_u#norm.cdf(u_, loc = mx1, scale = sqrt(dx1))

#p = [0]*(N_)
hi_stat=0
hi_stat_unknown=0


for i in range(N_-1):#######################
     p_un[i]=F_un(min(xi)+u*(i+1))-F_un(min(xi)+u*(i))
     #hi_stat_unknown += (vk_[k]-n*p_un[k])**2/(n*p_un[k])
p_un[N_-1]=F_un(min(xi)+u*10)-F_un(min(xi)+u*8)

for i in range(N_):##
     hi_stat_unknown += (vk_[i]-n*p_un[i])**2/(n*p_un[i])


text.write("veroyatnost' un "+(str(p_un))+ "\n")
    
text.write(str(hi_stat_unknown)+"\n")


if hi_stat_unknown >= hi2_un:
    text.write("гипотезу H0 отвергают\n")
else:
    text.write("гипотезу H0 принимают\n")
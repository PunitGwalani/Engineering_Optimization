from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from random import *
from sys import exit

def ackley_func(x,y):
    
    f=-20.0 * exp(-0.2 * sqrt(0.5 * (x**2 + y**2))) - exp(0.5 * (cos(2 * pi * x) + cos(2 * pi * y))) + e + 20
    return f


[lb_x,ub_x]=[-100,100] #domain
[lb_y,ub_y]=[-100,100] #domain
[x,y]=[10,10] #initial position
x_step=0.1
y_step=0.1
            
fitness_og=ackley_func(x,y)

n=1000 #Number of iterations
T_init=10 #Initial temperature
T_n=T_init
fitness=np.zeros(n)
count=np.arange(0,n,1)
fit_plot=np.array([])
alpha=0.8
beta=150
gamma=0.08
gompertz=np.zeros(n)
fitness_difference=np.zeros(n)
optima_list=np.array([])
threshold=10
new_optima=10
temp=15
var=0

for i in range(n):
        if(i%10==0):
           T_n=T_n*0.95
                
        fitness_current=ackley_func(x,y)
        
        fitness[i]=fitness_current
        fitness_new=ackley_func(x+x_step,y+y_step)
        fitness_diff=fitness_current-fitness_new
        fitness_difference[i]=fitness_diff
        k=-fitness_diff/T_n
        
        gompertz[i]=1-(alpha*exp(-beta*exp(-gamma*i)))
        if(fitness_diff>0):
            x_step=uniform(lb_x-x,ub_x-x)*gompertz[i]
            y_step=uniform(lb_y-y,ub_y-y)*gompertz[i]
                
            x+=x_step
            y+=y_step
          
        elif(pow(np.e,k)>randint(0,1)):
            x_step=uniform(lb_x-x,ub_x-x)*gompertz[i]
            y_step=uniform(lb_y-y,ub_y-y)*gompertz[i]
                
            x+=x_step
            y+=y_step
            
        if(fitness_current<threshold):
                x-=x_step
                y-=y_step
                for j in range(n-i):
                    switch=1
                    if(abs(fitness_diff)<0.1):
                        switch=2
                        print('fitness_diff = ',fitness_diff)
                        print('fitness for iteration ',i,' = ',fitness_current,' and Temp = ',T_n)
                        print('CONVERGED !!')
                        plt.plot(count,fitness)
                        threshold=fitness_current
                        exit()
                        
                    elif(fitness_diff>0):  
                        x_step=x_step/2
                        y_step=y_step/2
                        
                        fitness_current=ackley_func(x,y)
                        fitness[i+j]=fitness_current
                        fitness_new=ackley_func(x+x_step,y+y_step)
                        fitness_diff=fitness_current-fitness_new
                            
                        x+=x_step
                        y+=y_step
                        
                        print('fitness_diff = ',fitness_diff)
                        print('fitness for iteration ',i,' = ',fitness_current,' and Temp = ',T_n)
                        print(' At step = ',x_step,', ',y_step)

        print('fitness_diff = ',fitness_diff)
        print('fitness for iteration ',i,' = ',fitness_current,' and Temp = ',T_n)
        print(' At position = ',x,', ',y)
        print(' At step = ',x_step,', ',y_step)
        print()
                    
        
plt.plot(count,fitness)
plt.xlabel('x -->')
plt.ylabel('Gompertz Function -->')

print('NOT CONVERGED !!')

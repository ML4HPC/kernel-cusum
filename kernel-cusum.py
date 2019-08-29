import numpy as np
import matplotlib.pyplot as plt


# Function to generate a random sequence.
# The sequence may or may not have a change point,
# (mu1,sigma1) are the distribution parameters before the change,
# (mu2,sigma2) are the parameters after the change.
def random_sequence(length,
                    mu1 = 0, sigma1 = 1,
                    change_point=None,
                    mu2=1 ,sigma2 = 1): #ret[change_point] = 1st new pt 
    if change_point == None: change_point = length
    pre_change = sigma1*np.random.randn(change_point)  + mu1
    post_change = sigma2*np.random.randn(length-change_point) + mu2
    ret = np.concatenate((pre_change,post_change))
    return ret

# Demo of a random sequence
random_sequence(10,change_point=5)

def k(x,y,alpha=1,sigma2=2): #alpha=1: Laplace kernel. alpha=2: Gaussian
    return np.exp(-pow(np.linalg.norm(x-y),alpha)/sigma2)
#Notes on the kernel parameters: 
# Smaller sigma -> more noisy, rough signal.
# Bigger sigma -> sharper, smooth statistic
# Alpha smaller -> just seems to shrink the signal


def h(x1, y1, x2, y2): #See p. 729 in "a kernel two sample test"
    return k(x1, x2) + k(y1, y2) - k(x1, y2) - k(x2, y1)


## Two-sample test experiments

# Task: Determine whether two samples x_1,..,x_n and y_1,...,y_n are from the same distribution.

def linear_stat(x,y): #See p. 739 in "a kernel two sample test"
    m = len(x)/2
    z = 0.
    z_seq = []
    for i in range(0,m):
        z += h(x[2*i], y[2*i], x[2*i+1], y[2*i+1])
        z_seq.append(z/(i+1))
    return np.array(z_seq)


x = random_sequence(1000, mu1 = 5)

y = random_sequence(1000, mu1 = 10)

plt.plot(x)
plt.plot(y);

z = linear_stat(x,y); plt.plot(z);


y = random_sequence(1000, mu1 = 6)


plt.plot(x);plt.plot(y);


z = linear_stat(x,y); plt.plot(z);


## Change detection experiments

x = random_sequence(length = 1000, mu1 = 5, change_point=500, mu2= 7)

r = random_sequence(length = 1000, mu1 = 5) # Reference signal

plt.plot(x);

def cd_stat(x,r,delta = 0.1 ): #r for 'reference'
    m = len(x)/2
    z = 0.
    z_seq = []
    for i in range(0,m):
        z = max(0, z + h(x[2*i], r[2*i], x[2*i+1], r[2*i+1]) - delta)
        z_seq.append(z)
    return np.array(z_seq)

z = cd_stat(x,r); plt.plot(z);

x = random_sequence(length = 1000, mu1 = 5, change_point=500, mu2= 5, sigma2=0.1)


plt.plot(x);

z = cd_stat(x,r); plt.plot(z);

x = random_sequence(length = 1000, mu1 = 5, change_point=500, mu2= 5, sigma2=3)

plt.plot(x);

z = cd_stat(x,r); plt.plot(z);

r = random_sequence(length = 10000, mu1 = 5) # Reference signal

x = random_sequence(length = 10000, mu1 = 5) #no change

plt.plot(x);plt.plot(r);

z = cd_stat(x,r); plt.plot(z);plt.ylim(ymax=16)




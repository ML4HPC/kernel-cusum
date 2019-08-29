import numpy as np
import matplotlib.pyplot as plt

# Hyperparameters:

kalpha = 2 # Kernel exponent. kalpha=1: Laplace kernel. kalpha=2: Gaussian
ksigma = 1 # Kernel width
delta = 0.075 # Size of the 'don't care' region for changes.
thresh = 10 # Change is detected when statistic exceeds threshold.

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


# Kernel function.
def k(x,y): 
    return np.exp(-pow(np.linalg.norm(x-y),kalpha)/ksigma)

#Notes on the kernel parameters: 
# Smaller sigma -> more noisy, rough signal.
# Bigger sigma -> sharper, smooth statistic
# Alpha smaller -> just seems to shrink the signal


def h(x1, y1, x2, y2): #See p. 729 in "a kernel two sample test"
    return k(x1, x2) + k(y1, y2) - k(x1, y2) - k(x2, y1)


def cd_stat(x,r,delta = 0.075 ): #r for 'reference'
    m = len(x)/2
    z = 0.
    z_seq = []
    for i in range(0,m):
        z = max(0, z + h(x[2*i], r[2*i], x[2*i+1], r[2*i+1]) - delta)
        z_seq.append(z)
    return np.array(z_seq)

x = random_sequence(length = 1000, mu1 = 0, sigma1=1, change_point=500, mu2=0, sigma2=0.3)

r = random_sequence(length = 1000, mu1 = 0,sigma1=1) # Reference signal
plt.figure()
plt.plot(r,label='Reference');

plt.plot(x,label='Observations w/ change at 500');
plt.legend()

z = cd_stat(x,r);

plt.figure();
plt.plot(np.arange(0,1000,2),z,label='Statistic');
print( filter(lambda i: z[i] > thresh, range(len(z))))
detectedtime = 2*min( [len(z)] + filter(lambda i: z[i] > thresh, range(len(z))))


plt.title("KCUSUM Statistic")
plt.axvline(detectedtime,label='Detected time',color='orange')
plt.axvline(500,label='True time',color='brown')
plt.legend()
plt.show()





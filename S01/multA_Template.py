import numpy as np
import time
from matplotlib import pyplot as plt

def multA(a,x):
    """
    MULTA Return y = A*x for matrix A

    USAGE: y = multA(a,x);

    PARAMETERS:
    Input: a ... array of length n used to build A.
           x ... vector of length n.
    Output: y ... A*x
    """

    n = a.shape[0]
    # check for odd numbers:
    if ( n % 2 != 0 ):
        a[n/2] = a[n/2]/2.;

    # build matrix A
    A = np.diag(a,0);
    A = A + np.flipud(A);

    # perform matrix vector multiplication
    # complete here ...
    # y = ...
    y = np.zeros(n);
    for i in xrange(0, n):
        for j in xrange(0, n):
            y[i] += A[i,j] * x[j]

    return y

def multB(a,x):
 
    n = a.shape[0]
 
   # check for odd numbers:                                                                                                                                                                                       
    if ( n % 2 != 0 ):
        a[n/2] = a[n/2]/2.;

    # build matrix A
    B = np.diag(a,0);
    B = B + np.flipud(B);

    y = np.dot(B,x)

    return y


if __name__=='__main__':
    ## uncomment this for a simple test of your multA
    #
    a = np.array([1,2,3,4])
    x = np.array([1,1,1,1])
    print multA(a,x)

    # uncomment this for runtime comparisons
    nruns = 4	              # number of repetitions 
    ps = np.r_[5:13]              # problem sizes = 2^ps
    npoints = ps.size             # number of problems
    ts = 1e6*np.ones((npoints,2)) # container for running times

    # loop over different problem sizes
    cnt = 0                   
    for n in 2**ps:
        print 'Vector length: ', n, '\n'
        a = np.random.rand(n)
        x = np.random.rand(n)

        #loop for repetitions
        for k in np.r_[1:nruns]:
            t = time.time() # start time
            y0 = multA(a,x)
            t0 = time.time() - t # compute elapsed time
            ts[cnt,0]=np.minimum(ts[cnt,0],t0)

            t = time.time() # start time
            y1 = multB(a,x)
            t1 = time.time() - t # compute elapsed time
            ts[cnt,1]=np.minimum(ts[cnt,1],t1)

        cnt=cnt+1

    # create plot
    xval=2**ps

    # tricky way to place lines with slope 1 and slope 2                
    a1 = ts[np.ceil(npoints/2.),0]*np.sum(ts[:,0])/np.sum(xval)
    a2 = np.sum(ts[:,0])/np.sum(xval**2)

    b1 = ts[np.ceil(npoints/2.),1]*np.sum(ts[:,1])/np.sum(xval)
    b2 = np.sum(ts[:,1])/np.sum(xval**2)

    plt.figure()
    plt.loglog(xval,ts[:,0]   ,'r+-' \
              ,xval,a1*xval   ,'--' \
              ,xval,a2*xval**2,'--' \
              ,xval,ts[:,1]   ,'b+-' \
              ,xval,b1*xval   ,'--' \
              ,xval,b2*xval**2,'--')
    plt.legend(('multA','O(n)' \
               ,'O(n^2)','multB','O(n)','O(n^2)') \
              ,loc='upper left')
    plt.xlabel('{ problem size n}')
    plt.ylabel('{ average runtime (s)}')
    plt.savefig('a1_template.eps')
    plt.show()


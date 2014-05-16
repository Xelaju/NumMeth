from numpy import *
import time


def conv_fft(x, y):

    x_h = fft.fft(x)
    #print 'x_h: ', x_h
    #print 'size: ', size(x_h)
    y_h = fft.fft(y)
    z_h = x_h * y_h
    z = fft.ifft(z_h)
    return z


def poly_mult(a, b):
    n = a.size
    N = 2 * n - 1
    z = zeros(N + 1)
    for i in xrange(N + 1):
        for j in xrange(max(0, i - n + 1), min(i + 1, n)):
            z[i] += a[i - j] * b[j]
    return z[0:N]


def poly_multfft(a, b):

    n_a = size(a)
    n_b = size(b)
    n_c = 2 * max(n_a,n_b) - 1
    a = append(a,zeros(n_c - n_a))
    b = append(b,zeros(n_c - n_b))
    c = conv_fft(a,b)

    return c.real 


n = 1.e4

# a = ones(n)
a = random.random(n)

t0 = time.clock()
z1 = poly_mult(a, a)
print(time.clock() - t0)
t0 = time.clock()
z2 = poly_multfft(a, a)
print(time.clock() - t0)
print 'Welche der beiden Implementierungen ist effizienter? Die Rechenzeit spricht Baende!'

print(z1)
print(z2)

print linalg.norm(z2 - z1)

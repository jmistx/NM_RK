import matplotlib

import itertools

import array

import math

from pprint import pprint

def v_summ(v1, v2):
	return tuple((v1[i] + v2[i] for i in xrange(len(v1))))

def v_mult(a, v1):
	return tuple((a*v1[i] for i in xrange(len(v1))))	

def v_func(f, x, v):
	return tuple((f[i](x, v) for i in xrange(len(f))))

def v_abs(v):
	return tuple((abs(v[i]) for i in xrange(len(v))))	

def v_diff(v1, v2):
	return v_abs(v_summ(v1, v_mult(-1.0, v2)))

def rk2(func, xn, vn, hn): 
	return v_summ(vn, v_mult(hn/2.0, v_summ(v_func(func, xn, vn), v_func(func, xn + hn, v_summ(vn, v_mult(hn, v_func(func, xn, vn)))))))



x0 = 0.0
w0 = 0.0
phi0 = 1.0
phi_0 = 5.0
h = 0.01
k = 1.0
m = 1.0
n = 1.0
F = 100.0
J = 3.0
g = 1.0
b = 1.0

def _phi(v):
	return v[0]

def _phi1(v):
	return v[1]

def _w(v):
	return v[2]


#f0 = lambda x, v: _phi1(v)
#f1 = lambda x, v: n*n*_w(v)*math.sin(_phi(v))*math.cos(_phi(v)) - g*math.sin(_phi(v)) - (b/m)*+_phi1(v)
#f2 = lambda x, v: (k*math.cos(_phi(v)) - F) / J

f0 = lambda x, v: _phi(v)
f1 = lambda x, v: -2.0*_phi1(v)
f2 = lambda x, v: 1.0/math.sin(_w(v))

s0 = lambda x, v: math.exp(x)
s1 = lambda x, v: math.exp(-2*x) * 2.0
s2 = lambda x, v: math.acos(-x + math.cos(0.7))

func = (f0, f1, f2) 
sol = (s0, s1, s2)
#vn = (phi0, phi_0, w0)
vn = (1.0, 2.0, 0.7)

hn = 0.001
xn = 0

history = []


method = rk2

for i in xrange(100):
	vn_1 = method(func, xn, vn, hn)

	vn_half = method(func, xn, vn, hn/2.0)
	vn_wave = method(func, xn + hn/2.0, vn_half, hn/2.0)
	vn_1_vn_wave = v_summ(vn_wave, v_mult(-1.0, vn_1))

	p = 2
	S = v_mult(1.0 / (2.0**p - 1), vn_1_vn_wave)

	vn = vn_1
	xn += hn
	un = v_func(sol, xn, 0)
	un_vn = v_diff(vn, un)

	history.append({"n": i, 
					"x": xn, 
					"v": vn, 
					"h": hn, 
					"u": un, 
					"uv": un_vn, 
					"v2": vn_1_vn_wave,
					"S": S})

	# print xn, "\t", vn, 

for h in history:
	print h['n'], "\t", h['x'], '\t'
	print "      v:", h['v']
	print "      u:", h['u']
	print "      S:", h['S']
	print "     v2:", h['v2']
	print "|u - v|:", h['uv']
	
maxUv = tuple( max((h['uv'][i] for h in history)) for i in xrange(len(vn)) )
print 'max |u - v|:', maxUv
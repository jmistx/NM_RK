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
Nmax = 100
eps = 10.0**-5

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
f2 = lambda x, v: 3.0*_w(v)

s0 = lambda x, v: math.exp(x)
s1 = lambda x, v: math.exp(-2*x) * 2.0
s2 = lambda x, v: math.exp(3.0*x) * 0.7

func = (f0, f1, f2) 
sol = (s0, s1, s2)
#vn = (phi0, phi_0, w0)
vn = (1.0, 2.0, 0.7)

hn = 0.1
xn = 0
C1 = 0
C2 = 0

history = []

method = rk2
p = 2

history.append({"n": 0, 
				"x": xn, 
				"v": vn, 
				"v2": vn,
				"h": hn, 
				"u": vn, 
				"uv": v_diff(vn, vn), 
				"v1v2": v_diff(vn, vn),
				"S": v_diff(vn, vn),
				"Smax": 0,
				"C1": C1,
				"C2": C2})

for i in xrange(1, Nmax + 1):
	vn_1 = method(func, xn, vn, hn)

	vn_half = method(func, xn, vn, hn/2.0)
	vn_wave = method(func, xn + hn/2.0, vn_half, hn/2.0)

	vn_1_vn_wave = v_summ(vn_wave, v_mult(-1.0, vn_1))

	S = v_mult(1.0 / (2.0**p - 1), vn_1_vn_wave)
	S_max = max(v_abs(S))

	h_prev = hn
	xn_current = xn + hn
	#local error control
	if eps / (2.0**(p+1)) <= S_max <= eps:
		xn += hn
		vn = vn_1

	elif eps < S_max:
		hn = hn / 2.0
		C1 += 1

	else:
		xn += hn
		vn = vn_1
		hn = hn * 2.0
		C2 += 1

	#compare exact solution
	un = v_func(sol, xn, 0)
	un_vn = v_diff(vn, un)

	history.append({"n": i, 
					"x": xn_current, 
					"v": vn_1, 
					"v2": vn_wave,
					"h": h_prev, 
					"u": un, 
					"uv": un_vn, 
					"v1v2": vn_1_vn_wave,
					"S": S,
					"Smax": S_max,
					"C1": C1,
					"C2": C2})

	# print xn, "\t", vn, 

for h in history:
	print 'i:', h['n'], "\t", 'x:', h['x'], '\t','h:',h['h'], '\t','C1:',h['C1'], '\t','C2:', h['C2']
	print "       v:", h['v']
	print "      v2:", h['v2']
	print "       u:", h['u']
	print "|v - v2|:", h['v1v2']
	print "       S:", h['S']
	print "    Smax:", h['Smax']
	print " |u - v|:", h['uv']
	
maxUv = tuple( max((h['uv'][i] for h in history)) for i in xrange(len(vn)) )
print 'max |u - v|:', maxUv
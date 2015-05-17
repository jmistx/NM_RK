# -*- coding: utf-8 -*-

import matplotlib

import itertools

import array

import math

import yaml

import sys

import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

font = {'family': 'Verdana',
        'weight': 'normal'}
rc('font', **font)

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

def rkMenson(func, xn, vn, hn):
	k1 = v_func(func, xn, vn)
	k2 = v_func(func, xn + hn/3.0, v_summ(vn, v_mult(hn/3.0, k1)))
	k3 = v_func(func, xn + hn/3.0, v_summ(vn, v_summ(v_mult(hn/6.0, k1), v_mult(hn/6.0, k2))))
	k4 = v_func(func, xn + hn/2.0, v_summ(vn, v_summ(v_mult(hn/8.0, k1), v_mult(3.0*hn/8.0, k3))))
	k5 = v_func(func, xn + hn, v_summ(vn, v_summ(v_mult(hn/2.0, k1), v_summ(v_mult(-3.0*hn/2.0, k3), v_mult(2.0*hn, k4)))))
	vn1 = v_summ(vn, v_mult(hn/10.0, v_summ(k1, v_summ(v_mult(3.0, k3), v_summ(v_mult(4.0, k4), v_mult(2.0, k5))))))
	S = v_mult(hn / 30.0, v_summ(v_mult(2.0, k1), v_summ(v_mult(-9.0, k3), v_summ( v_mult(8.0, k4), v_mult(-1.0, k5) ))))
	return vn1, S

def show_S_plot(result_s_max, result_h):
	f = plt.figure(1)
	iters = range(len(result_s_max))

	plt.subplot(2, 1, 1)
	plt.title(u'Оценка локальной погрешности')
	plt.xlabel('i')
	plt.ylabel('Smax')
	plt.plot(iters, result_s_max)

	plt.subplot(2, 1, 2)
	plt.title(u'Размер шага')
	plt.xlabel('i')
	plt.ylabel('h')
	plt.plot(iters, result_h)
	f.show()

def show_Func_plot(result, captions):
	g = plt.figure(2)
	_, v = result[0]
	plotCount = len(v)

	x_result = [x for x, v in result]
	
	for i in xrange(plotCount):
		plt.subplot(plotCount, 1, i + 1)

		if i == 0:
			plt.title(u'Значения функций')

		plt.ylabel(captions[i])

		if i == plotCount - 1:
			plt.xlabel(u'x, время')

		
		plt.plot(x_result, [v[i] for x, v in result])

	g.show()

def show_Fase_portrait(result, captions, members):
	g = plt.figure(3)
	x_index, y_index = members

	x_caption = captions[x_index]
	y_caption = captions[y_index]

	x_axis = [v[x_index] for x, v in result]
	y_axis = [v[y_index] for x, v in result]

	plt.xlabel(x_caption)
	plt.ylabel(y_caption)

	plt.title(u'Фазовый портрет: ' + x_caption + ' ' + y_caption)
	plt.plot(x_axis, y_axis)

	g.show()

def show_exact_and_test_solution(result, captions, solution):
	g = plt.figure(4)
	
	funcCount = len(solution)
	plotCount = funcCount * 2;
	x_result = [x for x, v in result]

	for i in xrange(funcCount):
		plt.subplot(plotCount, 1, i * 2 + 1)

		if i == 0:
			plt.title(u'Сравнение точного и численного решения: ')

		exact_solution = [solution[i](x, 0) for x in x_result]
		numeric_solution = [v[i] for x, v in result]

		plt.ylabel(captions[i])
		plt.plot(x_result, exact_solution, color="g")
		plt.plot(x_result, numeric_solution, color="r")

		plt.subplot(plotCount, 1, i * 2 + 2)

		diff_solutions = [ abs(u - v) for u, v in itertools.izip(exact_solution, numeric_solution)]
		plt.plot(x_result, diff_solutions, color="b")
		plt.ylabel('|u - v|')
		if i == funcCount - 1:
			plt.xlabel(u'x, время: ')

	g.show()

if getattr(sys, 'frozen', False):
    application_path = os.path.dirname(os.path.realpath(sys.executable))
elif __file__:
    application_path = os.path.dirname(os.path.realpath(__file__))

cfg = yaml.load(open(os.path.join(application_path, 'params.yaml')))

#params
right_edge_stop_enabled = cfg['right_edge_stop_enabled']
table_output_enabled = cfg['table_output_enabled']
test_task = cfg['test_task']
method = rk2 if cfg['method'] == 'rk2' else rkMenson
p = cfg['p']

Xend = cfg['Xend']
XendEps = cfg['XendEps']

Nmax = cfg['Nmax']
eps = cfg['eps']

#task params
x0 = cfg['x0']
w0 = cfg['w0']
phi0 = cfg['phi0']
phi_0 = cfg['phi_0']
h = cfg['h']
k = cfg['k']
m = cfg['m']
n = cfg['n']
F = cfg['F']
J = cfg['J']
g = cfg['g']
b = cfg['b']

def _phi(v):
	return v[0]

def _phi1(v):
	return v[1]

def _w(v):
	return v[2]

if test_task:
	captions = ['phi', 'phi\'', 'w']

	f0 = lambda x, v: _phi(v)
	f1 = lambda x, v: -2.0*_phi1(v)
	f2 = lambda x, v: 3.0*_w(v)

	s0 = lambda x, v: math.exp(x)
	s1 = lambda x, v: math.exp(-2*x) * 2.0
	s2 = lambda x, v: math.exp(3.0*x) * 0.7

	vn = (1.0, 2.0, 0.7)

	sol = (s0, s1, s2)

else:
	captions = ['phi', 'phi\'', 'w']

	f0 = lambda x, v: _phi1(v)
	f1 = lambda x, v: n*n*_w(v)*_w(v)*math.sin(_phi(v))*math.cos(_phi(v)) - g*math.sin(_phi(v)) - (b/m)*+_phi1(v)
	f2 = lambda x, v: (k*math.cos(_phi(v)) - F) / J
	vn = (phi0, phi_0, w0)

func = (f0, f1, f2) 

hn = 0.1
xn = 0
C1 = 0
C2 = 0

history = []
result_s_max = [0]
result_h = [hn]

result = [(xn, vn)]

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

step = 0

if method == rk2:
	while(step < Nmax):
		step += 1
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
			result.append((xn, vn))
			result_s_max.append(S_max)
			result_h.append(h_prev)

		elif eps < S_max:
			hn = hn / 2.0
			C1 += 1

		else:
			xn += hn
			vn = vn_1
			hn = hn * 2.0
			C2 += 1
			result.append((xn, vn))
			result_s_max.append(S_max)
			result_h.append(h_prev)

		#compare exact solution
		if test_task:
			un = v_func(sol, xn, 0)
			un_vn = v_diff(vn, un)
		else:
			un = 0
			un_vn = 0

		history.append({"n": step, 
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

		if right_edge_stop_enabled:
			if Xend - xn < XendEps:
				break
if method == rkMenson:
	while(step < Nmax):
		step += 1
		vn_1, S = method(func, xn, vn, hn)

		S_max = max(v_abs(S))

		h_prev = hn
		xn_current = xn + hn
		#local error control
		if eps / (2.0**(p+1)) <= S_max <= eps:
			xn += hn
			vn = vn_1
			result.append((xn, vn))
			result_s_max.append(S_max)
			result_h.append(h_prev)

		elif eps < S_max:
			hn = hn / 2.0
			C1 += 1

		else:
			xn += hn
			vn = vn_1
			hn = hn * 2.0
			C2 += 1
			result.append((xn, vn))
			result_s_max.append(S_max)
			result_h.append(h_prev)

		#compare exact solution
		if test_task:
			un = v_func(sol, xn, 0)
			un_vn = v_diff(vn, un)
		else:
			un = 0
			un_vn = 0

		history.append({"n": step, 
						"x": xn_current, 
						"v": vn_1, 
						"v2": 0,
						"h": h_prev, 
						"u": un, 
						"uv": un_vn, 
						"v1v2": 0,
						"S": S,
						"Smax": S_max,
						"C1": C1,
						"C2": C2})

		if right_edge_stop_enabled:
			if Xend - xn < XendEps:
				break

if table_output_enabled:
	for h in history:
		print 'i:', h['n'], "\t", 'x:', h['x'], '\t','h:',h['h'], '\t','C1:',h['C1'], '\t','C2:', h['C2']
		print "       v:", h['v']
		print "      v2:", h['v2']
		print "       u:", h['u']
		print "|v - v2|:", h['v1v2']
		print "       S:", h['S']
		print "    Smax:", h['Smax']
		print " |u - v|:", h['uv']
	

if test_task:
	maxDiff = [0, 0, 0]
	for i in xrange(len(maxDiff)):
		x_result = [x for x, v in result]
		exact_solution = [sol[i](x, 0) for x in x_result]
		numeric_solution = [v[i] for x, v in result]
		maxDiff[i] = max((abs(u - v) for u, v in itertools.izip(exact_solution, numeric_solution)))
	print 'max |u - v|:', maxDiff

hmin = min((h['h'] for h in history))
hmax = min((h['h'] for h in history))

print 'C1:', C1
print 'C2:', C2
print 'Steps', step
print 'h min', hmin
print 'h max', hmax
print 'xn', xn



show_S_plot(result_s_max, result_h)
show_Func_plot(result, captions)
show_Fase_portrait(result, captions, (1, 0))

if test_task:
	show_exact_and_test_solution(result, captions, sol)

raw_input()
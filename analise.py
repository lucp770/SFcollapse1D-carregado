import os
import numpy as np
from matplotlib import pyplot as plt
import utilities

out_dir = './SFcollapse2/out/'

# initial amplitude and user parameters.
phi0 = 0.3;
Nx0 = 320 + 1;
c_final = 0.0001
w_q = 1.0/3.0
# precisamos -1<w_q <0 quando c>0.

# wq < -1(como -4/3) e c<0. TENTAR ISSO.

anisotropic_exponent = 2.0 + 3.0 * w_q

rmax = 16.0
sinhA             = rmax;
sinhW             = 0.2;
inv_sinhW         = 1.0/sinhW;
sinh_inv_W        = np.sinh(inv_sinhW);
A_over_sinh_inv_W = sinhA / sinh_inv_W;

x0_min = 0.0
x0_max = 1.0
dx0 = (x0_max - x0_min)/(Nx0-1.0);
Ngz0 = 0
inv_dx0 = 1.0/dx0


def define_x_array():

	x = np.zeros(Nx0)

	for j in range(0,len(x)):
		x[j] = x0_min + (j-Ngz0)*dx0
	return x

x = define_x_array()
# print(x[1]-x[0], x[10] - x[9]);

# calculate the correspondent r
r_ito_x0 = np.zeros_like(x)

for j in range(len(x)):
	r_ito_x0[j] = A_over_sinh_inv_W * np.sinh( x[j] * inv_sinhW );


# set up initial condition

R0 = 0
delta = 1
r = np.linspace(0, rmax, 1000)

factor  = (r_ito_x0-R0)/delta**2

expfactor = (r_ito_x0-R0)*factor;
initial_condition = phi0*np.exp(-expfactor);

def show_initial_condition():
	plt.plot(r,initial_condition, '-k')
	plt.show()

Pi = np.zeros_like(x);

# initial condition for Phi
Phi = np.zeros_like(x)

Phi = -2.0*factor*phi0*np.exp(-expfactor);
Phi[0] = 0.0 #reforcing

def show_Phi():
	plt.plot(r,Phi, '-k')
	plt.show()

# calculate a and alpha
a = np.zeros_like(x);
alpha = np.zeros_like(x);

a[0] = 1.0
alpha[0] = 1.0;
# populate the rest of the array (a-> Hamiltonian constraint; alpha -> polar slicing constraint)

def Hamiltonian_constraint(j,dbg = False):
	# Hamiltonian constraint for a iteration j, initiate in j=1
	A = np.log(a[j-1]);
	avgPhi  = .5*(Phi[j] + Phi[j-1]);
	avgPi = .5*(Pi[j] +Pi[j-1]);
	PhiSqr = avgPhi**2;
	PiSqr = avgPi**2;
	midx0 = 0.5 * ( x[j] + x[j-1] );
	
	# 
	PhiPiTerm = 2.0 * np.pi * ((sinhA**2) * inv_sinhW) * np.sinh(midx0*inv_sinhW) * np.cosh(midx0*inv_sinhW) /(sinh_inv_W**2) * ( PhiSqr + PiSqr );
	# print('PhiPiTerm',PhiPiTerm)
	half_invr = 0.5/( sinhW * np.tanh(midx0*inv_sinhW));
	# print('half_invr',half_invr)
	r_sinh = A_over_sinh_inv_W * np.sinh(inv_sinhW * midx0);

	ans_fluid_term = (3 * c_final * w_q)/(2*r_sinh**anisotropic_exponent);
	A_old = np.log(a[j-1]);

	tmp0 = half_invr * np.exp(A_old+A)
	# print('tmp0',tmp0);

	superior_limit = 100#before, when was 1000 i got overflow.
	inferior_limit = 0

	funcao = lambda x: inv_dx0 * (x - A) + tmp0 - half_invr - PhiPiTerm + ans_fluid_term * np.exp(x+A)
	# print('ans_fluid_term',ans_fluid_term);

	if dbg == True:
		c_final2 = np.array([-0.0001,-0.00001,-0.000001, -0.0000001])
		c_final3 = np.array([ -0.00000000001,-0.0000000000001])
		c_final4 = np.array([0.0001,0.00001])

		c_final5 = np.array([0.000001, 0.0000001, 0.00000000001,0.0000000000001])

		c_final6 = np.array([-1e-7])

		c_zero = np.array([0.0])
		# set the space of shots for the root.
		x_root = np.linspace(-.5,5,1000)
		ans_fluid_term = (3 * c_final6 * w_q)/(2*r_sinh**anisotropic_exponent);

		zero_axis = [0 for i in x_root]

		print('deveria ser um array: ans_fluid_term \n', ans_fluid_term)
		plt.figure()

		for i in range(len(ans_fluid_term)):
			y = inv_dx0 * (x_root - A) + tmp0 - half_invr - PhiPiTerm + ans_fluid_term[i] * np.exp(x_root+A)
			title = 'c = {}'.format(c_final6[i])

			plt.subplot(1,2,i+1)
			plt.plot(x_root,y,label = title)
			plt.plot(x_root, zero_axis, '-k')
			plt.legend()
		plt.show()

	a_j = utilities.bissection_method(inferior_limit, superior_limit, funcao,tolerancia = 1e-6)
	a_j = np.exp(a_j)
	return a_j


def Apply_Hamiltonian_constraint():
	# invoque the hamiltonian constraint in all the points.

	for j in range(1,len(a)):
		a[j] = Hamiltonian_constraint(j)

Apply_Hamiltonian_constraint()
# print(a)


alpha = np.zeros_like(x)
alpha[0] = 1.0

def polar_slicing(j):
	b = a[j] + a[j-1];
	c = a[j] - a[j-1];
	midway_r = sinhW *np.tanh( inv_sinhW * 0.5 * (x[j] + x[j-1]) );
	r_sinh = A_over_sinh_inv_W * np.sinh(inv_sinhW * (0.5 * (x[j] + x[j-1])));
	ans_fluid_term = (3 * c_final * w_q)/(2* (r_sinh**anisotropic_exponent));
	d = (1.0 - 0.25 * b**2)/( 2.0 * midway_r ) - inv_dx0 * c / b - 0.25 * ans_fluid_term *b**2;
	result = alpha[j-1]*( 1.0 - d*dx0 )/( 1.0 + d*dx0 );
	return result

def Apply_Polar_Slicing():
	for j in range(1, len(alpha)):
		alpha[j] = polar_slicing(j)
		
# Apply_Polar_Slicing()


def test_stability():
	functions = Hamiltonian_constraint(1,A=np.log(a[j-1]), dbg =True)
	print(len(functions))
	x = np.linspace(-2,2,1000)
	zero_line = [0 for i in x]

	print(functions)

	labels = [-0.0001,-0.00001,-0.000001, -0.0000001, -0.00000000001,-0.0000000000001,0.0001,0.00001,0.000001, 0.0000001, 0.00000000001,0.0000000000001]

	for i in range(len(functions)):
		y = functions[i](x)
		# print(functions[i](x)[:10])
		label = labels[i]

		plt.plot(x,y,label = label)

	plt.plot(x,zero_line)
	plt.legend()
	plt.show()


	# y=funcao(x)
	# plt.plot(x,y)
	# plt.plot(x,zero_line,'-r')
	
	# plt.scatter([0],[0])
	# plt.xlabel('chutes para a raiz')
	# plt.ylabel('valor da funcao')
	# plt.show()

	# print('valor em x=0: ', funcao(0))

	# print('funcao',y[:5])
# test_stability()

# plot the hamiltonian constraint for 3 steps.

def comparison_plot(filepath):
	y_file =[];
	x_file =[];
	r_file = [];

	with open(filepath, 'r', encoding = 'utf-8') as text:
			for line in text:
				l = line.split()
				y_file.append(float(l[2]))
				r_file.append(float(l[1]))
				x_file.append(float(l[0]))
	return [r_file, x_file, y_file]

def show_comparison_a_alpha():

	a_dir= out_dir+'a_00000000.dat'
	alpha_dir  = out_dir+'alpha_00000000.dat'

	r_file = comparison_plot(a_dir)[0]
	a_file = comparison_plot(a_dir)[2]
	alpha_file = comparison_plot(alpha_dir)[2]

	plt.figure()

	plt.subplot(2,2,1)
	plt.plot(r_file, a_file, '-r', label = 'a(0,t) - SFcollapse1d');
	plt.legend()

	plt.subplot(2,2,2)
	plt.plot(r_ito_x0,a,label = 'a(0,t)')
	plt.legend()

	plt.subplot(2,2,3)
	plt.plot(r_file, alpha_file, '-r', label = 'alpha(0,t) - SFcollapse1d');
	plt.legend()

	plt.subplot(2,2,4)
	plt.plot(r_ito_x0,alpha,label = 'alpha(0,t)')
	plt.legend()

	plt.show()

# show_comparison_a_alpha()
# 

def show_comparison_Phi_Pi():

	Phi_dir= out_dir+'Phi_00000000.dat'
	Pi_dir  = out_dir+'Pi_00000000.dat'

	r_file = comparison_plot(Phi_dir)[0]
	Phi_file = comparison_plot(Phi_dir)[2]
	Pi_file = comparison_plot(Pi_dir)[2]

	print(len(Phi_file), len(r_file), len(Pi), len(r_ito_x0))
	plt.figure()

	plt.subplot(2,2,1)
	plt.plot(r_file, Phi_file, '-r', label = 'Phi(0,t) - SFcollapse1d');
	plt.legend()

	plt.subplot(2,2,2)
	plt.plot(r_ito_x0,Phi,label = 'Phi(0,t)')
	plt.legend()

	plt.subplot(2,2,3)
	plt.plot(r_file, Pi_file, '-r', label = 'Pi(0,t) - SFcollapse1d');
	plt.legend()

	plt.subplot(2,2,4)
	plt.plot(r_ito_x0,Pi,label = 'Pi(0,t)')
	plt.legend()

	plt.show()
# show_comparison_Phi_Pi()



Hamiltonian_constraint(1,dbg = True)
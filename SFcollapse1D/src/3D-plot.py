"""
generate 3-d plots for the scalarfield, lapse function and a (radial component of the metric)

"""

from matplotlib import pyplot as plt
import numpy as np
import os


#1st step: print the names of files about a particular function
a_files,alpha_files,scalar_files,Phi_files,Pi_files =[],[],[],[],[]
lista = os.listdir('../out/')

print("coletando nome dos arquivos ....")
for i in lista:
	search_element = i[:2]
	if search_element == 'a_':
		a_files.append(i)
	elif search_element == 'al':
		alpha_files.append(i)
	elif search_element == 'Pi':
		Pi_files.append(i)
	elif search_element == 'Ph':
		Phi_files.append(i)
	elif search_element == 'sc':
		scalar_files.append(i)

print("\n nomes coletados \n")
def read_append_file(filename):
	pivot = []
	with open(filename, 'r', encoding = 'utf-8') as text:
		for line in text:
			line = line.split()
			pivot.append(float(line[2]))
		pivot = pivot[:int(6*len(pivot)/7)]#only take half of the inputs
	return pivot

#create arrays
a=[]
print("\n coletando valores de a .... \n")
for i in a_files:
	a.append(read_append_file('../out/'+i))
a = np.array(a)

# print("\n coletando valores de alpha .... \n")
# alpha=[]
# for i in alpha_files:
# 	alpha.append(read_append_file('../out/'+i))
# alpha = np.array(alpha)

# print("\n coletando valores de Pi .... \n")
# Pi=[]
# for i in Pi_files:
# 	Pi.append(read_append_file('../out/'+i))
# Pi = np.array(Pi)

# print("\n coletando valores de Phi .... \n")
# Phi=[]
# for i in Phi_files:
# 	Phi.append(read_append_file('../out/'+i))
# Phi = np.array(Phi)

# print("\n coletando valores do campo escalar .... \n")
# scalar=[]
# for i in scalar_files:
# 	scalar.append(read_append_file('../out/'+i))
# scalar = np.array(scalar)

print("\n coletando valores de r e t .... \n")
#create r and t arrays.


#check if the values are NaN

nan_values =np.isnan(a)
formato = nan_values.shape
for i in range(formato[0]):
	for j in range(formato[1]):
		if nan_values[i][j] == True: a[i][j] =2.5




print(a[3110])
r=[]
with open('../out/'+a_files[0], 'r', encoding = 'utf-8') as text:
	for line in text:
		line = line.split()
		r.append(float(line[1]))
r = r[:int(6*len(r)/7)]
r=np.array(r)

t = [i for (i,j) in enumerate(a_files)]
t = np.array(t)

#create 3-D plots:
rr,tt = np.meshgrid(r,t)


print("t:",t.shape,'\n',"r:",r.shape,'\n',"tt:",tt.shape,'\n','a: ', a.shape,'\n', )
print("a_files",len(a_files))
fig =plt.figure()
ax = plt.axes(projection = '3d')

ax.plot_surface(rr,tt,a,cmap ='viridis')
ax.set_xlabel('r')
ax.set_ylabel('t')
ax.set_zlim(0,3)
plt.show()

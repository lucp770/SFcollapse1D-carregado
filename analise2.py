import os
import numpy as np
from matplotlib import pyplot as plt

def generate_bissection_plot(filepath):
	y_file =[];
	x_file =[];

	with open(filepath, 'r', encoding = 'utf-8') as text:
			for line in text:
				l = line.split()
				y_file.append(float(l[1]))
				x_file.append(float(l[0]))
	return [x_file[:300], y_file[:300]]


filepath = './SFcollapse1D/bissection_output.dat'
# filepath = './SFcollapse1D/out/a_00005900.dat'
# filepath = './SFcollapse1D/out/a_00000173.dat'

array = generate_bissection_plot(filepath)

x_array = np.linspace(-30, 0,200)
y_array = 0 * x_array

plt.plot(array[0], array[1])
plt.xlabel('r')
plt.ylabel('$a (r, t_{max}) $')
plt.plot(x_array,y_array, '-r')
plt.show()
print(array[1][:10])

import numpy as np
from matplotlib import pyplot as plt

########### analisar os passos inicais do método de newton ##################


###########definir o método de newton
#definir uma funcao que pode ser aplicada em cada passo do programa.


def Newton(A_old, iter_max):




a0 = 1# a no ponto a(t,0) precisa satisfazer "elementary flatness"
N = 320#input de SFCollapse


a=np.linspace(a0,N,1000)

A_old = np.log(a[i-1])



for i in range(1,N):
	

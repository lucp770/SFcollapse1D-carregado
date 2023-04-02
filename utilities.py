import numpy as np


def bissection_method(inferior, superior, funcao,tolerancia):
	
	max_iteration = 10000

	f_inferior = funcao(inferior)
	f_superior = funcao(superior)

	error = (f_superior - f_inferior)
	ponto_medio =0

	iteration = 0

	while error > tolerancia and iteration < max_iteration:

		if f_inferior <= 0 and f_superior >0:

			ponto_medio = (inferior+superior)/2
			f_c = funcao(ponto_medio)

			if f_c > 0:
				superior = ponto_medio
				f_superior = f_c
			elif f_c < 0:
				inferior = ponto_medio
				f_inferior = f_c
			elif f_c ==0:
				return ponto_medio

		elif f_superior >= 0 and f_inferior <0:

			ponto_medio = (inferior +superior)/2
			f_c = funcao(ponto_medio);

			if f_c >0:
				superior = ponto_medio 
				f_superior = f_c

			elif f_c <0:
				inferior = ponto_medio
				f_inferior = f_c

			elif f_c == 0:
				return f_c
		else:
			print('nao é possivel resolver')
			x = np.linspace(inferior,superior,1000)
			y=funcao(x)
			plt.title('funcao que nao foi encontrada raiz: ')
			plt.plot(x,y, '-r')
			break

		iteration +=1 

	return ponto_medio

def Newton_method():
	pass


if __name__ == "__main__":
	funcao = lambda x: x-np.cos(x)
	raiz=bissection_method(0, 10, funcao,tolerancia = 1e-6)
	print(raiz)
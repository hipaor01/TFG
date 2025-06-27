from sage.all import *
from random import choice
import random
import argparse

#Función para convertir los elementos pasados por línea de comandos en un vector de elementos de un cuerpo finito, para el caso no primo
def procesarVectores(vector, a, subcuerpoPrimo):
	vectorField = []

	for elemento in vector:
		coeficientesBase = [subcuerpoPrimo(int(coefi)) for coefi in elemento.split(',')] #Los coeficientes deben ser elementos del subcuerpo primo
		elementoField = sum(coefi * a**i for i,coefi in enumerate(reversed(coeficientesBase)))
		vectorField.append(elementoField)

	return vectorField 


class GRS:
	#El parámetro alphaCuerpo/vCuerpo lo utilizamos para indicar si el vector alpha/v lo hemos pasado ya procesado como elementos del cuerpo correspondiente o no
	def __init__(self,n,p,polinomioBase,ordenSubcuerpoPrimo,k,alpha,v,alphaCuerpo=False,vCuerpo=False):
		self.p = p #Entero representando el orden del cuerpo finito sobre el cual construir el código
		#Cuerpo finito sobre el cual construir el código
		self.a = None
		self.subcuerpoPrimo = None
		if not is_prime(self.p):
			if not is_prime(ordenSubcuerpoPrimo):
				print("El orden del subcuerpo primo debe ser primo")
				sys.exit(1)

			self.subcuerpoPrimo = GF(ordenSubcuerpoPrimo)
			try:
				R = PolynomialRing(self.subcuerpoPrimo,'x')
				polinomioModulo = R(polinomioBase)
				self.field = GF(p, name='a', modulus=polinomioModulo)

				assert ordenSubcuerpoPrimo**(polinomioModulo.degree()) == self.p
			except:
				print("El polinomio y subcuerpo primo no pueden construir el cuerpo deseado")
				sys.exit(1)
			
			self.a = self.field.gen() #Raíz del polinomio generador
		else:	
			self.field = GF(p) 
		self.n = n #Longitud del código
		try:
			assert self.n <= self.p and self.n >= 0
		except AssertionError as e:
			print("La longitud del código debe ser menor o igual que el orden del cuerpo y no negativa")
			sys.exit(1)	
		self.k = k #Dimensión del código
		try:
			assert self.k <= self.n and self.k >= 0
		except AssertionError as e:
			print("La dimensión debe ser menor o igual que la longitud del código y no negativa")
			sys.exit(1)	

		#En el caso de que el cuerpo sea primo, como los elementos de alpha son enteros podemos procesarlos fácilmente
		alphaField = None
		try:
			if is_prime(self.p) or alphaCuerpo:
				alphaField = [self.field(alpha[i]) for i in range(len(alpha))]	
			#En caso contrario, como los elementos del cuerpo tenemos que expresarlos en términos de una raíz del polinomio generador, pasamos cada elemento de alpha
			#como un vector separado por comas, que representa los coeficientes del elemento en la base dada por dicha raíz	
			else:
				alphaField = procesarVectores(alpha, self.a, self.subcuerpoPrimo)
		except:
			print("No ha sido posible procesar el vector alpha")
			sys.exit(1)	
		self.alpha = vector(self.field,alphaField) #Vector alpha de la definición de Código GRS

		try:
			assert len(self.alpha) == self.n
		except AssertionError as e:
			print("La longitud de alpha debe coincidir con la longitud del código")
			sys.exit(1)

		try:
			assert len(self.alpha) == len(set(self.alpha))
		except AssertionError as e:
			print("El vector alpha no puede tener elementos repetidos")
			sys.exit(1)
		
		#Para v hacemos lo mismo que para alpha
		vField = None
		try:
			if is_prime(self.p) or vCuerpo:
				vField = [self.field(v[i]) for i in range(len(v))]
			else:
				vField = procesarVectores(v,self.a, self.subcuerpoPrimo)
		except:
			print("No ha sido posible procesar el vector v")
			sys.exit(1)		
		self.v = vector(self.field,vField) #Vector v de la definición de Código GRS

		try:
			assert len(self.v) == self.n
		except AssertionError as e:
			print("La longitud de v debe coincidir con la longitud del código")
			sys.exit(1)

		try:
			assert all(x != self.field.zero() for x in self.v)
		except AssertionError as e:
			print("El vector v no puede tener elementos nulos")
			sys.exit(1)	

		self.G = Matrix(self.field,[[self.v[j]*(self.alpha[j]**(i)) for j in range(n)] for i in range(k)]) #Una matriz generadora del código

	#Función que recibe un mensaje y devuelve la palabra código correspondiente	
	def codificar(self,mensaje,mensajeCuerpo=False):
		try:
			assert len(mensaje) == self.k
		except AssertionError as e:
			print("La longitud del mensaje debe coincidir con la dimensión del código")
			return ()
		
		#Hacemos lo mismo que para alpha
		mensajeField = None
		try:
			if is_prime(self.p) or mensajeCuerpo:
				mensajeField = [self.field(mensaje[i]) for i in range(len(mensaje))]
			else:
				mensajeField = procesarVectores(mensaje, self.a, self.subcuerpoPrimo)
		except:
			print("No ha sido posible procesar el mensaje")
			return ()	
		return vector(self.field,mensajeField)*self.G

	#Función que recibe una palabra código y devuelve el mensaje original, en caso de que no tenga errores	
	def decodificar(self,palabraCodigo, palabraCodigoCuerpo=False):
		try:
			assert len(palabraCodigo) == self.n
		except AssertionError as e:
			print("La longitud de la palabra de código debe coincidir con la longitud del código")
			return ()

		#Hacemos lo mismo que para alpha
		palabraCodigoField = None
		try:
			if is_prime(self.p) or palabraCodigoCuerpo:
				palabraCodigoField = [self.field(palabraCodigo[i]) for i in range(len(palabraCodigo))]
			else:
				palabraCodigoField = procesarVectores(palabraCodigo, self.a, self.subcuerpoPrimo)
		except:
			print("No ha sido posible procesar la palabra código")
			return ()	
		return vector(self.field,palabraCodigoField)*self.G.solve_right(identity_matrix(self.G.nrows()))	

	#Función que recibe una palabra código e introduce de forma aleatoria el número de errores indicado por numErrores
	@staticmethod #Se hace estático para que pueda ser llamado durante el cifrado en McEliece, sin que el usuario tenga que disponer de una instancia de la clase
	def introducirErrores(palabraCodigo,numErrores,n,field,subcuerpoPrimo,palabraCodigoCuerpo=False):
		try:
			assert len(palabraCodigo) == n
		except AssertionError as e:
			print("La longitud de la palabra de código debe coincidir con la longitud del código")
			return ()
		
		try:
			assert numErrores <= n and numErrores >= 0
		except AssertionError as e:
			print("El número de errores introducido no puede ser superior a la longitud del código ni negativo")
			return ()
		
		#Queremos que los elementos de errores sean no nulos
		errores = []
		while len(errores) < numErrores:
			error = field.random_element()
			if error != field.zero():
				errores.append(error)
		indicesAleatorios = random.sample(range(len(palabraCodigo)), numErrores)

		#Hacemos lo mismo que para alpha
		palabraCodigoField = None
		try:
			if palabraCodigoCuerpo:
				palabraCodigoField = [field(palabraCodigo[i]) for i in range(len(palabraCodigo))]
			else:
				palabraCodigoField = procesarVectores(palabraCodigo, field.gen(), subcuerpoPrimo)
		except:
			print("No ha sido posible procesar la palabra código")
			return ()		

		return [palabraCodigoField[i] + errores.pop() if i in indicesAleatorios else palabraCodigoField[i] for i in range(n)]


	#Función que recibe una palabra código y devuelve la palabra código corregida
	def corregirCodigo(self,palabraCodigo,palabraCodigoCuerpo=False):
		try:
			assert len(palabraCodigo) == self.n
		except AssertionError as e:
			print("La longitud de la palabra código debe coincidir con la longitud del código")
			return None,None
		
		#Hacemos lo mismo que para alpha
		palabraCodigoField = None
		try:
			if is_prime(self.p) or palabraCodigoCuerpo:
				palabraCodigoField = [self.field(palabraCodigo[i]) for i in range(len(palabraCodigo))]
			else:
				palabraCodigoField = procesarVectores(palabraCodigo, self.a, self.subcuerpoPrimo)
		except:
			print("No ha sido posible procesar la palabra código")
			return None,None	
		palabraCodigo = vector(self.field,palabraCodigoField)
		#Declaramos el anillo de polinomios sobre nuestro cuerpo finito
		R = PolynomialRing(self.field,'x')
		x = R.gen()

		#Declaramos el parámetro r
		r = self.n - self.k


		#El código no es capaz de corregir ningún error
		try:
			assert r > 0
		except AssertionError as e:
			return (None,list(palabraCodigo))

		#Declaramos el anillo cociente que usaremos en el cálculo del polinomio síndrome
		polinomioModulo = x**r
		RCociente = R.quotient(polinomioModulo,'z')

		#Construimos los polinomios nodales de Lagrange
		L = R(1)
		Li = [R(1) for j in range(len(self.alpha))]
		for i in range(len(self.alpha)):
			L *= (x - self.alpha[i])
			Li = [Li[j]*(x - self.alpha[i]) if j != i else Li[j]  for j in range(len(self.alpha))]

		#Construimos el vector u
		u = [(self.v[i]*Li[i](self.alpha[i]))**(-1) for i in range(len(self.alpha))]


		#Declaramos el polinomio síndrome
		S = RCociente(0)

		
		#Calculamos el polinomio síndrome
		for i in range(len(self.alpha)):
			S += RCociente(u[i] * palabraCodigo[i])*RCociente(1 - self.alpha[i]*x)**(-1)


		#Si el polinomio síndrome es 0, entonces la palabra código no contiene ningún error
		if S.is_zero():
			return (None,list(palabraCodigo))


		#Implementamos el algoritmo descrito en el Teorema 3.4.7
		a = polinomioModulo 
		b = S.lift()
		reminder = R(0)
		s_prev = R(1)
		t_prev = R(0)
		s_curr = R(0)
		t_curr = R(1)


		while True:
			quotient,reminder = a.quo_rem(b)

			s_new = s_prev - quotient * s_curr
			t_new = t_prev - quotient * t_curr

			if reminder.degree() < r/2:
				break

			a,b = b,reminder
			s_prev, t_prev = s_curr, t_curr
			s_curr, t_curr = s_new, t_new


		#Calculamos los polinomios localizador y evaluador de errores	
		if reminder.is_zero():
			sigma = t_new
			omega = 0
		else:
			sigma = t_new(0)**(-1)*t_new
			omega = t_new(0)**(-1)*reminder	
			
		
		#Calculamos las raíces del polinomio localizador
		raicesMultiplicidades = sigma.roots()
		raices = [raicesMultiplicidades[i][0] for i in range(len(raicesMultiplicidades))]

		#Calculamos las inversas de estas raíces y buscamos su índice en alpha, para así hallar el conjunto B
		raicesInversas = [raices[i]**(-1) for i in range(len(raices))]
		B = [i for i, raiz in enumerate(self.alpha) if raiz in raicesInversas]


		#Calculamos el vector de errores
		derivadaSigma = derivative(sigma,x)
		e = [(-self.alpha[i]*omega(self.alpha[i]**(-1)))*(u[i]*derivadaSigma(self.alpha[i]**(-1)))**(-1) if i in B else 0 for i in range(self.n)]

		#Si alguna de las componentes de alpha fuese 0, y el grado de sigma y omega coinciden, (punto 4 Proposición 3.4.4)
		#debemos calcular la correspondiente posición del vector de error de acuerdo
		#con el punto 5 de la Proposición 3.4.4
		indiceAlphaNulo = -1
		contador = 0
		for elemento in self.alpha:
			if elemento == self.field.zero():
				indiceAlphaNulo = contador
				break
			contador += 1

		if indiceAlphaNulo != -1 and omega.degree() == sigma.degree():
			e[indiceAlphaNulo] = omega.leading_coefficient()*(u[indiceAlphaNulo]**(-1))*(prod(-self.alpha[i] if i in B else 1 for i in range(self.n))**(-1))

		#Calculamos la palabra código que fue enviada
		c = [palabraCodigo[i] - e[i] for i in range(self.n)]
		return e,c
			


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--longitud', type=int, required=True, help='Longitud del código')
	parser.add_argument('--ordenCuerpo', type=int, required=True, help='Orden del cuerpo')
	parser.add_argument('--polinomioBase', nargs='+', type=int, help='Polinomio base del cuerpo')
	parser.add_argument('--ordenSubcuerpoPrimo', type=int, help='Orden del subcuerpo primo')
	parser.add_argument('--dimension', type=int, required=True, help='Dimensión del código')
	parser.add_argument('--alpha', nargs='+',required=True, help='Alpha')
	parser.add_argument('--v', nargs='+',required=True, help='V')
	parser.add_argument('--modo', type=str, required=True, help='Modo de ejecución')
	parser.add_argument('--mensaje', nargs='+',required=False, help='Mensaje')
	parser.add_argument('--palabraCodigo', nargs='+',required=False, help='Palabra de Código')
	parser.add_argument('--numErrores', type=int,required=False, help='Número de Errores')
	args = parser.parse_args()


	codigo = GRS(args.longitud,args.ordenCuerpo,args.polinomioBase,args.ordenSubcuerpoPrimo,args.dimension,args.alpha,args.v)


	if not is_prime(args.ordenCuerpo):
		if not args.polinomioBase:
			parser.error("Si el orden del cuerpo no es primo hay que pasar un polinomio base a utilizar")
			if not args.ordenSubcuerpoPrimo:
				parser.error("Si el orden del cuerpo no es primo hay que pasar el orden de su subcuerpo primo")	 

	if args.modo == "codificar":
		if not args.mensaje:
			parser.error("Es necesario pasar un mensaje en modo codificar")
		else:	
			print(codigo.codificar(args.mensaje))
	elif args.modo == "decodificar":
		if not args.palabraCodigo:
			parser.error("Es necesario pasar una palabra de código en modo decodificar")
		else:
			print(codigo.decodificar(args.palabraCodigo))	
	elif args.modo == "introducirErrores":
		if not args.palabraCodigo:
			parser.error("Es necesario pasar una palabra de código en modo introducirErrores")
		else:
			if not args.numErrores:
				parser.error("Es necesario pasar el número de errores en modo introducirErrores")
			else:
				print(GRS.introducirErrores(args.palabraCodigo,args.numErrores, codigo.n, codigo.field, codigo.subcuerpoPrimo, is_prime(args.ordenCuerpo)))
	elif args.modo == "corregirCodigo":
		if not args.palabraCodigo:
			parser.error("Es necesario pasar una palabra código en modo corregirCodigo")
		else:
			e,c = codigo.corregirCodigo(args.palabraCodigo)
			print("Vector de error: " + str(e))
			print("Palabra código corregida: " + str(c))
	else:
		print("El modo introducido no está disponible")
		
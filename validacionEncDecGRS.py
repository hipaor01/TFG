import argparse
import time
from encDecGRS import GRS
from sage.all import *
import random

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--numIteraciones', type=int, required=True)
    args = parser.parse_args()

    #El orden de los parámetros es longitud,orden del cuerpo, polinomio base, subcuerpo primo, dimensión
    parametros = [[20,81,[1,0,0,1,2],3,8],[60,127,None,None,28],[80,149,None,None,30], [96,193,None,None,45], [128,251,None,None,64], [160,257,None,None,80], [200,281,None,None,100], [240,293,None,None,120]]
    
    #Contador para que en cada combinación de parámetros registremos el número de casos exitosos
    numExitos = 0
    #Variable para ir acumulando los tiempos de ejecución y luego hacer la media para cada combinación de parámetros
    tiempoMedio = 0

    for parametro in parametros:
        numExitos = 0
        tiempoMedio = 0
        #Construimos el cuerpo finito sobre el que vamos a trabajar para poder generar los vectores alpha y v
        field = None
        polinomioModulo = None
        if not is_prime(parametro[1]):
            subcuerpoPrimo = GF(parametro[3])
            R = PolynomialRing(subcuerpoPrimo,'x')
            polinomioModulo = R(parametro[2])
            field = GF(parametro[1], name='a', modulus=polinomioModulo)
        else:
            field = GF(parametro[1]) 

        for i in range(args.numIteraciones):
            #Generamos aleatoriamente el vector alpha, teniendo en cuenta que no puede haber entradas repetidas
            alpha = set()
            while len(alpha) < parametro[0]:
                elemento = field.random_element()
                alpha.add(elemento)
            alpha = list(alpha)

            #Generamos aleatoriamente el vector v, teniendo en cuenta que no puede tener entradas nulas
            v = []
            while len(v) < parametro[0]:
                elemento = field.random_element()
                if elemento != 0:
                    v.append(elemento)

            #Creamos una instancia de la clase GRS
            codigo = GRS(parametro[0],parametro[1],parametro[2],parametro[3],parametro[4],alpha,v,True,True)

            #Generamos aleatoriamente un mensaje de longitud la dimensión del código
            mensaje = random_vector(codigo.field, codigo.k)

            #Iniciamos medición del tiempo
            inicio = time.perf_counter()

            #Codificamos el mensaje
            codificado = codigo.codificar(mensaje,True)

            #Generamos un número aleatorio de errores a introducir
            numErrores = random.randint(0, int((codigo.n-codigo.k)/2))

            #Introducimos el número de errores dado por numErrores
            codificadoErrores = GRS.introducirErrores(codificado,numErrores,codigo.n,codigo.field,codigo.subcuerpoPrimo,True)

            #Corregimos los errores introducidos en codificadoErrores, obteniendo la palabra corregida y el vector de errores introducidos
            e,codificadoCorregido = codigo.corregirCodigo(codificadoErrores,True)

            #Finalmente decodificamos codificadoCorregido y deberíamos obtener mensaje
            mensajeCandidato = codigo.decodificar(codificadoCorregido, True)

            #Finalizamos medición del tiempo
            final = time.perf_counter()

            if mensaje == mensajeCandidato:
                numExitos += 1
                      
            tiempoMedio += (final - inicio)

        print("Longitud del código: " + str(parametro[0]))
        print("Orden del cuerpo: " + str(parametro[1]))
        if parametro[2] is not None:
            print("Polinomio base: " + str(polinomioModulo))
        if parametro[3] is not None:
            print("Oden Subcuerpo primo: " + str(parametro[3]))
        print("Dimensión del código: " + str(parametro[4]))
        print("Número de iteraciones totales: " + str(args.numIteraciones))
        print("Número de éxitos: " + str(numExitos))
        print("Porcentaje de aciertos: " + str(100 * (numExitos/args.numIteraciones)))
        print(f"Tiempo medio de ejecución: {tiempoMedio/args.numIteraciones:.6f} segundos" + "\n")
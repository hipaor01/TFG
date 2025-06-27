import argparse
from SidelnikovShestakov import ataqueSidelnikovShestakov
from sage.all import *
import time

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
        for i in range(args.numIteraciones):
            #Iniciamos medición del tiempo
            inicio = time.perf_counter()
            pt, ptCandidato = ataqueSidelnikovShestakov(parametro[0], parametro[1], parametro[2], parametro[3], parametro[4])
            #Finalizamos medición del tiempo
            final = time.perf_counter()
            if(pt == ptCandidato):
                numExitos += 1
            tiempoMedio += (final - inicio)
        print("Longitud del código: " + str(parametro[0]))
        print("Orden del cuerpo: " + str(parametro[1]))
        if parametro[2] is not None:
            subcuerpoPrimo = GF(parametro[3])
            R = PolynomialRing(subcuerpoPrimo,'x')
            polinomioModulo = R(parametro[2])
            print("Polinomio base: " + str(polinomioModulo))
        if parametro[3] is not None:
            print("Subcuerpo primo: " + str(parametro[3]))
        print("Dimensión del código: " + str(parametro[4]))
        print("Número de iteraciones totales: " + str(args.numIteraciones))
        print("Número de éxitos: " + str(numExitos))
        print("Porcentaje de aciertos: " + str(100 * (numExitos/args.numIteraciones)))
        print(f"Tiempo medio de ejecución: {tiempoMedio/args.numIteraciones:.6f} segundos" + "\n")

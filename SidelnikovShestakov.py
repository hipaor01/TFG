from sage.all import *
from encDecGRS import GRS
from McElieceGRS import Gen,Enc
from itertools import permutations
import argparse
import copy

#Función que nos va a ayudar a calcular los alpha faltantes tras obtener alpha_j, para k+1<=j<=n
def calcularAlphas(k,EB,indice,candidatosAlpha,field):
    #Inicializamos la lista de alphas para indice, con el mismo tamaño que el número de candidatos (un valor None indicará que no se ha conseguido calcular
    #alpha_indice para ese candidato, y que por lo tanto podemos descartarlo)
    alphas = [None for i in range(len(candidatosAlpha))]

    contador = 0
    for candidato in candidatosAlpha:
        #Resolvemos el siguiente sistema para cada candidato, donde las incógnitas son cindice = cb1/cbindice y sindice = ci*alpha_indice
        A = Matrix(field,2,2)
        A[0,0] = candidato[0]
        A[0,1] = -1
        A[1,0] = candidato[1]
        A[1,1] = -1
        b = vector(field,[(EB[0,k]/EB[indice,k])*A[0,0], (EB[0,k+1]/EB[indice,k+1])*A[1,0]])
        soluciones = A.solve_right(b)
        #Deshacemos el cambio de variable, de manera que obtenemos los alpha mediante el cociente soluciones[1]/soluciones[0], siempre que 
        #soluciones[0] no sea nulo. Además, también ignoramos aquellos alphas que se encuentren ya en el candidato
        if soluciones[0] != field.zero():
            alpha = soluciones[1]/soluciones[0]
            if alpha not in candidato:
                alphas[contador] = alpha
        contador +=1       

    return alphas

#Función que va a ejecutar una simulación del ataque a partir de los parámetros dados
def ataqueSidelnikovShestakov(n, p, polinomioBaseArgumento, ordenSubcuerpoPrimo, k):

    if k < 2 or k > n-2:
        print("La dimensión debe ser mayor o igual que 2 y menor o igual que la longitud menos 2")
        sys.exit(1)

    #No puede pasarse un cuerpo negativo
    if p <= 0:
        print("El cuerpo debe ser de orden positivo")
        sys.exit(1)

    #Connstruimos un polinomio a partir de los coeficientes pasados por línea de comandos
    polinomioBase = None
    cuerpoPrimo = None
    if(polinomioBaseArgumento and ordenSubcuerpoPrimo):
        if not is_prime(ordenSubcuerpoPrimo):
            print("El orden del subcuerpo primo debe ser primo")
            sys.exit(1)
        try:
            cuerpoPrimo = GF(ordenSubcuerpoPrimo)
            R = PolynomialRing(cuerpoPrimo, 'x')
            x = R.gen()
            coeficientes = [cuerpoPrimo(coefi) for coefi in polinomioBaseArgumento]
            polinomioBase = sum(coefi * x**i for i,coefi in enumerate(reversed(coeficientes)))
        except:
            print("No ha sido posible procesar el polinomio pasado")
            sys.exit(1)

    #Si p no es primo, significa que queremos trabajar sobre una extensión de un cuerpo primo, por lo que usaremos ordenSubcuerpoPrimo y polinomioBase
    field = None
    if not is_prime(p):
        try:
            R = PolynomialRing(cuerpoPrimo, 'x')
            polinomioModulo = R(polinomioBase)
            field = GF(p, name='a', modulus=polinomioModulo)
        except:
            print("El polinomio y subcuerpo primo no pueden construir el cuerpo deseado")
            sys.exit(1)
    else:
        field = GF(p)

    #Generamos un texto plano aleatorio, para luego comprobar que somos capaces de descifrarlo
    pt = random_vector(field,k)

    #Generamos las claves del criptosistema de McEliece
    pk,_= Gen(n,k,p,ordenSubcuerpoPrimo,polinomioBase)
    #Este va a ser el punto de partida del ataque
    B = pk[0]
    
    #Generamos el texto cifrado
    ct = Enc(pk,pt,True)

    #En primer lugar, calculamos la forma escalonada de la matriz B, EB
    EB = B.echelon_form()

   
    #Calculamos los cocientes b1,j/b2,j
    bj = []
    for j in range(k,n):
        bj.append(EB[0,j]/EB[1,j])

    '''
    Como el cociente cb1/cb2 tiene que ser "adivinado", debemos hacer una fuerza bruta, recorriendo todos los elementos del cuerpo menos
    los contenidos en bj, para evitar divisiones por 0, y el 0, ya que todos los cbi son no nulos
    Guardamos en candidatosAlpha todas las listas resultantes de las distintas combinaciones
    '''
    candidatosAlpha = []
    for elemento in field:
        if elemento in bj or elemento == field.zero():
            continue
        listaAux = [elemento/(elemento - bj[j]) for j in range(len(bj))]
        candidatosAlpha.append(listaAux)

    #Vamos a eliminar aquellos candidatos donde haya duplicados, ya que no pueden formar un vector alpha de un código GRS
    candidatosAlpha = [candidato for candidato in candidatosAlpha if len(candidato) == len(set(candidato))]

    #Realizamos una copia para ir añadiendo los elementos que vamos encontrando sin que afecte al cálculo de los que quedan
    copiaCandidatosAlpha = copy.deepcopy(candidatosAlpha)

    #Añadimos alpha_1 = 0 y alpha_2 = 1 al comienzo de todos los candidatos
    for candidato in copiaCandidatosAlpha:
        candidato.insert(0,field(0))
        candidato.insert(1,field(1))        

    #Calculamos ahora los alpha restantes, a partir del índice 2
    for i in range(2,k):
        alphas = calcularAlphas(k,EB,i,candidatosAlpha,field)
        contador = 0
        for candidato in copiaCandidatosAlpha:
            candidato.insert(i,alphas[contador])
            contador += 1

    copiaCandidatosAlpha = [candidato for candidato in copiaCandidatosAlpha if len(candidato) == len(set(candidato)) and all(x is not None for x in candidato)]
    candidatosAlpha = copiaCandidatosAlpha

    #En la matriz B1 vamos a guardar las primeras k+1 columnas de la matriz B
    B1 = Matrix(field,k,k+1)
    B1[:] = B[:,0:k+1]
    #Resolvemos el sistema B1*c = 0, lo que es equivalente a encontrar el núcleo de B1, y nos quedamos con una base suya
    c = B1.right_kernel().basis()[0]
    
    #De la candidata a G sabemos de momento los posibles alphas, pero desconocemos los vs. Sin embargo, podemos formar un sistema
    #similar al anterior cambiando los cs obtenidos en el sistema anterior por los vs, quedando los vs como incógnita. En G1
    #almacenaremos las primeras k+1 columnas de G, con los cambios indicados
    for candidato in candidatosAlpha:
        G1 = Matrix(field, [[c[j]*(candidato[j]**i) for j in range(k+1)] for i in range(k)])
        vCandidato = G1.right_kernel().basis()[0]
        vCandidato = list(vCandidato)
        #Recalculamos G1
        G1 = Matrix(field, [[vCandidato[j]*(candidato[j]**i) for j in range(k+1)] for i in range(k)])
        #Podemos calcular ya las candidatas a S
        Scandidata = B[:,0:k]*(G1[:,0:k]**(-1))
        #Y finalmente la candidata a G
        Gcandidata = (Scandidata**(-1))*B
        #Como conocemos toda la primera fila de la candidata a G, podemos reconstruir el resto de posiciones de v
        vCandidato += [Gcandidata[0,j] for j in range(k+1, n)]
        #Si alguna de las posiciones de vCandidato es nula, lo descartamos
        if any(v == field.zero() for v in vCandidato):
            continue
        #Construimos un código GRS  a partir de los candidatos
        codigo = GRS(n,p,polinomioBase,ordenSubcuerpoPrimo,k,candidato,vCandidato,True,True)

        #Si salta una excepción al llamar a corregirCodigo, descartamos ese candidato
        try:
            e, correccion = codigo.corregirCodigo(ct,True)
        except Exception as e:
            continue

        #Decodificamos correccion
        ptCandidato = codigo.decodificar(correccion, True)
        #Finalmente multiplicamos por la inversa de SCandidata
        ptCandidato = ptCandidato*(Scandidata**(-1))
        break

    return pt,ptCandidato

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--longitud', type=int, required=True, help='Longitud del código')
    parser.add_argument('--ordenCuerpo', type=int, required=True, help='Orden del cuerpo')
    parser.add_argument('--polinomioBase', nargs='+', type=int, help='Polinomio base del cuerpo')
    parser.add_argument('--ordenSubcuerpoPrimo', type=int, help='Orden del subcuerpo primo')
    parser.add_argument('--dimension', type=int, required=True, help='Dimensión del código')
    args = parser.parse_args()

    pt,ptCandidato = ataqueSidelnikovShestakov(args.longitud, args.ordenCuerpo, args.polinomioBase, args.ordenSubcuerpoPrimo, args.dimension)
    
    print("Texto plano original: ", pt)
    print("Texto plano descifrado: ", ptCandidato)


    

    
    

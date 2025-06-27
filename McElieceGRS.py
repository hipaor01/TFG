from encDecGRS import GRS
from sage.all import *
from sage.combinat.permutation import Permutations
import argparse
import random

def Gen(n,k,p,ordenSubcuerpoPrimo, polinomioBase):
    #No puede pasarse un cuerpo negativo
    if p <= 0:
        print("El cuerpo debe ser de orden positivo")
        sys.exit(1)


    #Antes de llamar al constructor de GRS, definimos el cuerpo finito sobre el que se quiere trabajar para generar los vectores alpha y v
    #Si p no es primo, significa que queremos trabajar sobre una extensión de un cuerpo primo, por lo que usaremos ordenSubcuerpoPrimo y polinomioBase
    field = None
    if not is_prime(p):
        try:
            R = PolynomialRing(GF(ordenSubcuerpoPrimo), 'x')
            polinomioModulo = R(polinomioBase)
            field = GF(p, name='a', modulus=polinomioModulo)

            assert ordenSubcuerpoPrimo**(polinomioModulo.degree()) == p
        except:
            print("El polinomio y subcuerpo primo no pueden construir el cuerpo deseado")
            sys.exit(1)
    else:
        field = GF(p)

    #Para evitar que haya un bucle infinito en la generación de alpha, hay que comprobar que el orden del cuerpo sea mayor o igual que n
    if p < n : 
        print("El orden del cuerpo debe ser mayor o igual que la longitud")
        sys.exit(1)

    #Generamos aleatoriamente el vector v, teniendo en cuenta que no puede tener entradas nulas
    v = []
    while len(v) < n:
        elemento = field.random_element()
        if elemento != 0:
            v.append(elemento)

    #Generamos aleatoriamente el vector alpha, teniendo en cuenta que no puede haber entradas repetidas
    alpha = set()
    while len(alpha) < n:
        elemento = field.random_element()
        alpha.add(elemento)
    alpha = list(alpha)

    #Ya tenemos todos los parámetros para definir un código GRS
    codigo = GRS(n,p,polinomioBase,ordenSubcuerpoPrimo,k,alpha,v,True,True)
    if codigo is None:
        return

    #Generamos aleatoriamente la matriz S, que debe ser invertible
    S = random_matrix(field, k, algorithm='unimodular')

    #Generamos aleatoriamente la matriz P, que debe ser de permutación
    P = Matrix(field,Permutations(n).random_element().to_matrix())

    #Generamos la clave pública, junto con el cuerpo sobre el que estamos trabajando y el orden de su subcuerpo primo
    pk = [S * codigo.G * P,field,ordenSubcuerpoPrimo]

    #Generamos la clave privada
    sk = [codigo,S,P]

    return pk,sk

def Enc(pk,pt,ptCuerpo=False):
    #Desglosamos la clave pública en la matriz generadora "camuflada" y cuerpo finito sobre el que trabajamos
    G = pk[0]
    field = pk[1]
    ordenSubcuerpoPrimo = pk[2]
    subcuerpoPrimo = None
    if ordenSubcuerpoPrimo:
        subcuerpoPrimo = GF(ordenSubcuerpoPrimo)

    #A partir de las dimensiones de G, obtenemos la longitud y dimensión del código
    k = G.nrows()
    n = G.ncols()

    #Comprobamos que el texto plano pt es de longitud k
    try:
        assert len(pt) == k
    except AssertionError as e:
        print("La longitud del texto plano debe ser " + str(k))
        sys.exit(1)
    #Convertimos pt en un vector para poder operar con él
    #Si el cuerpo no es primo, tendremos que hacer un procesado previo de lo enviado por línea de comandos
    vectorField = []

    try:
        if is_prime(field.order()) or ptCuerpo:
            pt = [field(pt[i]) for i in range(len(pt))]
        else:
            for elemento in pt:
                coeficientesBase = [subcuerpoPrimo(int(coefi)) for coefi in elemento.split(',')] #Los coeficientes deben ser elementos del subcuerpo primo
                elementoField = sum(coefi * field.gen()**i for i,coefi in enumerate(reversed(coeficientesBase)))
                vectorField.append(elementoField)
            pt = vectorField
    except:
        print("No ha sido posible procesar el texto plano")
        sys.exit(1)
        
    pt = vector(field, pt)

    #Obtenemos el texto cifrado ct introduciendo hasta (n-k)/2 errores sobre ptG
    numErrores = random.randint(0, int((n-k)/2))
    ct = GRS.introducirErrores(pt*G, numErrores, n, field, subcuerpoPrimo, palabraCodigoCuerpo=True)

    return ct  

def Dec(sk,ct):
    #Desglosamos la clave privada en la instancia de GRS, matriz invertible y de permutación
    codigo = sk[0]
    S = sk[1]
    P = sk[2]
    n = codigo.n 
    field = codigo.field

    #Comprobamos que el texto cifrado ct es de longitud n
    try:
        assert len(ct) == n
    except AssertionError as e:    
        print("La longitud del texto cifrado debe ser " + str(n))
    #Convertimos ct en un vector para poder operar con él
    ct = [field(ct[i]) for i in range(len(ct))]
    ct = vector(field, ct)    
    
    #Corregimos los errores introducidos durante el cifrado
    errorPalabraCodigo = codigo.corregirCodigo(ct * P.inverse(), True)
    #Lo que hemos obtenido en el paso anterior es pt'G, por lo que hay que decodificarlo en pt'
    pt = codigo.decodificar(errorPalabraCodigo[1], True)

    #Finalmente, obtenemos pt multiplicando a la derecha por la inversa de S
    pt = pt * S.inverse()

    return list(pt)

#En esta función, realizamos el cifrado del texto plano dado por pt, lo desciframos y devolvemos el resultado de dicho descifrado
def CifradoDescifrado(n,p,polinomioBaseArgumento,ordenSubcuerpoPrimo,k,pt):
    #Construimos un polinomio a partir de los coeficientes pasados en polinomioBaseArgumento
    polinomioBase = None
    if(polinomioBaseArgumento and ordenSubcuerpoPrimo):
        if not is_prime(ordenSubcuerpoPrimo):
            print("El orden del subcuerpo primo debe ser primo")
            sys.exit(1)
        cuerpoPrimo = GF(ordenSubcuerpoPrimo)

        try:
            R = PolynomialRing(cuerpoPrimo, 'x')
            x = R.gen()
            coeficientes = [cuerpoPrimo(coefi) for coefi in polinomioBaseArgumento]
            polinomioBase = sum(coefi * x**i for i,coefi in enumerate(reversed(coeficientes)))
        except Exception:
            print("No ha sido posible procesar el polinomio pasado")
            sys.exit(1)

    #Generamos la clave pública y privada
    pk,sk = Gen(n, k, p, ordenSubcuerpoPrimo, polinomioBase)
    if pk is None or sk is None:
        sys.exit(1)

    #Si pt es nulo, generamos desde aquí pt de manera aleatoria ya que disponemos del cuerpo finito usado
    #Además, como las componentes de pt estarán ya en el cuerpo, ponemos ptCuerpo = True
    ptCuerpo = False
    if pt is None:
        pt = random_vector(pk[1], k)
        ptCuerpo = True

    #Ciframos el mensaje recibido usando la clave pública
    ct = Enc(pk, pt, ptCuerpo)
    
    #Desciframos el texto cifrado usando la clave privada
    ptCandidato = Dec(sk, ct)

    #Devolvemos también pt por si se hubiera generado dentro de esta función
    return ct,ptCandidato,pt
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--longitud', type=int, required=True, help='Longitud del código')
    parser.add_argument('--ordenCuerpo', type=int, required=True, help='Orden del cuerpo')
    parser.add_argument('--polinomioBase', nargs='+', type=int, help='Polinomio base del cuerpo')
    parser.add_argument('--ordenSubcuerpoPrimo', type=int, help='Orden del subcuerpo primo')
    parser.add_argument('--dimension', type=int, required=True, help='Dimensión del código')
    parser.add_argument('--mensaje', nargs='+',required=True, help='Mensaje')
    args = parser.parse_args()

    ct,ptCandidato,_ = CifradoDescifrado(args.longitud,args.ordenCuerpo,args.polinomioBase,args.ordenSubcuerpoPrimo,args.dimension,args.mensaje)

    print("El texto cifrado es: " + str(ct))
    print("El texto descifrado es: " + str(ptCandidato))
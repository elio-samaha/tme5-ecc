# Sorbonne Université LU3IN024 2021-2022
# TME 5 : Cryptographie à base de courbes elliptiques
#
# Etudiant.e 1 : Samaha Elio 21105733
# Etudiant.e 2 : NOM ET NUMERO D'ETUDIANT

from math import sqrt
import matplotlib.pyplot as plt
from random import randint

# Fonctions utiles

def exp(a, N, p):
    """Renvoie a**N % p par exponentiation rapide."""
    def binaire(N):
        L = list()
        while (N > 0):
            L.append(N % 2)
            N = N // 2
        L.reverse()
        return L
    res = 1
    for Ni in binaire(N):
        res = (res * res) % p
        if (Ni == 1):
            res = (res * a) % p
    return res


def factor(n):
    """ Return the list of couples (p, a_p) where p is a prime divisor of n and
    a_p is the p-adic valuation of n. """
    def factor_gen(n):
        j = 2
        while n > 1:
            for i in range(j, int(sqrt(n)) + 1):
                if n % i == 0:
                    n //= i
                    j = i
                    yield i
                    break
            else:
                if n > 1:
                    yield n
                    break

    factors_with_multiplicity = list(factor_gen(n))
    factors_set = set(factors_with_multiplicity)

    return [(p, factors_with_multiplicity.count(p)) for p in factors_set]


def inv_mod(x, p):
    """Renvoie l'inverse de x modulo p."""
    return exp(x, p-2, p)


def racine_carree(a, p):
    """Renvoie une racine carrée de a mod p si p = 3 mod 4."""
    assert p % 4 == 3, "erreur: p != 3 mod 4"

    return exp(a, (p + 1) // 4, p)


# Fonctions demandées dans le TME

def est_elliptique(E):
    """
    Renvoie True si la courbe E est elliptique et False sinon.

    E : un triplet (p, a, b) représentant la courbe d'équation
    y^2 = x^3 + ax + b sur F_p, p > 3
    """
    p , a , b = E
    return (-16*(4*exp(a , 3 , p) + 27 * exp(b , 2 , p)) )% p != 0


def point_sur_courbe(P, E):
    """Renvoie True si le point P appartient à la courbe E et False sinon."""
    if P == () : return True
    p , a , b = E
    x , y = P
    return exp(y , 2 , p) == ( (exp(x , 3 , p) + a*x + b) % p)


def symbole_legendre(a, p):
    """Renvoie le symbole de Legendre de a mod p."""
    
    return exp(a , (p-1)//2 , p) #car vaut 0 si p divise a : direct.   
    #Vaut 1 si il existe b t.q a = b**2 mod p car par fermat b**(2 * (p - 1)/2) = b ** (p-1) = 1 mod p (fermat p premier)
    #Vaut -1 sinon car supposons a ** ((p-1)/2) != -1 ---> a ** 2((p-1)/2) != (-1)**2 != 1 mod p (ce qui contredit le theoreme de fermat d'ou si a n est ni un carre ni un multiple de p alors il ne reste que a n est pas un carre qui est qssocie a la valeur -1 )

def def_value(): 
    return set()
      
from collections import defaultdict 

def cardinal(E):
    """Renvoie le cardinal du groupe de points de la courbe E."""
    res = 1
    (p , a , b) = E
    dic = defaultdict(def_value)
    dic[0].add(0)
    for i in range(1 , p//2 + 1):
        y2 = exp(i , 2 , p)
        dic[y2].add(i) 
        dic[y2].add(-i) 

    for j in range(p):
        y2 = ((exp(j , 3 , p) + a*j + b) % p)
        if y2 in dic:
            res += len(dic[y2])

    return res


def liste_points(E):
    """Renvoie la liste des points de la courbe elliptique E."""
    p, a, b = E

    assert p % 4 == 3, "erreur: p n'est pas congru à 3 mod 4." #a**2((p+1)/4) = a**2((p-1+2)/4) = a**((p-1)/2) * a = 1 * a mod p = a mod p (avant derniere egalité par quest 3) 
    #(p+1)/4 est un entier car p est congru a 3 mod 4 donc il existe k >0 t.q p = 4k + 3 d'ou p + 1 = 4k + 4 ---> (p+1)/4 = k + 1

    li = [()]
    for x in range(p):
        y2 = ((exp(x , 3 , p) + a*x + b) % p)
        
        if y2 == 0:
            li.append((x , 0))
        elif symbole_legendre(y2 , p) == 1:
            rac = exp(y2 , (p+1)//4 , p)
            if rac == 0:
                li.append(x , rac)
            else:
                li.append((x , p - rac))        
                li.append((x , rac))

    return li

#Q6 : Soit p un nombre premier, p +1−2√p⩽|Ea,b| ⩽ p+1+2√p

import math
from itertools import product

def cardinaux_courbes(p):
    """
    Renvoie la distribution des cardinaux des courbes elliptiques définies sur F_p.

    Renvoie un dictionnaire D où D[i] contient le nombre de courbes elliptiques
    de cardinal i sur F_p.
    """
    D = defaultdict(lambda : 0)
    
    for a , b in product(range(p) , range(p)): #on essaye tous les couples (a,b) possibles et on calcule le cardinale de E(a,b)
        if est_elliptique((p , a , b)):
            card = cardinal((p , a , b))
            D[card] += 1
    
    return D


def dessine_graphe(p):
    """Dessine le graphe de répartition des cardinaux des courbes elliptiques définies sur F_p."""
    bound = int(2 * sqrt(p))
    C = [c for c in range(p + 1 - bound, p + 1 + bound + 1)]
    D = cardinaux_courbes(p)

    plt.bar(C, [D[c] for c in C], color='b')
    plt.show()

#dessine_graphe(5)

def moins(P, p):
    """Retourne l'opposé du point P mod p."""
    if P == ():
        return () 
    p1 , p2 = P
    return (p1 , p - p2)


def est_egal(P1, P2, p):
    """Teste l'égalité de deux points mod p."""

    if P1 == ():
        return P2 == ()
    if P2 == ():
        return P1 == ()

    x1 , y1 = P1
    x2 , y2 = P2

    return ((x1 % p) == (x2 % p)) and ((y1 % p) == (y2 % p)) 


def est_zero(P):
    """Teste si un point est égal au point à l'infini."""

    return P == ()


def addition(P1, P2, E):
    """Renvoie P1 + P2 sur la courbe E."""
    p , a , b = E
    
    if est_zero(P1):
        return P2
    if est_zero(P2):
        return P1

    if est_egal(P1 , moins(P2 , p) , p):
        return ()

    x1 , y1 = P1
    x2 , y2 = P2

    if est_egal(P1 , P2, p):
        lbda = ((3 * exp(x1 , 2 , p) + a) * inv_mod(2*y1 , p)) % p
    else :
        lbda = ((y2 - y1) * inv_mod(x2 - x1 , p)) % p
    
    x3 = (pow(lbda , 2 , p) - x1 - x2) % p
    y3 = (lbda * (x1- x3) - y1) % p

    return (x3 , y3)

def multiplication_scalaire_aux(k, P, E, lookup): #memoization inutile mais bon

    p , a , b = E

    if k in lookup : return lookup[k]

    if k in lookup:
        return lookup[k]

    if k < 0 :
        lookup[k] = moins(multiplication_scalaire_aux(-k , P , E , lookup) , p)
        return lookup[k]

    if k % 2 == 0: #k = 2k' --> res = k' * P + k' * P
        k_sur_2_P = multiplication_scalaire_aux(k // 2, P, E, lookup)
        lookup[k] = addition(k_sur_2_P, k_sur_2_P, E)
    else: #k = 2k' + 1 --> res = k' * P + k' * P + P
        k_sur_2_P = multiplication_scalaire_aux(k // 2, P, E, lookup)
        even_part = addition(k_sur_2_P, k_sur_2_P, E)
        lookup[k] = addition(even_part, P, E)

    lookup[k//2] = k_sur_2_P

    return lookup[k]

def multiplication_scalaire(k, P, E):
    """Renvoie la multiplication scalaire k*P sur la courbe E."""
    p , _ , _ = E
    lookup = {0 : ()}
    return multiplication_scalaire_aux(k, P, E , lookup)



def ordre(N, factors_N, P, E):
    """Renvoie l'ordre du point P dans les points de la courbe E mod p. 
    N est le nombre de points de E sur Fp.
    factors_N est la factorisation de N en produit de facteurs premiers."""
    p , a , b = E
    if P == () : return 1

    power_ranges = []

    for p, vp in factors_N:     
        power_ranges.append(range(vp + 1)) #on cree une liste de liste chacune contenant [0 , vp] representant les choix de puissances qu on a pour chaque p

    divisor_powers = list(product(*power_ranges)) #on prend toute les combinaisons de puissance

    div = []

    for powers in divisor_powers:
        divisor = 1
        for (p, power) in zip(factors_N, powers): #on parcours nos combinaison en tenant compte du comptage de chaque puissance pour chaque p grace au zip 
            prime = p[0]
            divisor *= prime ** power   #on cree un nouveau diviseur en faisant un nouveau produit de sous diviseurs qui sont tous premiers
        div.append(divisor)

    div.sort()  #on commence par le plus petit vers le plus grand et on s'arrete des que [d]P s annule et on renvoie le diviseur "d" correspondant soit l'ordre de l'element P.
    for d in div :
        if est_zero(multiplication_scalaire(d , P , E)):
            return d

    return -1 #controle d erreur


def point_aleatoire_naif(E):
    """Renvoie un point aléatoire (différent du point à l'infini) sur la courbe E."""
    
    return 


def point_aleatoire(E):
    """Renvoie un point aléatoire (différent du point à l'infini) sur la courbe E."""

    return


def point_ordre(E, N, factors_N, n):
    """Renvoie un point aléatoire d'ordre N sur la courbe E.
    Ne vérifie pas que n divise N."""

    return

def keygen_DH(P, E, n):
    """Génère une clé publique et une clé privée pour un échange Diffie-Hellman.
    P est un point d'ordre n sur la courbe E.
    """
    sec = None # A remplacer
    pub = None # A remplacer
    
    return (sec, pub)

def echange_DH(sec_A, pub_B, E):
    """Renvoie la clé commune à l'issue d'un échange Diffie-Hellman.
    sec_A est l'entier secret d'Alice et pub_b est l'entier public de Bob."""

    return

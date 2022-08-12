from itertools import combinations, product
import numpy as np
from multiprocessing import Pool
from multiprocessing import cpu_count
import copy
import yaml



def operation(x,y):
    """
    Vectorlike merger
    
    Return
    ------
    gcd: greatest common divisor
    result: result of the operation, simplified and sorted from smallest to largest
    """
    op = np.sum(x * y**2) * x - np.sum(x**2 * y) * y
    op = np.sort(op)  ##Ordering the solution (both chiral and vectorlike)
    gcd = np.gcd.reduce(op)  ##Greatest common divisor of the solution
    
    if np.all((op == 0)): ##To avoid a nan array
        result = op
    else:
        result = op/gcd
    return gcd, result
    
def m_value(n):
    """
    m value for both even and odd cases
    """
    if n%2==0:
        m = n/2 -1
    else:
        m = (n-3)/2
    return m
    
def eq_satisfied(x, gcd):
    """
    This function verifies that equations (z1³+z2³+ ... + zn³ = 0) and (z1+z2+ ... + zn = 0) are 
    satisfied at the same time, and if it is true, then return the input parameters
    
    Parameters
    ----------
    x: it correspond to the solutions z1, z2, ..., zn
    gcd: the greatest common divisor
    
    Return
    ------
    x: it correspond to the solutions z1, z2, ..., zn
    gcd: the greatest common divisor
    """
    eq_3 = np.sum(x**3)
    eq_1 = np.sum(x)
    
    if eq_3 ==0 and eq_1 ==0:
        #print(f'{x} satisfies the equations with gcd = {gcd}')
        return gcd, x
        
def apply_merger(l, k, n):
    """
    Create the vectorlike for both even and odd cases, and apply the merger operation
    
    Parameters
    ----------
    l: array for construct the vectorlike
    k: array for construct the vectorlike
    n: number of elements in the solution (Z)
    """
    if n%2==0:

        v_p = np.append([l[0]], k)
        v_n = np.append([0, 0], l)

        v_p = np.append(v_p, -np.append([l[0]], k))
        v_n = np.append(v_n, -np.array(l))
        z_sol = operation(v_p, v_n)  ##Solution
    
    ##Odd
    else:

        u_p = np.append([0], k)
        u_n = np.append(l, [k[0], 0])

        u_p = np.append(u_p, -np.array(k))
        u_n = np.append(u_n, -np.append(l, [k[0]]))
        z_sol = operation(u_p, u_n) ##Solution
        
    z_sol = eq_satisfied(z_sol[1], z_sol[0])
    
    return z_sol
    
    
def chiral_solution(n):
    m = int(m_value(n))
    
    ##A vectorlike can contain elements between (-30,30)
    sample_list = np.arange(-30,30)
    
    ##Even case
    if n%2==0:
    ##Se eligen "m" valor de la lista que contiene enteros entre (-30,30)
    ##Para esta elección se usa una combinatoria de modo que no se tengan en cuenta casos repetidos, por ejemplo, [1,2] y [2,1]
        ls = np.array(list(combinations(sample_list, m)))
        ks = np.array(list(combinations(sample_list, m)))

    else:
    ##Se eligen "m" y "m+1" valores de la lista que contiene enteros entre (-30,30)
        ls = np.array(list(combinations(sample_list, m)))
        ks = np.array(list(combinations(sample_list, m+1)))

    ##Once the "ls" and "ks" lists (list of lists) have been constructed, with all the 
    ##possible combinations for a specific "m" value, the next step is to create the 
    ##combination between the "ls" and "ks" lists. For this, the indexes or the number 
    ##of lists within the "ls" and "ks" lists are taken into account. These indexes are 
    ##the ones that are now going to be combined using meshgrid
    index_ls = np.arange(len(ls))
    index_ks = np.arange(len(ks))
    mesh = np.array(np.meshgrid(index_ls, index_ks))
    
    ##It is a list of lists. Each list corresponds to an index for "l" and another for "k", 
    ##and contains all possible combinations.
    combs = mesh.T.reshape(-1, 2)
    
    
    ##Para almacenar los z y el valor del máximo común divisor
    z_results = []
    gcd_results = []
    
    for comb in combs: 
        l = ls[comb[0]]
        k = ks[comb[1]]
        
        dummy = apply_merger(l,k,n)
        z_dummy = dummy[1]
        gcd_dummy = dummy[0]
        #Las siguientes líneas verifican si una solución es de la forma "vectorlike", si lo es, la rechaza
        c = 0 ##Contador para "vectorlike" 
        for i in range(1,(len(z_dummy)//2)+1):
        ##La lista con la solución está ordenada, de modo que verificar los elementos en los extremos, y así sucesivamente hacía adentro, 
        ##nos ayuda a identificar si son o no "vectorlike". Si lo son, entonces el contador "c" incrementa. Al final, el que no es
        ##"vectorlike" tiene c=0.
            if z_dummy[i-1]+z_dummy[-i] != 0:
                c += 1
        if c != 0:  #Solo guarda los que no son "vectorlike"
            z_results.append(z_dummy.tolist())
            gcd_results.append(gcd_dummy.tolist())


    return (n, gcd_results, z_results)


####################################################################
##Multiprocessing
####################################################################
    
n_s = np.arange(5,6)  ##"n" a verificar
processes = cpu_count() ##Para saber cuántos procesadores tengo 
use = 7  ##Tengo 8, entonces usaré 7
print(f"I have {processes} cores here. Using: {use}")
pool = Pool(processes=use)
results = pool.map(chiral_solution, n_s) ##Calculando los z quirales 
pool.close()
pool.join()

####################################################################              

##Diccionario para guardar los resultados
values = {"gcd": '', "z_chiral": ''}  ##Estructura interna del diccionario
dict_solutions = {str(n):copy.deepcopy(values) for n in n_s}


for line in results:
    n_x, gcd_x, z_x = line

    dict_solutions[str(n_x)]['gcd'] = gcd_x
    dict_solutions[str(n_x)]['z_chiral'] = z_x
        

print (dict_solutions['5']) 

#Dump data in a yaml file
with open('Chiral_solutions.yaml', 'w') as outfile:
    yaml.dump(dict_solutions, outfile)
                  
     

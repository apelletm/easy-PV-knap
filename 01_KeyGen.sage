############################################
## This file contains code that generates ##
## the secret and public information      ##
## of some PV-knap instance               ##
############################################

from sage.rings.polynomial.cyclotomic import cyclotomic_coeffs

def param_set(param_choice, m = 32, q = 'NULL', t = 'NULL'):
  ## param_choice can take 6 values: 'LWA1', 'LZA2', 'HPSSW1', 'HPSSW2', 'HPSSW3', 'HPSSW4' and 'customized'
  ## the first 6 choices correspond to parameters suggested in the LWA and HPSSW papers
  ## for the "customized' choice, the values of m, q and t can be chosen freely 
  ## (provided as input of the function). 
  ## If q is unspeficied, it will be chosen as the smallest prime = 1 mod m and larger than m*phi(m)
  ## If t is unspecified, it will be chosen as m/2
  if param_choice == 'LZA1': 
    m = 1024
    q = 2^16+1
    t = 256
  elif param_choice == 'LZA2':
    m = 2048
    q = 2^16+1
    t = 512
  elif param_choice == 'HPSSW1': 
    m = 433
    q = 775937
    t = 200
  elif param_choice == 'HPSSW2':
    m = 577
    q = 743177
    t = 280
  elif param_choice == 'HPSSW3':
    m = 769
    q = 1047379
    t = 386
  elif param_choice == 'HPSSW4':
    m = 1153
    q = 968521
    t = 600
  elif param_choice == 'customized':
    if q == 'NULL':
      q = m*euler_phi(m)+1
      while not is_prime(q):
        q = q+m
    if not is_prime(q) or not q%m == 1:
      print("warning, bad choice of parameters, q should be prime and =1 mod m")
    if t == 'NULL':
      t = euler_phi(m)//2
  else:
    print("\n !! Warning, undefined parameter choice !! \n")
  return (m,q,t)
  
#######################
## Some auxilliary functions related 
## to order of elements mod m
#######################

def order(x,m):
  # compute the order of x mod m (multiplicatively)
  if gcd(x,m) != 1:
    return 0
  res = 1
  x_tmp = x
  while not (x_tmp - 1)%m == 0:
    res += 1
    x_tmp = x_tmp*x
  return res
  
def all_orders(m):
  dic_orders = {}
  for i in range(m):
    o = order(i,m)
    if not o in dic_orders:
      dic_orders[o] = [i]
    else:
      dic_orders[o] += [i]
  return dic_orders
    
def max_order_bounded(dic_orders,b):
  #return the maximal order of an element mod m smaller than b
  max_o = -1
  for o in dic_orders:
    if o > max_o and o <= b:
      max_o = o
  return max_o
  
################################
## Functions generating the set
## Omega (worst case or average case)
###############################

def worst_case_Omega(m,t, size_stable_subset = 'NULL', verbose = True):
  ## Computes a set Omega that is stabilized by a cyclic subgroup
  ## of Z/mZ* of cardinality size_stable_subset
  ## (for simplicity, we only consider cyclic subgroups and not arbitrary subgroups)
  ## if size_stable_subset is unspecified, it is set to be the largest possible one smaller than t
  ## (this is the optimal case for the attack)
  dic_orders = all_orders(m)
  if size_stable_subset == 'NULL':
    size_stable_subset = max_order_bounded(dic_orders,t)
  actual_size = max_order_bounded(dic_orders, size_stable_subset)
  if actual_size != size_stable_subset:
    if verbose:
      print("There were no cyclic subsets of the requested size, the largest size below was chosen: ",actual_size)
  j0 = dic_orders[actual_size][0] # j0 has order actual_size mod m
  
  Omega = []
  for _ in range(t//actual_size):
    l = 1
    while (l < m) and ((l in Omega) or (gcd(l,m) != 1)):
      l = l+1
    Omega += [(l*j0^i)%m for i in range(actual_size)]
  
  if len(Omega) < t: ## this simply means that we discard some public information
                     ## this can only increase the hardness of the problem
    if verbose:
      print("Omega has size ",len(Omega)," instead of the allowed size ",t)
  if len(Omega) > t: ## this should not happen
    print("Warning, Omega has size ",len(Omega)," bigger than t = ",t)
  
  return Omega
  
def random_Omega(m,t, verbose = True):
  ## Computes a uniformly random set Omega of size t
  d = euler_phi(m)
  Omega = []
  if t > d:
    t = d
    if verbose:
      print("warning, t was too big, it was reduced to phi(m)")
  while len(Omega) < t:
    i = ZZ.random_element(m)
    if (gcd(i,m) == 1) and (not i in Omega):
      Omega += [i]
  return Omega
  
##################################
## Generating a PV-knap instance
##################################

def generate_PV_knap(param_choice, m = 32, q = 'NULL', t = 'NULL', Omega_sampling = 'worst-case', Omega = 'NULL', size_stable_subset = 'NULL', verbose = True):
  ## Generates a PV-knap instance
  ## param_choices, m, q and t are the same as in param_set function
  ## Omega_sampling can be set to 
  ## - 'worst-case' if one wants a worst-case choice 
  ## (size_stable_subset is as in the function worst_case_Omega)
  ## - 'random' for a random choice of Omega
  ## - 'given' if one wants to provide an already generated Omega
  ## in this case, Omega should be provided in the variable Omega
  ## as a list of size <= t of elements coprime to m (in [0,m-1])
  
  PV_knap_instance = {}
  
  ## generating the parameters
  (m,q,t) =  param_set(param_choice, m, q, t)
  d = euler_phi(m)
  Fq = FiniteField(q)
  Poly_q = Fq['x']
  Phi_q = Poly_q(cyclotomic_coeffs(m))
  
  PV_knap_instance["m"] = m
  PV_knap_instance["q"] = q
  PV_knap_instance["t"] = t
  PV_knap_instance["d"] = d
  
  ## generating the set Omega
  if Omega_sampling == 'worst-case':
    Omega = worst_case_Omega(m,t, size_stable_subset, verbose = verbose)
  elif Omega_sampling == 'random':
    if param_choice in ['HPSSW1', 'HPSSW2', 'HPSSW3', 'HPSSW4']:
      Omega = random_Omega(m,t-1, verbose = verbose) ## removing 1 evaluation point (cf explanation below)
    else:
      Omega = random_Omega(m,t, verbose = verbose)
  elif Omega_sampling == 'given':
    if Omega == 'NULL':
      print("Warning, Omega should be provided in generate_PV_knap")
      return -1
  else:
    print("Warning, invalid input in generate_PV_knap")
    return -1
    
  PV_knap_instance["Omega"] = Omega
  
  ## generating omega, a root of unity mod q
  omega = (Phi_q.roots())[-1][0]
  if Phi_q(omega) != Fq(0):
    print("Warning! omega is not a good primitive m-th root mod q")
  PV_knap_instance["omega"] = omega

  ## generating the secret key
  if param_choice in ['HPSSW1', 'HPSSW2', 'HPSSW3', 'HPSSW4']:
    ## For the parameters in HPSSW, the polynomials are over the
    ## cyclic ring Z[X]/(X^m-1) (m prime) instead of over the cyclotomic field
    ## we reduce these instance to elements over the field by forgetting
    ## the evaluation at 1 (which means that t might get reduced by 1 in the random case)
    ## The generation of the secret polynomial also happens mod X^m-1,
    ##  and we then reduce it in the field 
    e = [ZZ.random_element(-1,2) for i in range(m)]
  else:
    e = [ZZ.random_element(-1,2) for i in range(d)]  
  PV_knap_instance["secret"] = vector(ZZ,e)
  
  ## generating the public key
  poly_e = Poly_q(e)%Phi_q
  pk = {i:poly_e(omega^i) for i in Omega}
  if len(pk) != len(Omega):
    print("Warning, something went wrong in generate_PV_knap")
  PV_knap_instance["public"] = pk
  return PV_knap_instance
  
  
  
  

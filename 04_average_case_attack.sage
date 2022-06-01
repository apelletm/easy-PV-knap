#####################################
## This file contains functions    ##
## to test the efficiency of the   ##
## attack when the set Omega is    ##
## sampled uniformly at random     ##
#####################################

load("03_worst_case_attacks.sage")

######################################
## Some functions that computes
## an optimal subset of Omega that is
## stable by as many automorphisms as
## possible and is not too small
######################################

def best_subset_Omega(Omega,m,order, dic_orders):
  # Input: a set Omega, dic_orders = all_orders(m) and order is an element in dic_orders
  # Output: a subset Omega' of Omega of maximal size that is stabilized by
  # an element of order 'order'
  optimal_subset = []
  for j in dic_orders[order]:
    fixed = []
    for i in Omega:
      good = True
      for l in range(order):
        if (i*j^l)%m not in Omega:
          good = False
      if good:
        fixed += [i]
    if len(fixed) > len(optimal_subset):
      optimal_subset = fixed
  return optimal_subset
  
def best_subset_not_too_small(Omega, m, dic_orders, min_size):
  # Input: Omega and dic_orders as in best_subset_Omega, min_size an integer <= len(Omega)
  # Output: a subset Omega' of Omega of size at least min_size, which is stabilized
  # by an element of maximal possible order
  if min_size > len(Omega):
    min_size = len(Omega)
    print("Warning, the resquested size for the subset was larger than the set's size")
  best_order = 1
  best_subset = Omega
  for order in range(2,m):
    if order in dic_orders:
      optimal_subset = best_subset_Omega(Omega,m,order, dic_orders)
      if len(optimal_subset) >= min_size:
        best_order = order
        best_subset = optimal_subset
      else:
        return best_order
  return best_order


########################################
## Estimating empirically the maximal
## size of v in q*I_Omega^{-1} that
## can be used to distinguish between
## random and PV-knap
########################################


def nb_small(g,d,q):
  ## compute the number of coefficients of g mod q smaller than q/4 (in absolute value)
  nb_small = 0
  
  for i in range(d):
    if abs(center_lift(g[i],q)) <= q/4:
      nb_small += 1
  return RR(nb_small/d)
  

def distinguish_proba_fixed_f(f, m, q, nb_e, is_HPSSW = False):
  ## Use nb_small to distinguish e small from e random
  ## The probability that the distinguishing algorithm output 1
  ## is computed empirically by sampling nb_e small e and nb_e random ones
  
  PolyZ = ZZ['x']
  y = PolyZ.gen()
  Phi_m = PolyZ(cyclotomic_coeffs(m))
  d = euler_phi(m)
  output_1_small = 0
  output_1_random = 0
  for _ in range(nb_e):
    #small e
    if is_HPSSW:
      e = PolyZ(0)
      for i in range(m):
        e = PolyZ(e+ZZ.random_element(-1,2)*y^i)
      #e = ZZ['x']([ZZ.random_element(-1,2) for i in range(m)]) ## faster but memory leakage here
      e = e%Phi_m
    else:
      e = PolyZ(0)
      for i in range(d):
        e = PolyZ(e+ZZ.random_element(-1,2)*y^i)
      #e = ZZ['x']([ZZ.random_element(-1,2) for i in range(d)]) ## faster but memory leakage here
    g = (e*f)%Phi_m
    nb_small_coeff = nb_small(g,d,q)
    if nb_small_coeff > 1/2:
      output_1_small += 1 
    
    #random e
    e = PolyZ(0)
    for i in range(d):
      e = PolyZ(e+ZZ.random_element(0,q)*y^i)
    #e = PolyZ([ZZ.random_element(0,q) for i in range(d)]) ## faster but memory leakage here
    g = (e*f)%Phi_m
    nb_small_coeff = nb_small(g,d,q)
    if nb_small_coeff > 1/2:
      output_1_random += 1 
  return (RR(output_1_small/nb_e), RR(output_1_random/nb_e))
  
def Hoeffding_bound(dist, nb_samples):
  # Compute an upper boud on the probability that the theoretical mean is at distance
  # more that dist from the empirical mean when computed with nb_samples samples
  
  #Hoeffding's inequality
  # Pr(|emp_mean-actual_mean| >= t/nb_samples) <= 2 exp(-2t^2/nb_samples)
  t = RR(nb_samples*dist)
  return min(1,RR(2*exp(-2*t^2/nb_samples)))
  
def advantage_fixed_f(f, m, q, nb_e, is_HPSSW = False, error_proba = 0.01):
  ## Computes the distinguishing advantage for a given polynomial f
  ## outputs [Adv,meaningful] where Adv is a lower bound on the advantage
  ## and meaningful = True/False is a boolean saying whether the result
  ## is meaningful or not, given the error_proba allowed
  (emp_p1,emp_p2) = distinguish_proba_fixed_f(f, m, q, nb_e, is_HPSSW)
  dist = abs(emp_p1-emp_p2)*9/20 # assume that the actual mean are at distance at most 9/20 of the distance emp_p1-emp_p2
                                 # so that the distance between the theoretical means is >= 1/10*(emp_p1 - emp-p2)
  meaningful = (Hoeffding_bound(dist, nb_e) <= error_proba)
  return [abs(emp_p1-emp_p2)/10, meaningful]
  
def max_size_non_negl_advantage(m, q, min_f = 1, max_nb_e = 5000, error_proba = 0.01, is_HPSSW = False):
  ## Try random f of increasing size until one does not manage to distinguish anymore between 
  ## small secrets and random secrets
  ## it outputs the infinity norm of the largest f that enables distinguishing
  size_f = min_f
  nb_e = 100
  stop = False
  while (not stop) and size_f < q:
    f = ZZ['x']([ZZ.random_element(-size_f, size_f+1) for i in range(euler_phi(m))])
    [Adv, meaningful] = advantage_fixed_f(f,m,q,nb_e,is_HPSSW, error_proba)
    if meaningful:
      size_f += q//100
    else:
      if nb_e >= max_nb_e:
        stop = True
      else:
        nb_e = min(max_nb_e, nb_e*5)
  
  # reducing the gap q//100 to q//400 by dichotomy
  min_best_size = max(min_f,size_f - q//100)
  max_best_size = size_f
  for i in range(2):
    size_f = (max_best_size+min_best_size)//2
    f = ZZ['x']([ZZ.random_element(-size_f, size_f+1) for i in range(euler_phi(m))])
    [Adv, meaningful] = advantage_fixed_f(f,m,q,nb_e,is_HPSSW, error_proba)
    if meaningful:
      min_best_size = size_f
    else:
      max_best_size = size_f
  
  return min_best_size
        
#######################################
## Estimating the minimal size Omega'
## that can be used to distinguish
## the BDD input from uniform
#######################################

def minimal_size_Omega(m, q, min_f = 1, max_nb_e = 5000, error_proba = 0.01, is_HPSSW = False, verbose = True):
  
  # Computing the shortest infinity norm that allows distinguishing
  if verbose:
    print("\nComputing the empirical length of shortest vector needed for the attack (this can take a few minutes)...")
  max_allowed_norm = max_size_non_negl_advantage(m, q, min_f, max_nb_e, error_proba, is_HPSSW)
  if verbose:
    print("...done")
    print("\nMax infinity norm of short vector in q*I_Omega^{-1} required:", max_allowed_norm)
  
  # We estimate the infinity norm of the minimal vector in qI_Omega'^{-1}
  # by N(q*I_Omega'^{-1})^{1/d} = q^((d-len(Omega'))/d)
  min_size_Omega = euler_phi(m)*(1-log(max_allowed_norm)/log(q))
  if verbose:
    print("\nMinimal size of Omega' required for distinguishing:", ZZ(ceil(min_size_Omega)))
  return ZZ(ceil(min_size_Omega))
  
  
######################################
## Empirical probability that a subset
##  of Omega of a given size if fixed 
## by a subgroup of given order
######################################
  
def proba_best_subgroup_Omega(param_choice, m = 32, q = 'NULL', t = 'NULL', min_f = 1, max_nb_e = 5000, nb_Omega = 1000, error_proba = 0.01, verbose = True):

  # generating parameters
  (m,q,t) = param_set(param_choice, m, q, t)
  t = min(t,euler_phi(m))
  dic_orders = all_orders(m)
  is_HPSSW = False
  if param_choice in ['HPSSW1', 'HPSSW2', 'HPSSW3', 'HPSSW4']:
    is_HPSSW = True
    
  # compute minimal size Omega'
  min_size_subset = minimal_size_Omega(m, q, min_f, max_nb_e, error_proba, is_HPSSW, verbose)
  
  # Generate nb_Omega random sets Omega and compute the maximal order that stabilize a subset of size at least min_size_subset
  if verbose:
    print("Computing the empirical probability that Omega contains a large enough subset of given order")
  probabilities = {}
  for _ in range(nb_Omega):
    Omega = random_Omega(m,t, verbose = True)
    res = best_subset_not_too_small(Omega, m, dic_orders, min_size_subset)
    if not res in probabilities:
      probabilities[res] = 1
    else:
      probabilities[res] += 1
    
  for order in probabilities:
    probabilities[order] = RR(probabilities[order]/nb_Omega)
    
  return probabilities
  
##########################################
## Actually running the distinguishing
## attack
##########################################

def distinguishing_attack_random(param_choice, m = 32, q = 'NULL', t = 'NULL', nb_Omega = 3000, verbose = True, reduced_dim = 2, block_size = 55, save_in_file = True, max_hours = 10, max_nb_e = 500000):
  # run the distinguishing attack
  # nb_Omega corresponds to the number of times we are allowed to re-sample Omega before 
  # reduced_dim is the reduction in the dimension that we are aiming at
  
  # Parameters
  (m,q,t) = param_set(param_choice, m, q, t)
  t = min(t,euler_phi(m))
  omega = generate_PV_knap(param_choice,m,q,t, verbose = False)["omega"]
  is_HPSSW = False
  if param_choice in ['HPSSW1', 'HPSSW2', 'HPSSW3', 'HPSSW4']:
    is_HPSSW = True
  dic_orders = all_orders(m)
  
  # Generating a good Omega
  if param_choice in ["LZA1", "LZA2"] and reduced_dim == 2:
    (best_Omega, best_subset) = best_Omega_LZA(m,t, nb_Omega, verbose = verbose) # imposing a subfield with half the coefficients = 0
  else:
    best_subset = []
    best_Omega = []
    for _ in range(nb_Omega):
      Omega = random_Omega(m,t, verbose = True)
      subset = best_subset_Omega(Omega,m,reduced_dim, dic_orders)
      if len(subset) > len(best_subset):
        best_Omega = Omega
        best_subset = subset
        
  Omega = best_Omega # to prevent stupid mistakes
  if len(best_subset) == 0:
    print("Could not find a subset reducing the dimension as required\nAborted")
    return -1
  if verbose:
    expected_shortest = ZZ(round(q^((euler_phi(m)-len(best_subset))/euler_phi(m))))
    print("\nThe best Omega found has a stable subset of size", len(best_subset))
    print("From this, we expect a shortest vector of infinity norm", expected_shortest,"\n")
  if save_in_file:
    output_name = "data/params_"+str(euler_phi(m))+"_BKZ_"+str(block_size)+"_time_"+str(max_hours)+"h.sage"
    f = open(output_name, 'w')
    f.write("param_choice = '"+param_choice+"'\nm = "+str(m)+"\nq = "+str(q)+"\nt = "+str(t)+"\nnb_Omega = "+str(nb_Omega)+"\nOmega = "+str(best_Omega)+"\nsubset = "+str(best_subset))
    f.close()
    
  # Computing a short vector in q*I_(best_subset)^{-1} (included in q*I_{Omega}^{-1})
  v = small_vector_I_Omega_inverse(best_subset,omega,q,m,fast = True, k = 2, block_size = block_size, auto_abort = True, max_hours = max_hours, float_type = 'NULL', precision = 60, save_result_in_file = save_in_file, verbose = verbose)
  if verbose:
    print("\nShortest vector computed has infinity norm:", v.norm(infinity),"\n(to be compared with the shortest infinity norm needed)\n")
  
  # double checking that v is indeed in I_{Omega}^(-1)
  # by checking that the attack indeed produces v*e mod q, where e is the small PV-knap secret
  old_values = 'NULL'
  poly_v = (FiniteField(q)['y'])(list(v))
  
  for _ in range(2):
    PV_knap_instance = generate_PV_knap(param_choice,m,q,t, Omega_sampling = 'given', Omega = best_Omega, verbose = False)
    (target, old_values) = create_target(PV_knap_instance["public"],omega, q, m, old_values)
    r_pub = centered_product(m,q,target,poly_v)
    poly_e = (FiniteField(q)['y'])(list(PV_knap_instance["secret"]))
    r_s = centered_product(m,q,poly_e,poly_v)
    if r_pub != r_s:
      print(r_s)
      print(r_pub)
      print("\n\n!! Warning !! The short vector computed is not in I_{Omega}^(-1)... Aborted")
      return -1
  
  # Estimate the success probability of the attack
  # In order to run tests faster, we only compute small e and random e, and try to distinguish e*v mod q
  # (the above verification ensures that this is indeed what the attack would obtain on input the public key pk)
  nb_e = 100
  meaningful = False
  f = (ZZ['x'])(list(v))
  while nb_e < max_nb_e and not meaningful:
    nb_e = min(10*nb_e, max_nb_e)
    (Adv, meaningful) = advantage_fixed_f(f, m, q, nb_e, is_HPSSW, error_proba = 0.01)
  if meaningful:
    print("\n\n=== Average-case distinguishing attack: Success ===")
    print("After chosing among ", nb_Omega, " public keys, it was able to distinguish PV-knap instance from random with probability at least ", RR(Adv))
    return 1
  else:
    print("\n\n=== Average-case distinguishing attack: Failure ===")
    print("You can try running it again with larger block size or more time")
    if save_in_file:
      print("(The reduced basis of the ideal computed so far and the parameters have been writen in two files.")
      print("You may want to use the function distinguishing_attack_random_a_bit_more)")
    return 0
  
  
  
  
def distinguishing_attack_random_a_bit_more(param_input_file, mat_input_file, block_size = 55, save_in_file = True, max_hours = 10, max_nb_e = 500000, auto_abort = True, verbose = True):

  # recovering reduced basis and parameters
  M = IntegerMatrix.from_file(mat_input_file)
  load(param_input_file)
  if save_in_file:
    output_file_name = mat_input_file+"_"+str(block_size)+"_"+str(max_hours)+"h"
  is_HPSSW = False
  if param_choice in ['HPSSW1', 'HPSSW2', 'HPSSW3', 'HPSSW4']:
    is_HPSSW = True
  
  # Reducing the basis even more
  M_red = reducing_basis(M, block_size = block_size, auto_abort = auto_abort, max_hours = max_hours, save_result_in_file = save_in_file, output_file_name = output_file_name, verbose = verbose)
  
  # short vector
  v = best_vector(M_red)
  if verbose:
    d = euler_phi(m)
    expected_shortest = RR(q^((d-len(subset))/d)*sqrt(d/(2*pi*e))) #expected length of a shortest vector in M_small
    print("\nExpected length of shortest vector: ", expected_shortest)
    print("Actual length of computed vector: ", RR(v.norm(2)))
    
  # double checking that v is indeed in I_{Omega}^(-1)
  # by checking that the attack indeed produces v*e mod q, where e is the small PV-knap secret
  old_values = 'NULL'
  poly_v = (FiniteField(q)['y'])(list(v))
  for _ in range(2):
    PV_knap_instance = generate_PV_knap(param_choice,m,q,t, Omega_sampling = 'given', Omega = Omega, verbose = False)
    (target, old_values) = create_target(PV_knap_instance["public"],PV_knap_instance["omega"], q, m, old_values)
    r_pub = centered_product(m,q,target,poly_v)
    poly_e = (FiniteField(q)['y'])(list(PV_knap_instance["secret"]))
    r_s = centered_product(m,q,poly_e,poly_v)
    if r_pub != r_s:
      print(r_s)
      print(r_pub)
      print("\n\n!! Warning !! The short vector computed is not in I_{Omega}^(-1)... Aborted")
      return -1
  
  # Estimate the success probability of the attack
  # In order to run tests faster, we only compute small e and random e, and try to distinguish e*v mod q
  # (the above verification ensures that this is indeed what the attack would obtain on input the public key pk)
  nb_e = 100
  meaningful = False
  f = (ZZ['x'])(list(v))
  while nb_e < max_nb_e and not meaningful:
    nb_e = min(10*nb_e, max_nb_e)
    (Adv, meaningful) = advantage_fixed_f(f, m, q, nb_e, is_HPSSW, error_proba = 0.01)
  if meaningful:
    print("\n\n=== Average-case distinguishing attack: Success ===")
    print("After chosing among ", nb_Omega, " public keys, it was able to distinguish PV-knap instance from random with probability at least ", RR(Adv))
    return 1
  else:
    print("\n\n=== Average-case distinguishing attack: Failure ===")
    print("You can try running it again with larger block size or more time")
    if save_in_file:
      print("(The reduced basis of the ideal computed so far has been writen in a file,")
      print("you may want to use the function distinguishing_attack_random_a_bit_more)")
    return 0


#########################################
## Special functions for power-of-two
## cyclotomic fields, to try to optimize
## the attack: consider only the degree
## 2 subfield with half the coefficients
## that are 0 (so that we gain a factor
## two on the size of the vectors, this
## seems to help in practice)
##########################################

def best_Omega_LZA(m,t, nb_Omega = 3000, verbose = True):
  ## Compute the best Omega among nb_Omega random ones
  ## specific to power-of-two cyclotomic fields
  ## only consider the subgroup of order 2 corresponding
  ## to the cyclotomic field of degree d/2
  if verbose:
    print("Using a special subgroup of order 2 for the LZA parameter sets")
  best_Omega = []
  best_subset = []
  
  for _ in range(nb_Omega):
    Omega = random_Omega(m,t, verbose = True)
    
    j = euler_phi(m)+1
    if order(j,m) != 2:
      print("!! warning !! The element d+1 has not order 2 modulo m")

    fixed = []
    for i in Omega:
      if (i*j)%m in Omega:
        fixed += [i]
        
    if len(fixed) > len(best_subset):
      best_subset = fixed
      best_Omega = Omega
  return (best_Omega, best_subset)




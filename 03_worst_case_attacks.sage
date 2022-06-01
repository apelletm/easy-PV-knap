#####################################
## This file contains functions    ##
## to test the attack in the worst ##
## case scenario, i.e., when       ##
## Omega is specifically chosen to ##
## be stable by a large subgroup   ##
#####################################

load("02_short_vector_inverse.sage")

##############################################
## Transform the PV-knap instance into a
## CVP instance in the ideal I_Omega
## (computes the target of this CVP instance)
##############################################

def create_target(pk,omega, q, m, old_values = 'NULL'):
  ## create the linear system to solve
  
  # old_values = [coprime_m,M_inv] contains potentially already computed values
  # of coprime_m and M^(-1) (they depend only on m and omega and can be reused accross tests)
  # if those values where never computed before, they are set to NULL
  Fq = FiniteField(q)
  if old_values == 'NULL': 
    # list of integers coprime to m
    coprime_m = []
    for i in range(m):
      if gcd(i,m) == 1:
        coprime_m += [i]
    # Vandermonde matrix M corresponding to evaluation at Omega^j (j coprime to m)
    list_M = []
    for i in coprime_m:
      list_M += [[Fq(omega)^(i*k) for k in range(euler_phi(m))]]
    M = Matrix(Fq,list_M)
    M_inv = M^(-1)
    old_values = [coprime_m, M_inv]
   
  else: 
    coprime_m = old_values[0]
    M_inv = old_values[1]  
    
  # Now the part that depends on pk
  list_b = []
  for i in coprime_m:
    if i in pk:
    ## evaluation of secret polynomial at omega^i (i in Omega)
      list_b += [Fq(pk[i])] 
    else:
    ## imposing arbitrary values to evaluations in Omega_c    
      list_b += [Fq(0)]
      
  b = vector(Fq,list_b)
  z = vector(M_inv*b)
  t = (Fq['y'])(list(z))
  for i in pk:
    if t(omega^i) != pk[i]:
      print("Warning! the lift polynomial t is not good")
  return (t, old_values)
  

def centered_product(m,q,t,s):
  # Given t the target vector which is a BDD in I_Omega
  # and s a short vector in q*I_Omega^{-1}
  # computes the product t*s with centered coefficients (mod q)
  
  r = (t*s)%(FiniteField(q)['y'])(cyclotomic_coeffs(m))
  
  # if v was sufficiently small, r mod q is equal to f*e over ZZ 
  # (with e the secret small polynomial)
  # even if it is not, it might still be distinguishable from
  # uniformly random modulo q
  lift_r = ZZ['x'](list(poly_to_vec_center(r,euler_phi(m),q)))
  return lift_r
  
###########################
## Key recovering attack
###########################
  
def test_successful_key_recovery(param_choice, guess_e, e, m):
  # Test whether a guessed polynomial guess_e actually corresponds to the secret e of the PV_knap instance
  # This is almost immediate except for the HPSSW parameters, which need more care because of the difference between
  # the solution in the cyclotomic field and the actual secret in the ring Z[X]/(X^m-1)
  if param_choice in ['HPSSW1', 'HPSSW2', 'HPSSW3', 'HPSSW4']:
    # we need to recontruct the ternary secret e from guess_e = e mod Phi_m(X) (if the guess is correct)
    coeff_e = vector(QQ,m)
    for i in range(len(list(guess_e))):
      coeff_e[i] = guess_e[i]
    if -2 in coeff_e:
      for i in range(len(coeff_e)):
        coeff_e[i] += 1
    elif 2 in coeff_e:
      for i in range(len(coeff_e)):
        coeff_e[i] += -1
        
  else:
    coeff_e = vector(QQ,euler_phi(m))
    for i in range(len(list(guess_e))):
      coeff_e[i] = guess_e[i]
  
  return (coeff_e == e)
      
      

def key_recovering_attack(param_choice, m = 32, q = 'NULL', t = 'NULL', size_stable_subset = 'NULL', nb_tests = 20, fast =  True, verbose = True):
  # Parameters are as in the generate_PV_knap function
  # the extra parameter nb_tests determines how PV-knap instance we generate
  
  successful_key_recovery = 0
  
  # generate Omega and compute a short vector in I_Omega^{-1}
  if verbose:
    print("\nComputing parameters and the short vector in the inverse ideal...")
  PV_knap_instance = generate_PV_knap(param_choice, m, q, t, Omega_sampling = 'worst-case', Omega = 'NULL', size_stable_subset = size_stable_subset, verbose = verbose)
  (m,q,Omega,omega) = (PV_knap_instance["m"],PV_knap_instance["q"],PV_knap_instance["Omega"],PV_knap_instance["omega"])
  if verbose:
    print("\nParameters: m = ", m, ", q = ", q, ", t = ", PV_knap_instance["t"])
    
  ## Computing a short vector in I_Omega^{-1}
  v = small_vector_I_Omega_inverse(Omega,omega,q,m,fast = fast, k = 2, block_size = 20, auto_abort = True, max_hours = 1, float_type = 'NULL', precision = 60, save_result_in_file = False, verbose = verbose)
  short_vec = (FiniteField(q)['y'])(list(v))
  lift_short_vec = ZZ['x'](list(v))
  
  ## Definitions
  Phi_m = ZZ['x'](cyclotomic_coeffs(m))
  K = QuotientRing(QQ[x], Phi_m)
  old_values = 'NULL'
  
  if verbose:
    print("...done")
    
  ## Testing key recovering for different PV-knap instances (with the same worst-case Omega)
  if verbose:
    print("\n\nTesting the key recovery success ", nb_tests, " times...")
  for _ in range(nb_tests):
    PV_knap_instance = generate_PV_knap(param_choice, m, q, t, Omega_sampling = 'given', Omega = Omega, verbose = False)
    
    (target, old_values) = create_target(PV_knap_instance["public"],omega, q, m, old_values)
    r = centered_product(m,q,target,short_vec)
    guess_e = K(r)/K(lift_short_vec)
    if test_successful_key_recovery(param_choice, guess_e, PV_knap_instance["secret"], m):
      successful_key_recovery += 1
      
  if verbose:
    print("...done")
  print("\n\n===== Worst-case Key recovery attack =====")
  print("successful with probability ", RR(successful_key_recovery/nb_tests),"\n")
  
###########################
## Distinguishing attack
###########################

def test_small(r,q,m):
  ## test whether the coefficients of r are all smaller than q*0.45
  ## (in absolute value) or not
  for ri in r:
    if abs(center_lift(ri,q)) > q*0.45:
      return False
  return True

def distinguishing_attack(param_choice, m = 32, q = 'NULL', t = 'NULL', size_stable_subset = 'NULL', nb_tests = 20, fast = True, verbose = True):
  # Parameters are as in the generate_PV_knap function
  # the extra parameter nb_tests determines how PV-knap instance we generate
  
  successful_distinguishing = 0
  
  # generate Omega and compute a short vector in I_Omega^{-1}
  if verbose:
    print("\nComputing parameters and the short vector in the inverse ideal...")
  PV_knap_instance = generate_PV_knap(param_choice, m, q, t, Omega_sampling = 'worst-case', Omega = 'NULL', size_stable_subset = size_stable_subset, verbose = verbose)
  (m,q,Omega,omega) = (PV_knap_instance["m"],PV_knap_instance["q"],PV_knap_instance["Omega"],PV_knap_instance["omega"])
  if verbose:
    print("\nParameters: m = ", m, ", q = ", q, ", t = ", PV_knap_instance["t"])
    
  ## Computing a short vector in I_Omega^{-1}
  v = small_vector_I_Omega_inverse(Omega,omega,q,m,fast = fast, k = 2, block_size = 20, auto_abort = True, max_hours = 1, float_type = 'NULL', precision = 60, save_result_in_file = False, verbose = verbose)
  short_vec = (FiniteField(q)['y'])(list(v))
  lift_short_vec = ZZ['x'](list(poly_to_vec_center(short_vec,euler_phi(m),q)))
  
  ## Definitions
  old_values = 'NULL'
  
  if verbose:
    print("...done")
    
  ## Testing key recovering for different PV-knap instances (with the same worst-case Omega)
  if verbose:
    print("\n\nTesting the distinguishing success ", nb_tests, " times...")
  for _ in range(nb_tests):
  
    # Sampling a coin to compute challenge
    coin = ZZ.random_element(0,2)
    if coin:
      PV_knap_instance = generate_PV_knap(param_choice, m, q, t, Omega_sampling = 'given', Omega = Omega, verbose = False)
      pk = PV_knap_instance["public"]
    else:
      pk = {}
      for i in Omega:
        pk[i] = ZZ.random_element(q)
        
    # Distinguishing small vs large
    (target, old_values) = create_target(pk,omega, q, m, old_values)
    r = centered_product(m,q,target,short_vec)
    if test_small(r,q,m):
      guess_coin = 1
    else:
      guess_coin = 0
    
    if coin == guess_coin:
      successful_distinguishing += 1
      
  if verbose:
    print("...done")
  print("\n\n===== Worst-case distinguishing attack =====")
  print("successful with probability ", RR(successful_distinguishing/nb_tests),"\n")




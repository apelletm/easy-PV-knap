#########################################
## This file contains the code used to ##
## run the experiments reported on in  ##
## the article                         ##
#########################################

reset()
load("04_average_case_attack.sage")
global_verbose = False
RRshort = RealField(15)

###########################################
## Experiments regarding worst-case
## attacks
###########################################

## LZA1
print("\n\n========================================")
print("== Testing worst case attacks on LZA1 ==")
print("========================================")
set_random_seed(42)
key_recovering_attack("LZA1", verbose = global_verbose) # 5 seconds
set_random_seed(42)
distinguishing_attack("LZA1",verbose = global_verbose) # 4 seconds


## LZA2
print("\n\n========================================")
print("== Testing worst case attacks on LZA2 ==")
print("========================================")
set_random_seed(42)
key_recovering_attack("LZA2", verbose = global_verbose) # 25 seconds
#not taking the best possible Omega, but one that reduces the dimension "only" by 16
set_random_seed(42)
key_recovering_attack("LZA2", verbose = global_verbose, size_stable_subset = 16) # 2 minutes
print("(Omega stable by a subgroup of size 16 only)")
set_random_seed(42)
distinguishing_attack("LZA2",verbose = global_verbose) # 18 seconds


## HPSSW1
print("\n\n==========================================")
print("== Testing worst case attacks on HPSSW1 ==")
print("==========================================")
set_random_seed(42)
key_recovering_attack("HPSSW1", verbose = global_verbose) # 15 seconds
set_random_seed(42)
distinguishing_attack("HPSSW1",verbose = global_verbose) # 4 seconds

## HPSSW2
print("\n\n==========================================")
print("== Testing worst case attacks on HPSSW2 ==")
print("==========================================")
set_random_seed(42)
key_recovering_attack("HPSSW2", verbose = global_verbose) # 33 seconds
set_random_seed(42)
distinguishing_attack("HPSSW2",verbose = global_verbose) # 7 seconds

## HPSSW3
print("\n\n==========================================")
print("== Testing worst case attacks on HPSSW3 ==")
print("==========================================")
set_random_seed(42)
key_recovering_attack("HPSSW3", verbose = global_verbose) # 1 minute
set_random_seed(42)
distinguishing_attack("HPSSW3",verbose = global_verbose) # 12 seconds

## HPSSW4
print("\n\n==========================================")
print("== Testing worst case attacks on HPSSW4 ==")
print("==========================================")
set_random_seed(42)
key_recovering_attack("HPSSW4", verbose = global_verbose) # 2'30 minutes
#not taking the best possible Omega, but one that reduces the dimension "only" by 16
set_random_seed(42)
key_recovering_attack("HPSSW4", verbose = global_verbose, size_stable_subset = 16) # 4 minutes
print("(Omega stable by a subgroup of size 16 only)")
set_random_seed(42)
distinguishing_attack("HPSSW4",verbose = global_verbose) # 30 seconds

###########################################
## Experiments regarding average-case
## attacks
###########################################

######################################
## Computing the reduced_dimension
######################################

## LZA1
print("\n\n=================================================")
print("== Estimating cost average case attack on LZA1 ==")
print("=================================================")
set_random_seed(42)
probabilities = proba_best_subgroup_Omega("LZA1", verbose = global_verbose) # 2'30 minutes
for o in probabilities:
  print("Probability to reduce the degree by a factor ",o,": ",RRshort(probabilities[o]))
  
## LZA2
print("\n\n=================================================")
print("== Estimating cost average case attack on LZA2 ==")
print("=================================================")
print("\n(This one is quite long, roughly 3h)")
set_random_seed(42)
probabilities = proba_best_subgroup_Omega("LZA2", verbose = global_verbose, nb_Omega = 50000, max_nb_e = 100000) # 3h
for o in probabilities:
  print("Probability to reduce the degree by a factor ",o,": ",RRshort(probabilities[o]))
  
## HPSSW1
print("\n\n===================================================")
print("== Estimating cost average case attack on HPSSW1 ==")
print("===================================================")
set_random_seed(42)
probabilities = proba_best_subgroup_Omega("HPSSW1", verbose = global_verbose, nb_Omega = 10000) # 2'10 minutes
for o in probabilities:
  print("Probability to reduce the degree by a factor ",o,": ",RRshort(probabilities[o]))
  
## HPSSW2
print("\n\n===================================================")
print("== Estimating cost average case attack on HPSSW2 ==")
print("===================================================")
set_random_seed(42)
probabilities = proba_best_subgroup_Omega("HPSSW2", verbose = global_verbose, nb_Omega = 10000) # 3'50 minutes
for o in probabilities:
  print("Probability to reduce the degree by a factor ",o,": ",RRshort(probabilities[o]))
  
## HPSSW3
print("\n\n===================================================")
print("== Estimating cost average case attack on HPSSW3 ==")
print("===================================================")
set_random_seed(42)
probabilities = proba_best_subgroup_Omega("HPSSW3", verbose = global_verbose) # 4'40 minutes
for o in probabilities:
  print("Probability to reduce the degree by a factor ",o,": ",RRshort(probabilities[o]))
  
## HPSSW4
print("\n\n===================================================")
print("== Estimating cost average case attack on HPSSW4 ==")
print("===================================================")
set_random_seed(42)
probabilities = proba_best_subgroup_Omega("HPSSW4", verbose = global_verbose, nb_Omega = 10000) # 9'30 minutes
for o in probabilities:
  print("Probability to reduce the degree by a factor ",o,": ",RRshort(probabilities[o]))
  
  
###################################################
## Running the average case attack on LZA1
###################################################

print("\n\n=========================================")
print("== Running average case attack on LZA1 ==")
print("=========================================")

set_random_seed(42)
print("\n Starting the attack with block-size 40\n")
b40 = distinguishing_attack_random("LZA1", nb_Omega = 3000, block_size = 40, max_hours = 10, verbose = True, max_nb_e = 1000) # 2h

if b40 == 0:
  print("\n The attack did not work, let us increase the block-size to 45\n")
  b45 = distinguishing_attack_random_a_bit_more("data/params_512_BKZ_40_time_10h.sage", "data/mat_512_BKZ_40_time_10h", block_size = 45, max_hours = 10, verbose = True, max_nb_e = 1000) # 1h30
  
  if b45 == 0:
    print("\n The attack did not work, let us increase the block-size to 50\n")
    b50 = distinguishing_attack_random_a_bit_more("data/params_512_BKZ_40_time_10h.sage", "data/mat_512_BKZ_40_time_10h_45_10h", block_size = 50, max_hours = 10, verbose = True, max_nb_e = 1000000) # 7h30
    
    if b50 == 0:
      print("\n The attack did not work, let us increase the block-size to 55\n")
      b55 = distinguishing_attack_random_a_bit_more("data/params_512_BKZ_40_time_10h.sage", "data/mat_512_BKZ_40_time_10h_45_10h_50_10h", block_size = 55, max_hours = 10, verbose = True, max_nb_e = 10000000) # 8h 


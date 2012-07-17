# TODO: Add comment
# 
# Author: Sophie Baillargeon
###############################################################################

# À faire : test d'un cas avec bornes initiales qui ne bougent pas.
# trouvé le 29 juin 2012

test <- strata.LH(x = Sweden$P75, n = 80, Ls = 2, alloc = c(0.5, 0, 0), takeall = 0, initbh = 51.8073353879545, 
    algo = "Kozak", algo.control = list(minsol = 2, maxiter = 1e+06, minNh = 2, maxstep = 3, 
        maxstill = 100, method = "original", rep = 1L, idopti = "nh"))

# Attention : on ne peut pas tester si les bornes finales produites sont égales aux bornes initiales
# car même si l'algo ne bouge pas, les bornes finales sont une version modifée des bornes
# initiales : on travaille avec les positions sur le vecteur des valuers possibles et en sortie
# on prend la val moyenne entre deux valeurs différentes consécutives.


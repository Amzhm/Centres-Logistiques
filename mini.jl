using JuMP
using CPLEX
using CPUTime
using Random

include("IO_Placement.jl")

############################################******PLNE******#####################################################
function PLNE_P_median(filename,P)
    n=0
    tabX=Float64[]
    tabY=Float64[]
    f= Float64[]

    # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)

     # Calcul de la matrice des distances
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]

    m = Model(CPLEX.Optimizer)
    @variable(m, x[1:n,1:n], Bin)

    @objective(m, Min, sum(d[i,j]*x[i,j] for i in 1:n, j in 1:n))
    @constraint(m,sum(x[j,j] for j in 1:n)<=P)
    @constraint(m, [i=1:n], sum(x[i,j] for j in 1:n) == 1)
    @constraint(m, [i=1:n, j=1:n], x[i,j] <= x[j,j])

    optimize!(m)
    status = termination_status(m)

    if status == JuMP.MathOptInterface.OPTIMAL
        println("Valeur optimale = ", objective_value(m))
    else
        println("Problème non optimalement résolu.")
        return
    end

    S = zeros(Int, n)
    centres=Int[]
    for j in 1:n
        if (value(x[j,j])) == 1.0
            S[j] = 1
            #println("\t x[",j,",",j,"] = ", value(x[j,j]))
            push!(centres,j)
        end
    end

    #Dessine_Solution_Placement(filename, n, tabX, tabY, S)
    return centres,S

end

function PLNE_P_median_with_f(filename,P)
    n=0
    tabX=Float64[]
    tabY=Float64[]
    f= Float64[]

    # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
    
    # Calcul de la matrice des distances
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]

    m = Model(CPLEX.Optimizer)
    @variable(m, x[1:n,1:n], Bin)

    @objective(m, Min, sum(d[i,j]*x[i,j] for i in 1:n, j in 1:n)+sum(f[j]*x[j,j] for j in 1:n))
    @constraint(m,sum(x[j,j] for j in 1:n)<=P)
    @constraint(m, [i=1:n], sum(x[i,j] for j in 1:n) == 1)
    @constraint(m, [i=1:n, j=1:n], x[i,j] <= x[j,j])

    optimize!(m)
    status = termination_status(m)

    if status == JuMP.MathOptInterface.OPTIMAL
        println("Valeur optimale = ", objective_value(m))
    else
        println("Problème non optimalement résolu.")
        return
    end

    total_distance = 0.0
    total_installation = 0.0
    for i in 1:n
        for j in 1:n
            total_distance += d[i, j] * value(x[i, j])
        end
    end
    for j in 1:n
        total_installation += f[j] * value(x[j, j])
    end

    S = zeros(Int, n)
    centres=Int[]
    for j in 1:n
        if (value(x[j,j])) == 1.0
            S[j] = 1
            #println("\t x[",j,",",j,"] = ", value(x[j,j]))
            push!(centres,j)
        end
    end
    #Dessine_Solution_Placement(filename, n, tabX, tabY, S)
    return centres,S

end

function PLNE_P_centre(filename,P)
    n=0
    tabX=Float64[]
    tabY=Float64[]
    f= Float64[]

    # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
    
    # Calcul de la matrice des distances
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]

    m = Model(CPLEX.Optimizer)
    @variable(m, x[1:n,1:n], Bin)
    @variable(m, z>=0)

    @objective(m, Min, z)
    @constraint(m, sum(x[j,j] for j in 1:n) <= P)
    @constraint(m, [i=1:n], sum(x[i,j] for j in 1:n) == 1)
    @constraint(m, [i=1:n], sum(d[i,j] * x[i,j] for j in 1:n) <= z)
    @constraint(m, [i=1:n, j=1:n], x[i,j] <= x[j,j])

    optimize!(m)
    status = termination_status(m)

    if status == JuMP.MathOptInterface.OPTIMAL
        println("Valeur optimale = ", objective_value(m))
    else
        println("Problème non optimalement résolu.")
        return
    end

    S = zeros(Int, n)
    centres=Int[]
    for j in 1:n
        if (value(x[j,j])) == 1.0
            S[j] = 1
            #println("\t x[",j,",",j,"] = ", value(x[j,j]))
            push!(centres,j)
        end
    end
    #Dessine_Solution_Placement(filename, n, tabX, tabY, S)
    return centres,S
end

############################################******Glouton******#####################################################
function P_median_greedy_random(filename, P)
    n = 0
    tabX = Float64[]
    tabY = Float64[]
    f = Float64[]
    # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
     # Calcul de la matrice des distances
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]
    # Tirage aléatoire de P centres différents
    centres = randperm(n)[1:P]
    S = zeros(Int, n)
    for j in centres
        S[j] = 1
    end
    return centres,S
end

function P_centre_greedy_random(filename, P; l=3)
    n = 0
    tabX = Float64[]
    tabY = Float64[]
    f = Float64[]
    max_min_dist=0
    # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
     # Calcul de la matrice des distances
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]

    # Initialisation avec un centre aléatoire
    centre_initial = rand(1:n)
    S = [centre_initial]

    while length(S) < P
        dist_min = Dict{Int, Float64}()

        for j in setdiff(1:n, S)
            dist_min[j] = minimum(d[j, k] for k in S)
        end

        # Sélection des l villes les plus mal desservies
        pires = sort(collect(keys(dist_min)), by = j -> -dist_min[j])
        RCL = pires[1:min(l, length(pires))]

        # Choix aléatoire dans la liste restreinte
        choisi = rand(RCL)
        push!(S, choisi)
    end
    solution = zeros(Int, n)
    for j in S
        solution[j] = 1
    end
    return S,solution
end

function P_median_greedy(filename, P)
    tabX = Float64[]
    tabY = Float64[]
    f=Float64[]
    # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
    
    # Calcul de la matrice des distances
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]
    
    S = Int[]
    meilleur_cout = 1e10
    ordre_villes = randperm(n)  # Permutation aléatoire des indices de villes

    while length(S) < P
        meilleur = -1
        meilleur_cout_local = 1e10

        #Recherche de la ville qui améliore la solution
        for j in ordre_villes
            if j in S
                continue
            end
            tempS = push!(copy(S), j)
            #Calcul du cout de la solution
            total = sum(minimum(d[i, k] for k in tempS) for i in 1:n)
            if total < meilleur_cout_local
                meilleur_cout_local = total
                meilleur = j
            end
        end
        if meilleur != -1
            push!(S, meilleur)
            meilleur_cout = meilleur_cout_local
        else
            break
        end
    end

    solution = zeros(Int, n)
    for j in S
        solution[j] = 1
    end

    return S, solution
end

function P_centre_choix_villes_eloignées(filename, P)
    n = 0
    tabX = Float64[]
    tabY = Float64[]
    f = Float64[]
    max_min_dist = -1.0

    # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
    
    # Calcul de la matrice des distances
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]
    
    # choisir les 2 villes les plus éloignées
    max_dist = -1.0
    a, b = 1, 2
    for i in 1:n
        for j in i+1:n
            if d[i,j] > max_dist
                max_dist = d[i,j]
                a, b = i, j
            end
        end
    end
    S = [a, b]
    while length(S) < P
        # Trouver la ville la plus éloignée des centres actuels
        max_min_dist = -1.0
        meilleur = -1
        #Calcul de la distance Max entre une ville est son centre
        for j in setdiff(1:n, S)
            dist_min = minimum(d[j, k] for k in S)
            if dist_min > max_min_dist
                max_min_dist = dist_min
                meilleur = j
            end
        end
        push!(S, meilleur)
    end

    solution = zeros(Int, n)
    for j in S
        solution[j] = 1
    end
    return S,solution
end



############################################******METAHEURISTIQUE******#####################################################
function P_median_greedy_meta(filename, P,maxIter)
    n = 0
    tabX = Float64[]
    tabY = Float64[]
    f = Float64[]
    nbIter = 0
    nb_stable=0
    min_dist_init = 0
    min_dist_total = 0
    min_dist_initT = 0
    min_dist_totalT = 0
    to_be_replaced=0

    # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
    
    # Calcul de la matrice des distances
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]
    
    # Tirage aléatoire de P centres différents via une heuristique random
    centres, S = P_median_greedy_random(filename, P)

    #Calcul du cout de la solution actuelle
    for i in 1:n
        min_dist_init = minimum(d[i, k] for k in centres)
        min_dist_total += min_dist_init
    end
        while (nbIter <= maxIter)
            #Si la solution s'est stabilisée alors on tente de trouver une meilleure solution
            if(nb_stable>(1/10)*maxIter)
                centresT, STemp = P_median_greedy_random(filename, P)
                #Calcul du cout de la nouvelle solution
                for i in 1:n
                    min_dist_initT = minimum(d[i, k] for k in centresT)
                    min_dist_totalT += min_dist_initT
                end
                #Meilleure solution trouvée
                if(min_dist_total>min_dist_totalT)
                    centres=centresT
                    S=STemp
                    min_dist_total=min_dist_totalT
                end
            end
            #Amélioration de la solution actuelle: descente stochastique
            if !isempty(setdiff(1:n, centres))
                centres_copy = copy(centres)
                #Selectionner aléatoirement une ville candidate
                nouvelle_ville = rand(setdiff(1:n, centres_copy))
                #Selectionner aléatoirement un centre 
                to_be_replaced = rand(1:length(centres))
                #Supression du centre et ajout de la ville à la solution
                centres_copy[to_be_replaced] = nouvelle_ville
                S_copy = zeros(Int, n)

                for j in centres_copy
                    S_copy[j] = 1
                end
                #Calcul du cout de la nouvelle solution
                dist = 0
                total = 0
                for i in 1:n
                    dist = minimum(d[i, k] for k in centres_copy)
                    total += dist
                end
                #Meilleure solution trouvée
                if (total < min_dist_total)
                    centres=centres_copy
                    S=S_copy
                    min_dist_total = total
                    nb_stable=0
                else
                    nb_stable+=1
                end

            end
            nbIter += 1
            
        end
    return centres,S
    
end

function P_centre_greedy_meta(filename, P,maxIter)
    n = 0
    tabX = Float64[]
    tabY = Float64[]
    f = Float64[]
    nbIter = 0
    nb_stable=0
    to_be_replaced=0
    max_min_dist=0

    # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
    
    # Calcul de la matrice des distances
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]
    
    # Tirage aléatoire de P centres différents via une heuristique random
    centres, S = P_centre_greedy_random(filename, P)

    max_min_dist = maximum([minimum(d[i, j] for j in centres) for i in 1:n])
        while (nbIter <= maxIter)
            #Si la solution s'est stabilisée alors on tente de trouver une meilleure solution
            if(nb_stable>(1/10)*maxIter)
                centresT, STemp = P_centre_greedy_random(filename, P)
                #Calcul du cout de la nouvelle solution
                max_min_distTemp = maximum([minimum(d[i, j] for j in centresT) for i in 1:n])
                #Meilleure solution trouvée
                if(max_min_dist>max_min_distTemp)
                    centres=centresT
                    S=STemp
                    max_min_dist=max_min_distTemp
                end
            end
            #Amélioration de la solution actuelle: descente stochastique
            if !isempty(setdiff(1:n, centres))
                centres_copy = copy(centres)
                #Selectionner aléatoirement une ville candidate
                nouvelle_ville = rand(setdiff(1:n, centres_copy))
                #Selectionner aléatoirement un centre 
                to_be_replaced = rand(1:length(centres))
                #Supression du centre et ajout de la ville à la solution
                centres_copy[to_be_replaced] = nouvelle_ville
                S_copy = zeros(Int, n)

                for j in centres_copy
                    S_copy[j] = 1
                end
                #Calcul du cout de la nouvelle solution
                max_min_distTemp = maximum([minimum(d[i, j] for j in centres_copy) for i in 1:n])
                #Meilleure solution trouvée
                if (max_min_distTemp < max_min_dist)
                    centres=centres_copy
                    S=S_copy
                    max_min_dist = max_min_distTemp
                    nb_stable=0
                else
                    nb_stable+=1
                end

            end
            nbIter += 1
            
        end
    return centres,S
    
end
############################################******METAHEURISTIQUE-Descente-stochastiqique******#####################################################

function P_median_descente_stoch_iteree(filename, P, maxIter, nbLancements)
    n = 0
    tabX = Float64[]
    tabY = Float64[]
    f = Float64[]
     # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
    
    dist_tot = 10e10
    best_centres = []
    best_S = []

    for essai in 1:nbLancements
        centres, S = P_median_greedy_meta(filename, P, maxIter)
        dist_loc = calcul_distance_totale(centres, tabX, tabY, n)
        #Meilleure solution trouvée
        if dist_loc < dist_tot
            dist_tot = dist_loc
            best_centres = centres
            best_S = S
        end
    end
    return best_centres, best_S
end

function P_centre_descente_stoch_iteree(filename, P, maxIter, nbLancements)
    n = 0
    tabX = Float64[]
    tabY = Float64[]
    f = Float64[]
   # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)

    min_dist = 10e10
    best_centres = []
    best_S = []

    for essai in 1:nbLancements
        centres, S= P_centre_greedy_meta(filename, P, maxIter)
        min_dist_loc=calcul_max_distance(centres,tabX,tabY,n)
        #Meilleure solution trouvée
        if min_dist_loc < min_dist
            min_dist = min_dist_loc
            best_centres = centres
            best_S = S
        end
    end

    return best_centres, best_S
end

############################################******Relaxation-lineaire******#####################################################

function PLNE_P_median_relaxation(filename,P)
    n=0
    tabX=Float64[]
    tabY=Float64[]
    f= Float64[]
    # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
    
    # Calcul de la matrice des distances
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]
    
    m = Model(CPLEX.Optimizer)
    @variable(m, 0 <= x[1:n, 1:n] <= 1)  # relaxation linéaire ici

    @objective(m, Min, sum(d[i,j]*x[i,j] for i in 1:n, j in 1:n))
    @constraint(m,sum(x[j,j] for j in 1:n)<=P)
    @constraint(m, [i=1:n], sum(x[i,j] for j in 1:n) == 1)
    @constraint(m, [i=1:n, j=1:n], x[i,j] <= x[j,j])

    optimize!(m)
    status = termination_status(m)

    if status == JuMP.MathOptInterface.OPTIMAL
        println("Borne inférieure z_b (relaxation) = ", objective_value(m))
    else
        println("La relaxation n’a pas trouvé de solution optimale.")
        return
    end
    return objective_value(m)

end

function PLNE_P_centre_relaxation(filename,P)
    n=0
    tabX=Float64[]
    tabY=Float64[]
    f= Float64[]

    # Lecture des données depuis le fichier
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
    
    # Calcul de la matrice des distances
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]
    
    m = Model(CPLEX.Optimizer)
    @variable(m, 0 <= x[1:n, 1:n] <= 1)  # relaxation linéaire ici
    @variable(m, z>=0)

    @objective(m, Min, z)
    @constraint(m, sum(x[j,j] for j in 1:n) <= P)
    @constraint(m, [i=1:n], sum(x[i,j] for j in 1:n) == 1)
    @constraint(m, [i=1:n], sum(d[i,j] * x[i,j] for j in 1:n) <= z)
    @constraint(m, [i=1:n, j=1:n], x[i,j] <= x[j,j])

    optimize!(m)
    status = termination_status(m)

    if status == JuMP.MathOptInterface.OPTIMAL
         println("Borne inférieure z_b (relaxation) = ", objective_value(m))
    else
        println("La relaxation n’a pas trouvé de solution optimale.")
        return
    end

   return objective_value(m)
end

############################################******Fonctions auxiliaires******#####################################################
function calcul_distance_totale(centres,tabX,tabY,n)
    min_dist_total=0
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]
    for i in 1:n
        min_dist_init = minimum(d[i, k] for k in centres)
        min_dist_total += min_dist_init
    end
    return min_dist_total
end

function calcul_max_distance(centres,tabX,tabY,n)
    max_min_dist=0
    if isempty(centres)
        return 0
    end
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]
    max_min_dist = maximum([minimum(d[i, j] for j in centres) for i in 1:n])
    return max_min_dist
end

function calcul_gap(z_heuristique, z_b)
    if z_heuristique == 0
        println("Attention : le coût heuristique est nul, division impossible.")
        return 
    end
    gap = (z_heuristique - z_b) / z_heuristique
    return gap
end

function cout_economique(centres, tabX, tabY, f, n)
    # Calcul du coût total de transport (somme des distances au centre le plus proche)
    d = [dist(tabX[i], tabY[i], tabX[j], tabY[j]) for i in 1:n, j in 1:n]
    total_distance = 0.0
    for i in 1:n
        min_dist = minimum(d[i, k] for k in centres)
        total_distance += min_dist*10
    end
    cout_transport = total_distance * 0.1078  # coût essence/km

    # Calcul du coût total de construction
    cout_construction = sum(f[j] * 10000 for j in centres)

    # Coût global
    cout_total = cout_transport + cout_construction

    println("Coût total transport      : ", round(cout_transport*2, digits=2), " €")
    println("Coût total installation  : ", round(cout_construction, digits=2), " €")
    println("Coût global              : ", round(cout_total, digits=2), " €")
    return cout_total
end  








            
    




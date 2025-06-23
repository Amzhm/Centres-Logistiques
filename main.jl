include("mini.jl")

function experimentation(filename, centres, tabX, tabY, n, P,type)
    if isempty(centres)
        return
    end

    if type == "P1"
        z_b = PLNE_P_median_relaxation(filename, P)
        dist_tot = calcul_distance_totale(centres, tabX, tabY, n)
        println("Borne relaxation : ", z_b)
        println("Gap : ", calcul_gap(dist_tot, z_b))
    elseif type == "P2"
        z_b = PLNE_P_centre_relaxation(filename, P)
        dist_max = calcul_max_distance(centres, tabX, tabY, n)
        println("Borne relaxation : ", z_b)
        println("Gap : ", calcul_gap(dist_max, z_b))
    end
end


function main()
    maxIter=100
    n = 0
    tabX = Float64[]
    tabY = Float64[]
    f = Float64[]
    dist_tot=0
    min_max_dist=0
    prob="P1"
    
    println("=== Problème de localisation de centres ===")
    println("Entrez le nom du fichier (ex: inst_X.flp) :")
    user_input = readline()
    filename = "Instances_Placement/" * strip(user_input)
    println("Lecture du fichier: ", filename)
    n = Lit_fichier_Placement(filename, tabX, tabY, f)
    println("Le fichier contient ", n, " villes")

    println("Combien de centres (P) voulez-vous ouvrir ?")
    P = parse(Int, readline())

    println("\nQuelle méthode souhaitez-vous utiliser ?")
    println("==== Méthodes P-MEDIAN ====")
    println("1. PLNE P-median")
    println("2. PLNE P-median avec coût f")
    println("3. P-median glouton aléatoire")
    println("4. P-median glouton insertion")
    println("5. P-median métaheuristique")
    println("6. P-median descente stochastique itérée")
    println("==== Méthodes P-CENTRE ====")
    println("7. PLNE P-centre")
    println("8. P-centre glouton (villes les plus éloignées)")
    println("9. P-centre glouton randomisé")
    println("10. P-centre métaheuristique")
    println("11. P-centre Descente stochastique itérée")
    print("Votre choix (1-11) : ")
    choix = parse(Int, readline())

    centres = Int[]
    S = zeros(Int, 0)

    if choix == 1
    println("  Problème : P-médian")
    println("Objectif : ouvrir P centres pour minimiser la somme des distances entre chaque ville et son centre.")
    println("Méthode : résolution exacte avec un PLNE via CPLEX.")
    t1 = time()
    centres, S=PLNE_P_median(filename, P)
    t2 = time()
    println("Temps d'exécution du PLNE: ", round(t2-t1, digits=3), " s")
    return

    elseif choix == 2
        println("  Problème : P-médian avec coût d'installation")
        println("Objectif : minimiser la somme des distances + les coûts d'ouverture f[j].")
        println("Méthode : PLNE avec objectif combiné, résolu par CPLEX.")
        t1 = time()
        centres, S=PLNE_P_median_with_f(filename, P)
        t2 = time()
        println("Temps d'exécution du PLNE: ", round(t2-t1, digits=3), " s")
        return

    elseif choix == 3
        prob="P1"
        println("  Problème : P-médian")
        println("Objectif : minimiser la somme des distances.")
        println("Méthode : gloutonne aléatoire. On tire P centres au hasard et on affecte les villes.")
        t1 = time()
        centres, S = P_median_greedy_random(filename, P)
        t2 = time()
        dist_tot=calcul_distance_totale(centres,tabX,tabY,n)
        min_max_dist = calcul_max_distance(centres,tabX,tabY,n)
        println("Temps d'exécution de l'algorithme glouton: ", round(t2-t1, digits=3), " s")
        println("La distance total: ",dist_tot)
        println("La maximum entre une ville est un centre: ",min_max_dist)
    elseif choix == 4
        prob="P1"
        println("  Problème : P-médian")
        println("Objectif : minimiser la somme des distances.")
        println("Méthode : gloutonne par insertion. On ajoute les centres un par un, en minimisant le coût à chaque étape.")
        t1 = time()
        centres, S = P_median_greedy(filename, P)
        t2 = time()
        dist_tot=calcul_distance_totale(centres,tabX,tabY,n)
        min_max_dist = calcul_max_distance(centres,tabX,tabY,n)
        println("Temps d'exécution de l'algorithme glouton: ", round(t2-t1, digits=3), " s")
        println("La distance total: ",dist_tot)
        println("La maximum entre une ville est un centre: ",min_max_dist)
    elseif choix == 5
        prob="P1"
        println("  Problème : P-médian")
        println("Objectif : minimiser la somme des distances.")
        println("Méthode : métaheuristique locale. On modifie les centres aléatoirement et on accepte si c’est mieux.")
        t1 = time()
        centres, S = P_median_greedy_meta(filename, P, maxIter)
        t2 = time()
        dist_tot=calcul_distance_totale(centres,tabX,tabY,n)
        min_max_dist = calcul_max_distance(centres,tabX,tabY,n)
        println("Temps d'exécution de l'algorithme glouton: ", round(t2-t1, digits=3), " s")
        println("La distance total: ",dist_tot)
        println("La maximum entre une ville est un centre: ",min_max_dist)
    elseif choix == 6
        prob="P1"
        println("  Problème : P-médian")
        println("Objectif : minimiser la somme des distances.")
        println("Méthode : Descente stochastique itérée pour P-médian")
        t1 = time()
        centres, S = P_median_descente_stoch_iteree(filename, P, maxIter, 10)
        t2 = time()
        dist_tot=calcul_distance_totale(centres,tabX,tabY,n)
        min_max_dist = calcul_max_distance(centres,tabX,tabY,n)
        println("Temps d'exécution de l'algorithme glouton: ", round(t2-t1, digits=3), " s")
        println("La distance total: ",dist_tot)
        println("La maximum entre une ville est un centre: ",min_max_dist)
    elseif choix == 7
        println("  Problème : P-centre")
        println("Objectif : minimiser la plus grande distance entre une ville et son centre.")
        println("Méthode : PLNE exact avec variable z représentant la distance maximale.")
        t1 = time()
        centres, S=PLNE_P_centre(filename, P)
        t2 = time()
        println("Temps d'exécution du PLNE: ", round(t2-t1, digits=3), " s")
        return

    elseif choix == 8
        prob="P2"
        println("  Problème : P-centre")
        println("Objectif : minimiser la distance maximale.")
        println("Méthode : gloutonne. On sélectionne les villes les plus éloignées des centres déjà choisis.")
        t1 = time()
        centres, S = P_centre_choix_villes_eloignées(filename, P)
        t2 = time()
        dist_tot=calcul_distance_totale(centres,tabX,tabY,n)
        min_max_dist = calcul_max_distance(centres,tabX,tabY,n)
        println("Temps d'exécution de l'algorithme glouton: ", round(t2-t1, digits=3), " s")
        println("La distance total: ",dist_tot)
        println("La maximum entre une ville est un centre: ",min_max_dist)
    elseif choix == 9
        prob="P2"
        println("  Problème : P-centre")
        println("Objectif : minimiser la distance maximale.")
        println("Méthode : gloutonne randomisée. À chaque étape, on ajoute une des villes les plus mal desservies, au hasard.")
        t1 = time()
        centres, S = P_centre_greedy_random(filename, P)
        t2 = time()
        dist_tot=calcul_distance_totale(centres,tabX,tabY,n)
        min_max_dist = calcul_max_distance(centres,tabX,tabY,n)
        println("Temps d'exécution de l'algorithme glouton: ", round(t2-t1, digits=3), " s")
        println("La distance total: ",dist_tot)
        println("La maximum entre une ville est un centre: ",min_max_dist)
    elseif choix == 10
        prob="P2"
        println("  Problème : P-centre")
        println("Objectif : minimiser la distance maximale.")
        println("Méthode : métaheuristique locale. On remplace un centre par un autre si cela améliore la couverture.")
        t1 = time()
        centres, S = P_centre_greedy_meta(filename, P, maxIter)
        t2 = time()
        dist_tot=calcul_distance_totale(centres,tabX,tabY,n)
        min_max_dist = calcul_max_distance(centres,tabX,tabY,n)
        println("Temps d'exécution de l'algorithme glouton: ", round(t2-t1, digits=3), " s")
        println("La distance total: ",dist_tot)
        println("La maximum entre une ville est un centre: ",min_max_dist)
    elseif choix == 11
        prob="P2"
        println("  Problème : P-centre")
        println("Objectif : minimiser la distance maximale.")
        println("Méthode : Descente stochastique itérée pour P-centre")
        t1 = time()
        centres, S = P_centre_descente_stoch_iteree(filename, P, maxIter, 10)
        t2 = time()
        dist_tot=calcul_distance_totale(centres,tabX,tabY,n)
        min_max_dist = calcul_max_distance(centres,tabX,tabY,n)
        println("Temps d'exécution de l'algorithme glouton: ", round(t2-t1, digits=3), " s")
        println("La distance total: ",dist_tot)
        println("La maximum entre une ville est un centre: ",min_max_dist)
    else
        println("Choix invalide.")
        return
    end
    if(choix!= 1 && choix != 2 && choix != 7)
        println("Voulez-vous lancer l'expérimentation ? (o/n)")
        reponse = readline()
        if lowercase(reponse) == "o"
            experimentation(filename, centres, tabX, tabY, n, P, prob)
        end
    end
    # Affichage du dessin de la solution
    Dessine_Solution_Placement(filename, n, tabX, tabY, S)
end

# Lancement automatique
main()

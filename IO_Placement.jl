# Le package Plots va nous permettre de dessiner
# mais il faut auparavant avoir taper: import Pkg; Pkg.add("Plots")

ENV["GKSwstype"] = "nul"  # redirige les sorties de Plots vers les fichiers et non à l'écran
using Plots


# Fonction qui lit un fichier de données et renvoie le nombre de données et la liste des données
# Retourne: n: nombre de villes
# Remplit les tableaux
#    tabX, tabY: coordonnées des villes
#    f: poids des villes



function Lit_fichier_Placement(nom_fichier, tabX, tabY, f)

  n=0
  
  open(nom_fichier) do fic

        for (i,line) in enumerate(eachline(fic)) 
              lg = split(line," ")      # découpe la ligne suivant les espaces et crée un tableau 
              if (i<=1) 
                  n= parse(Int,lg[1])
              end
              if (i>=2) && (i<=n+1) 
                 push!(tabX,parse(Float64,lg[3]))
                 push!(tabY,parse(Float64,lg[4]))	
                 push!(f,parse(Float64,lg[5]))               
              end              
        end    
  end
  return n
end

function Lit_coordonnees_ville(nom_fichier, coordVillesX, coordVillesY)
  
  open(nom_fichier) do fic

        for (i,line) in enumerate(eachline(fic)) 
              lg = split(line," ")      # découpe la ligne suivant les espaces et crée un tableau 
              push!(coordVillesX,parse(Float16,lg[5]))
              push!(coordVillesY,parse(Float16,lg[6]))	
        end    
  end
end


# Fonction Dessine instance
function Dessine_Instance_Placement(nom_fichier)

    n=0
    tabX=Float64[]
    tabY=Float64[]
    f= Float64[]
    
    println("Lecture du fichier: ", nom_fichier)

    n= Lit_fichier_Placement(nom_fichier, tabX, tabY, f)

    println("Le fichier contient ",n, " villes")
  
    nom_fichier_en_deux_morceaux = split(nom_fichier,".") # découpe le nom du fichie d'entrée sans l'extension
    nom_fichier_avec_pdf_a_la_fin= nom_fichier_en_deux_morceaux[1]*".pdf"

    println("Création du fichier pdf de l'instance: ", nom_fichier_avec_pdf_a_la_fin)

    Plots.plot(tabX,tabY,seriestype = :scatter)
    Plots.savefig(nom_fichier_avec_pdf_a_la_fin)  # Ecrit la courbe créée à la ligne précédente dans un fichier .pdf
    
end


# Fonction Dessine solution
function Dessine_Solution_Placement(nom_fichier, n, tabX, tabY, S; couleur_centre=:green, couleur_segment=:red)
    X = Float64[]
    Y = Float64[]	

    nom_fichier_en_deux_morceaux = split(nom_fichier,".") # découpe le nom du fichie d'entrée sans l'extension
    nom_fichier_avec_pdf_a_la_fin = nom_fichier_en_deux_morceaux[1] * "_sol.pdf"

    println("Création du fichier pdf de la solution: ", nom_fichier_avec_pdf_a_la_fin)

    # Affichage des centres et des autres villes avec des couleurs différentes
    Plots.scatter(tabX[S .== 1], tabY[S .== 1], color=couleur_centre, label="Centres")
    Plots.scatter!(tabX[S .!= 1], tabY[S .!= 1], color=:blue, label="Villes")

    for i = 1:n
        min = 10e10
        minj = 0
        # Trouve la ville la plus proche de i dans S
        for j = 1:n
            if ( (S[j] == 1) && (min > dist(tabX[i], tabY[i], tabX[j], tabY[j])) )
                min = dist(tabX[i], tabY[i], tabX[j], tabY[j])                          
                minj = j             
            end
        end
        if (i != minj)
            empty!(X)
            empty!(Y)			
            push!(X, tabX[i])
            push!(X, tabX[minj])	
            push!(Y, tabY[i])
            push!(Y, tabY[minj])	          
            # dessine le segment reliant i à son centre
            Plots.plot!(X, Y, color=couleur_segment, legend=false)
        end
    end

    Plots.savefig(nom_fichier_avec_pdf_a_la_fin)  # Écrit la courbe dans un fichier .pdf
end

#Fonction qui calcule la distance euclidienne entre deux points
function dist(x1,y1,x2,y2)
  return sqrt((x1-x2)^2+(y1-y2)^2)
end









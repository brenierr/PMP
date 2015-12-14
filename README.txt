Ensimag 3A - IF

    BOUANANI Aida 
    BRENIER Roxane 
    MESSAOUDI Samia 
    TONNOIR Alexis 

Projet de Modélisation et Programmation

Compiler le projet : 

	1) créer un répertoire build à la racine du projet : mkdir build 
	2) déplacez-vous dans build : cd build
	3) lancer la commande : cmake -DPNL_ROOT=/matieres/5MMPMP5/pnl-1.7.2/build ..
	4) lancer la compile toujours dans build : make

Dossiers : 
	-> src/ contient toutes les sources .cpp et .hpp

	-> data/ contient des fichiers de tests au format .dat
	
	-> tests/ 
		- pricer : main à utiliser en spécifiant le fichier en entrée, option -c pour profit and loss
		- global_test_prices : 
			teste le prix et l'intervalle de confiance pour les fichiers de data/ et un call
			avec des tolérances passées en paramètres 
		- global_test_delta : 
			teste les deltas en 0 et en t!= 0 pour les fichiers data/basket_1.dat, data/basket_2.dat et un call
			avec des tolérances passées en paramètres 
		- global_test_couverture : 
			teste la couverture des options contenues dans les fichiers de data/ et un call
			avec un H un h et une tolérance passés en paramètres
		- tests unitaires : 
			test_asian_parser
			mc-pricer
			test_delta
			test_shift
			test_memory => (à utiliser avec valgrind) teste la libération de la mémoire sur chaque fonction individuellement
			test_couverture => teste sur un call, un asian et un basket, un nombre N de simu de P&L et permet de récupérer les valeurs dans un fichiers pour tracer un histogramme (très long)
		- run_tests.sh :
		  	Lance les global_tests 

	Exemples d'utilisation : ./global_test_prices 0.1 0.1 
				 avec :
				   1) la tolerance pour le prix
				   2) la tolerance pour l'intervalle de confiance 


				./global_test_delta 0.01 0.01 0.01 0.1 
				avec : 
				   1) la tolerance pour le prix
				   2) la tolerance pour l'intervalle de confiance 
				   3) la tolerance pour les deltas
				   4) fdStep = h 

				./global_test_couverture 0.3 0.001 
				avec : 
				   1) la tolerance relative pour profit and loss 
				   2) fdStep = h 

	-> doc/ 
		- compiler la documentation : lancer la commande "make" dans le dossier doc/ puis ouvrir index.html avec votre navigateur. 
		- gantt_pmp.pdf : Diagramme de GANTT du projet.


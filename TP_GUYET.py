#!/usr/bin/python2.7
# -*-coding:Utf-8 -*
#GUYET Ulysse
#lancement du script:
#TP_GUYET.py <nom_fichier.fasta> <seuil>

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SearchIO
import matplotlib.pyplot as plt
from pandas import *
import networkx as nx
import sys, os

if __name__=="__main__":
	#récupération des variables entrées en paramètre (fichier fasta et seuil)
	input_fasta_file = sys.argv[1];
	seuil = sys.argv[2];

	#vérification de la présence de l'existance du fichier entré en paramètre
	if (not os.path.isfile(input_fasta_file)):
		print(input_fasta_file +": no such file!");

	#lancement de la commande BLAST (en query et subject: le fichier fasta entré en paramètre, en sortie: un fichier XML, une e-value fixé à 0.001)
	blastn_cline = NcbiblastnCommandline(query=input_fasta_file, subject=input_fasta_file, evalue=0.001, outfmt=5, out="results.xml")
	#récupération du résultat du blast en sortie standard
	stdout, stderr = blastn_cline()

	#initialisation de la matrice de scores
	matr = {}
	#initialisation du dictionnaire de labels (pour le graphe)
	dict_labels = {}
	i = 0
	#initialisation des colonnes pour la matrice de scores et ajout de données dans le dictionnaire de labels
	for k in SeqIO.parse(input_fasta_file, 'fasta'):
		matr["%s"%k.id] = {}
		dict_labels.update({i:"%s"%k.id})
		i += 1

	#parser le résultat Blast contenu dans le fichier XML
	blast_records = SearchIO.parse('results.xml', 'blast-xml')
	#Pour chaque alignement on récupère le score, et le nom des 2 séquences alignées (subject et query)
	for record in blast_records:
		infos = record[0]
		score = infos.seq_len / infos[0].bitscore
		query = infos.query_id
		subject = infos.description
		subject_modifie = subject.split(' ', 1)[0]
		#on vérifie que l'alignement ne se fait pas avec la même séquence ou que le score ne soit pas en dessous du seuil entré par l'utilisateur
		# si ces conditions sont validées, on enregistre le score dans la matrice sinon on met 0 en valeur de score dans la matrice
		if ((query != subject_modifie) and (score > float(seuil))):
			matr["%s"%query]["%s"%subject_modifie] = score
		else:
			matr["%s"%query]["%s"%subject_modifie] = 0


	#conversion de la matrice de scores en matrice Pandas
	df = DataFrame(matr).T.fillna(0)

	#afficher la matrice
	print df

	#création du graphe
	G=nx.Graph()
	G=nx.from_numpy_matrix(df.values)
	pos=nx.spring_layout(G)
	colors=range(len(G.edges()))
	nx.draw(G,pos,node_color='#FF00FF',edge_color=colors,width=4,edge_cmap=plt.cm.PuRd, with_labels=True, labels=dict_labels, font_size=16)
	#enregistrement en PNG
	plt.savefig("alignement.png")
	#affichage à l'écran
	plt.show()

#sources:
#http://biopython.org/
#http://stackoverflow.com/
#https://networkx.github.io/
#https://www.biostars.org/
#http://matplotlib.org/
#http://pandas.pydata.org/
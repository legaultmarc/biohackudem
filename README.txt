Le fichier train.tsv contient 2000 exemples (lignes) pour entrainer votre classifieur. 
Il est constitué de 5 colonnes séparées par des tabulations:

Col 1 : le nom de l'échantillon
Col 2 : le tissue d'origine de la tumeur (sa classe)
Col 3 : une liste des gènes contenant des mutations somatiques (missense, nonsense ou frameshift), séparés par des "|".
Col 4 : une liste des régions chromosomiques dupliquées (#copies > 2N, autosomes seulement), séparées par des "|".
Col 5 : une liste des régions chromosomiques délétées (#copies < 2N, autosomes seulement), séparées par des "|".

Le fichier validation_masked.tsv est celui que vous utiliserez pour faire vos prédictions. 
Il est organisé de la même façon, sauf que la classe a été masquée. 
Vous devez donc remplacer les "?" de la colonne 2 par vos prédictions (sans altérer le reste du fichier).
À la fin de l'activité, vous m'enverrez ce fichier, renommé en ajoutant le nom de votre équipe, à l'addresse suivante : 
mathieu.lajoie@gmail.com 

Amusez-vous bien!

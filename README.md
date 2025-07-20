# Concordancier multilingue


Cet outil produit des concordances multilingues à partir de corpus structurés en XML-TEI, alignés au segment et annotés lexicalement. 

Des tables de concordance au format HTML sont produites pour chaque requête à partir d'un texte cible ou d'un texte source.

## Structure des fichiers d'entrée

L'outil fonctionne actuellement sur des corpus binaire source/cible. Il requiert deux fichiers TEI préalablement alignés.
L'unité d'alignement, c'est-à-dire l'élément qui contient le segment aligné, doit être indiqué. Par défault, il s'agit de
l'élément `cl` (clause). L'alignement est réalisé par un jeu d'attributs `@xml:id` et `@corresp`. La gestion des alignements 
de `1 > n`, `n > 1` et `n > n` segments est possible.

## Fonctionnement

Peuvent être interrogés les formes ou les lemmes: 

`python3 get_translations.py -s test_data/Val_S.xml -t test_data/Rome_W.xml -q "politia" -w 1 -o test_results -me cl`

Cette commande produit une table de concordance en prenant le texte `Val_S.xml` comme source, avec un contexte de 1 segment à gauche et à droite pour la source et la cible, sur la requête `politia`


## Sortie

Sont produites des tables au format HTML, CSV et LaTeX: ![alt text](img/exemple.png)
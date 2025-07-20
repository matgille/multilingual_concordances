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


---

# Multilingual concordancer



This tool produces multilingual concordances from corpora structured in XML-TEI, segment-aligned and lexically annotated. 


HTML-formatted concordance tables are produced for each query from target or source text.


## Input file structure


The tool currently works on binary source/target corpora. It requires two pre-aligned TEI files.

The alignment unit, i.e. the element containing the aligned segment, must be specified. The `cl` (clause) is the default element. Alignment must be indicated by a set of `@xml:id` and `@corresp` attributes. 
The tool can manage alignments of `1 > n`, `n > 1` and `n > n` segments.


## How it works


Forms and lemmas can be queried: 


`python3 get_translations.py -s test_data/Val_S.xml -t test_data/Rome_W.xml -q "politia" -w 1 -o test_results -me cl`


This command produces a concordance table taking the text `Val_S.xml` as source, with a context of 1 segment left and right for source and target, on the lemma `politia`.



## Output


HTML, CSV and LaTeX formatted tables are produced: ![alt text](img/exemple.png)

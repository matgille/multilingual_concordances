# Concordancier multilingue


Cet outil produit des concordances multilingues à partir de corpus structurés en XML-TEI, alignés au segment et annotés lexicalement. 

Des tables de concordance au format HTML sont produites pour chaque requête à partir d'un texte cible ou d'un texte source. Peuvent être interrogés les formes ou les lemmes: 

`python3 get_translations.py -s data/Val_S.aligned.xml -t data/Rome_W.aligned.xml -q "politia" -w 1 -o test_result`

Cette commande produit une table de concordance en prenant le texte Val_S.aligned.xml comme source, avec un contexte de 1 segment à gauche et à droite pour la source et la cible, sur la requête `politia`

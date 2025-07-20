import json

import lxml.etree as ET
import subprocess
import re
import os
import collections
import torch
import argparse
import pandas as pd

tei_namespace = 'http://www.tei-c.org/ns/1.0'
ns_decl = {'tei': tei_namespace}


def tree_to_dict_of_nodes(tree: ET._ElementTree, search_lemmas, window, minimal_element) -> dict:
    dictionnary = {}
    all_clauses = tree.xpath(f"descendant::tei:{minimal_element}", namespaces=ns_decl)
    for idx, clause in enumerate(all_clauses):
        localisation: str = clause.xpath("ancestor::tei:div/@n", namespaces=ns_decl)[-1]
        ident: str = clause.xpath("@xml:id")[0]
        corresp: list = clause.xpath("@corresp")[0].replace("#", "").split()
        current_tokens = " ".join(
            clause.xpath("descendant::node()[self::tei:pc or self::tei:w]/descendant::text()", namespaces=ns_decl))
        if idx == 0:
            tokens = current_tokens
        else:
            previous_tokens = " | ".join([" ".join(all_clauses[index].xpath("descendant::node()[self::tei:pc or self::tei:w]/descendant::text()",
                                             namespaces=ns_decl)) for index in range(idx - window, idx)])
            try:
                next_tokens = " | ".join([" ".join(all_clauses[index].xpath("descendant::node()[self::tei:pc or self::tei:w]/descendant::text()",
                                             namespaces=ns_decl)) for index in range(idx + 1, idx + window + 1)])
            except IndexError:
                next_tokens = ""
            tokens = previous_tokens + " " + current_tokens + " " + next_tokens
        if search_lemmas:
            searchable_tokens = " ".join(clause.xpath("descendant::node()[self::tei:pc or self::tei:w]/@lemma",
                                             namespaces=ns_decl))
        else:
            searchable_tokens = current_tokens
        dictionnary[ident] = {"corresp": corresp, "tokens": tokens, "searchable_tokens": searchable_tokens,
                              "localisation": localisation}
    with open("/home/mgl/Documents/test.json", "w") as output_json:
        json.dump(dictionnary, output_json)
    return dictionnary


def match_all_nodes(source_nodes, target_nodes, query, out_dir="test_results"):
    print(f"Searching for {query}")
    corresp_sents = {}
    result = []
    for idx, (ident, node) in enumerate(source_nodes.items()):
        if query in node['searchable_tokens']:
            correponding_sents = [target_nodes[corr] for corr in node['corresp']]
            corresp_sents[idx] = (node['corresp'], node['tokens'], correponding_sents)
            print("---")
            print(node['localisation'])
            print(f"Segment: {ident}\n{node['tokens']}")
            print([sent['tokens'] for sent in correponding_sents])
            print(f"Corresp segments: {', '.join(node['corresp'])}")
            result.append([node['localisation'],
                           ident,
                           node['tokens'],
                           " | ".join([sent['tokens'] for sent in correponding_sents]),
                           ", ".join(node['corresp'])])
    print(f"{len(result)} results found for query '{query}'.")
    with open(f"{out_dir}/{query}.tsv", "w") as result_file:
        result_file.write("Localisation\tId. source\tSegment source\tSegment cible\tId. cible\n")
        for item in result:
            result_file.write("\t".join(item) + "\n")

    with open(f"{out_dir}/{query}.tex", "w") as result_file_tex:
        table_string = pd.DataFrame([[ligne[0]] + ligne[2:4] for ligne in result]).to_latex(index=False, header=False)
        table_string = table_string.replace("\\\\", "\\\\\hline")
        table_string = table_string.replace("\\toprule", "")
        table_string = table_string.replace("\\bottomrule", "")
        table_string = table_string.replace("lll", "|p{1cm}|p{6.5cm}|p{6.5cm}|")
        table_string = table_string.replace("\midrule", "\hline")
        result_file_tex.write(table_string)
        
    tsv2html(f"{out_dir}/{query}.tsv", f"{out_dir}/{query}.html")


def tsv2html(tsv_file_name, html_file_name):
    df = pd.read_csv(tsv_file_name, sep='\t', header=0)
    old_width = pd.get_option('display.max_colwidth')
    # pd.set_option('display.max_colwidth', -1)
    with open(html_file_name, 'w') as html_file:
        html_file.write(df.to_html(index=False))
    pd.set_option('display.max_colwidth', old_width)


def normalize_spelling(text_file):
    '''
    On normalise avec des expressions régulières avant d'envoyer à la lemmatisation, pour adapter au mieux
    aux normes attendues par Freeling.
    IN: *tokenized.txt
    OUT: *.tokenized.normalized.txt
    '''

    # Attention à l'ordre.
    table_normalisation = collections.OrderedDict({
        re.compile(r'^rr(.*)'): r'r\g<1>',
        re.compile(r'(.*)mm(.*)'): r'\g<1>m\g<2>',
        # Provoque ne erreur dans `sujudgado`
        # re.compile(r'uj([aeiouyáéíóúý])'): r'vi\g<1>',
        re.compile(r'lv([aeiouyáéíóúý])'): r'lu\g<1>',
        re.compile(r'ç([ei])'): r'c\g<1>PROUTPROUT',
        re.compile(r'ç([ao])'): r'z\g<1>',
        re.compile(r'ç'): r'z',
        re.compile(r'([aeiouyáéíóúý])v([aeiouyáéíóúý])'): r'\g<1>u\g<2>',
    })

    normalisation = {'á': 'a',
                     'é': 'e',
                     'í': 'i',
                     'ó': 'o',
                     'ú': 'u',
                     'ý': 'y',
                     "⁊": "e",
                     "ç": "ç",
                     'çe': 'ce',
                     'ça': 'za',
                     'çi': 'ci',
                     'sci': 'ci',
                     'sce': 'ce',
                     'ꝑten': 'perten',
                     'ꝑte': 'parte',
                     'ꝑt': 'part',
                     'ꝑ': 'per',
                     'ço': 'zo',
                     "ç": "z"}

    with open(text_file, "r") as input_text_file:
        list_of_words = {index: line.replace('\n', '') for index, line in
                         enumerate(input_text_file.readlines())}

    text = "\n".join([form for index, form in list_of_words.items()])
    for orig, reg in normalisation.items():
        text = text.replace(orig, reg)

    with open(text_file.replace('.txt', '.normalized.txt'), "w") as output_text_file:
        output_text_file.write(text)
        output_text_file.write("\n")
    print("Done")


def lemmatisation(fichier, langue):
    """
        Lemmatisation du fichier XML et réinjection dans le document xml originel.
        :param temoin: le temoin à lemmatiser
        :param division: la division à traiter
        """
    moteur_transformation = "saxon9he.jar"
    fichier_sans_extension = fichier.split("/")[-1].replace(".xml", "")
    print(fichier_sans_extension)
    fichier_xsl = "transformation_pre_lemmatisation.xsl"
    chemin_vers_fichier = fichier
    fichier_entree_txt = f'data/{fichier_sans_extension}.tokenized.txt'
    fichier_normalized = f'data/{fichier_sans_extension}.tokenized.normalized.txt'
    param_sortie = f"sortie={fichier_entree_txt}"
    subprocess.run(["java", "-jar",
                    moteur_transformation,
                    chemin_vers_fichier,
                    fichier_xsl,
                    param_sortie])

    parser = ET.XMLParser(load_dtd=True,
                          resolve_entities=True)
    f = ET.parse(fichier, parser=parser)
    root = f.getroot()
    groupe_words = f"//node()[self::tei:w|self::tei:pc]"
    tokens = root.xpath(groupe_words, namespaces=ns_decl)
    for token in tokens:
        # Résolution d'un pb de tokens vides ou ne contenant que des espaces
        if token.text == " ":
            print("removing")
            token.getparent().remove(token)
    print(fichier_entree_txt)
    tokens = root.xpath(groupe_words, namespaces=ns_decl)
    with open(fichier_entree_txt, "w") as output_file:
        output_file.write("\n".join([token.text.replace("\n", "") for token in tokens]))

    if langue == "es":
        normalize_spelling(fichier_entree_txt)
        fichier_lemmatise = f'data/{fichier_sans_extension}.lemmatized.txt'
        texte_lemmatise = txt_to_liste(fichier_lemmatise)
        fichier_lemmatise = fichier
        n = 1
        for index, mot in enumerate(tokens):
            print("Replacing")
            # Ça marche bien si la lemmatisation se fait
            # sans retokenisation. Pour l'instant, ça bloque avec les chiffre (ochenta mill est fusionné). Voir
            # avec les devs de Freeling.
            try:
                _, lemme_position, pos_position, *autres_analyses = texte_lemmatise[index]
            except Exception as ecxp:
                print(f"Error in file {fichier}: \n {ecxp}. \n Last token of the file must be a "
                      f"punctuation mark. Otherwise, you should search for nested tei:w. Be careful of not adding "
                      f"any tei:w in the header !")
                print("You should check for empty tei:w in tokenized files.")
                print(mot.xpath("@xml:id"))
                print(f"{mot}, {[previous_token.text for previous_token in tokens[index - 10: index]]}")
            mot.set("lemma", lemme_position)
            mot.set("pos", pos_position)

    elif langue == "la":
        modele_latin = "model.tar"
        device = torch.device("cuda") if torch.cuda.is_available() else "cpu"
        cmd = f"pie tag --device {device} {fichier_entree_txt} " \
              f"<{modele_latin},lemma,pos,Person,Numb,Tense,Case,Mood>  --batch_size 2048"
        print(cmd)
        # subprocess.run(cmd.split())
        fichier_seul = os.path.splitext(fichier_entree_txt)[0]
        fichier_lemmatise = str(fichier_seul) + "-pie.txt"
        maliste = txt_to_liste(fichier_lemmatise)
        # Nettoyage de la liste
        maliste.pop(0)  # on supprime les titres de colonne
        test_list = list(zip([token.text for token in tokens], maliste))
        for idx, (token, analysis) in enumerate(test_list):
            if token.lower().strip() != analysis[0]:
                print(idx)
                print(token, analysis)
                exit(0)
        for index, mot in enumerate(tokens):
            liste_correcte = maliste[index]
            _, cas, mode, nombre, personne, temps, lemme, pos, *autres_arguments = liste_correcte

            # on nettoie la morphologie pour supprimer les entrées vides
            morph = f"CAS={cas}|MODE={mode}|NOMB.={nombre}|PERS.={personne}|TEMPS={temps}"
            morph = re.sub("((?!\|).)*?_(?=\|)", "", morph)  # on supprime les pipes non renseignés du milieu
            morph = re.sub("^\|*", "", morph)  # on supprime les pipes qui commencent la valeur
            morph = re.sub("(\|)+", "|", morph)  # on supprime les pipes suivis
            morph = re.sub("\|((?!\|).)*?_$", "", morph)  # on supprime les pipes non renseignés de fin
            morph = re.sub("(?!\|).*_(?!\|)", "", morph)  # on supprime les pipes non renseignés uniques
            #
            mot.set("lemma", lemme)
            mot.set("pos", pos)
            if morph:
                mot.set("morph", morph)

    with open(fichier, "w+") as sortie_xml:
        a_ecrire = ET.tostring(root, pretty_print=True, encoding='utf-8', xml_declaration=True).decode(
            'utf8')
        sortie_xml.write(str(a_ecrire))
    print("%s ✓" % fichier)


def txt_to_liste(filename):
    """
    Transforme le fichier txt produit par Freeling ou pie en liste de listes pour processage ultérieur.
    :param filename: le nom du fichier txt à transformer
    :return: une liste de listes: pour chaque forme, les différentes analyses
    """
    output_list = []
    fichier = open(filename, 'r')
    for ligne in fichier.readlines():
        if not re.match(r'^\s*$',
                        ligne):  # https://stackoverflow.com/a/3711884 élimination des lignes vides (séparateur de phrase)
            ligne_splittee = re.split(r'\s+', ligne)
            output_list.append(ligne_splittee)
    return output_list


def main(source, target, query, lang, lemmatize, search_forms, window, out_dir, minimal_element):
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass
    if lemmatize:
        lemmatisation(source, lang)
    source_as_tree = ET.parse(source)
    target_as_tree = ET.parse(target)
    search_lemmas = not search_forms
    nodes_source = tree_to_dict_of_nodes(source_as_tree, search_lemmas=search_lemmas, window=window, minimal_element=minimal_element)
    nodes_target = tree_to_dict_of_nodes(target_as_tree, search_lemmas=False, window=window, minimal_element=minimal_element)
    match_all_nodes(nodes_source, nodes_target, query, out_dir)
    # print(nodes_source)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--target", default="",
                        help="Document cible (pour lequel on recherche les traductions)")
    parser.add_argument("-s", "--source", default="",
                        help="Document source (duquel on cherche les traductions).")
    parser.add_argument("-q", "--query", default="",
                        help="Requête (lemme ou partie de lemme).")
    parser.add_argument("-w", "--window", default=1,
                        help="Fenêtre de contexte à droite et à gauche.")
    parser.add_argument("-l", "--lang", default="",
                        help="Fenêtre de contexte à droite et à gauche.")
    parser.add_argument("-o", "--out_dir", default="",
                        help="Répertoire de sortie des résultats.")
    parser.add_argument("-me", "--minimal_element", default="cl",
                        help="Élément servant à l'alignement du corpus.")
    parser.add_argument('--lemmatize', action=argparse.BooleanOptionalAction)
    parser.add_argument('--search_forms', action=argparse.BooleanOptionalAction)
    args = parser.parse_args()
    source = args.source
    out_dir = args.out_dir
    target = args.target
    query = args.query
    minimal_element = args.minimal_element
    window = int(args.window)
    lang = args.lang
    lemmatize = args.lemmatize
    search_forms = args.search_forms
    main(source=source,
         target=target,
         query=query,
         lang=lang,
         lemmatize=lemmatize,
         search_forms=search_forms,
         window=window,
         out_dir=out_dir,
         minimal_element=minimal_element)

import json
import random
import lxml.etree as ET
import subprocess
import re
import os
import collections
import torch
import argparse
import pandas as pd
import copy
import corpus_query_language as CQL


tei_namespace = 'http://www.tei-c.org/ns/1.0'
ns_decl = {'tei': tei_namespace}


def tree_to_dict_of_nodes(tree: ET._ElementTree, window, minimal_element, save_as) -> dict:
    dictionnary = {'ids': [], 'idents': {}}
    all_clauses = tree.xpath(f"descendant::tei:{minimal_element}", namespaces=ns_decl)
    if window != 0:
        delimiter = "|"
    else:
        delimiter = ""
    for idx, clause in enumerate(all_clauses):
        localisation: str = clause.xpath("ancestor::tei:div/@n", namespaces=ns_decl)[-1]
        ident: str = clause.xpath("@xml:id")[0]
        dictionnary["ids"].append(ident)
        corresp: list = clause.xpath("@corresp")[0].replace("#", "").split()
        all_tokens = clause.xpath("descendant::node()[self::tei:pc or self::tei:w]", namespaces=ns_decl)
        lemmas = [token.xpath("@lemma")[0] for token in all_tokens]
        pos = [token.xpath("@pos")[0] for token in all_tokens]
        morph = [token.xpath("@morph")[0] if len(token.xpath("@morph")) != 0 else '' for token in all_tokens]
        text = [token.text for token in all_tokens]
        annotations = []
        for word, lemma, pos, morph in list(zip(text, lemmas, pos, morph)):
            annotations.append({"word": word, "lemma": lemma, "pos": pos, "morph": morph})
        current_tokens = " ".join(
            clause.xpath("descendant::node()[self::tei:pc or self::tei:w]/descendant::text()", namespaces=ns_decl))
        dictionnary["idents"][ident] = {"corresp": corresp, "sent": current_tokens, "annotations": annotations,
                              "localisation": localisation}

    with open(save_as, "w") as output_json:
        json.dump(dictionnary, output_json)
    return dictionnary




def match_all_nodes(source_nodes, target_nodes, query, window, out_dir="test_results", filter="", negative_filter="", reduce=1, debug=False):
    results_name = f"{query}"
    if filter != "":
        results_name += f"-{filter}"
    if negative_filter != "":
        results_name += f"-not({negative_filter})"
    corresp_sents = {}
    result = []
    CQLParser = CQL.core.CQLEngine()
    for idx, (ident, node) in enumerate(source_nodes["idents"].items()):
        corresp_segments = ', '.join(node['corresp'])
        if CQLParser.match(node['annotations'], query, verbose=False, debug=False):
            corresponding_sents = [target_nodes["idents"][corr] for corr in node['corresp']]
            # Adapter la fonction de filtre
            all_corresp_annotations = [annotation for corr in node['corresp'] for annotation in target_nodes["idents"][corr]['annotations']]
            if filter != "" and not CQLParser.match(all_corresp_annotations, filter, verbose=False, debug=False):
                continue
            if negative_filter != "" and CQLParser.match(all_corresp_annotations, negative_filter, verbose=False, debug=False):
                continue

            first_source_node_id = ident
            index = next(idx for idx, id in enumerate(source_nodes["ids"]) if id == first_source_node_id)
            previous_source_nodes = " ".join([source_nodes["idents"][source_nodes["ids"][rolling_idx]]["sent"] for rolling_idx in range(index - window, index)])
            try:
                next_source_nodes = " ".join([source_nodes["idents"][source_nodes["ids"][rolling_idx]]["sent"] for rolling_idx in range(index + 1, index + 1 + window)])
            except IndexError:
                next_source_nodes = ""

            source_nodes_with_context = previous_source_nodes + " <b> " + node['sent'] + " </b> " + next_source_nodes

            # On doit gérer différemment la cible, car il y a fusion possible (ce qui n'est pas le cas dans la source
            # (l'outil ne gère pas le 2 > 1 segments pour l'instant)
            try:
                first_target_node_id = node['corresp'][0]
            except IndexError:
                result.append([node['localisation'],
                               ident,
                               source_nodes_with_context,
                               "ø",
                               "ø"])
                continue
            last_target_node_id = node['corresp'][-1]
            first_corresp_index = next(idx for idx, id in enumerate(target_nodes["ids"]) if id == first_target_node_id)
            last_corresp_index = next(idx for idx, id in enumerate(target_nodes["ids"]) if id == last_target_node_id)
            previous_target_nodes = " ".join([target_nodes["idents"][target_nodes["ids"][rolling_idx]]["sent"] for rolling_idx in range(first_corresp_index - window, first_corresp_index)])
            try:
                next_target_nodes = " ".join([target_nodes["idents"][target_nodes["ids"][rolling_idx]]["sent"] for rolling_idx in range(last_corresp_index + 1, last_corresp_index + 1 + window)])
            except IndexError:
                next_target_nodes = ""
            corresponding_sents_as_txt = " | ".join([sent["sent"] for sent in corresponding_sents])
            corresponding_sents_with_context = previous_target_nodes + " <b> " + corresponding_sents_as_txt + " </b> " + next_target_nodes



            corresp_sents[idx] = (node['corresp'], node["sent"], corresponding_sents)
            if debug:
                print(node['localisation'])
                print(f"Segment: {ident}\n{source_nodes_with_context}")
                print(corresponding_sents_with_context)
                print(f"Corresp segments: {corresp_segments}")
                print("---")
            result.append([node['localisation'],
                           ident,
                           source_nodes_with_context,
                           corresponding_sents_with_context,
                           ", ".join(node['corresp'])])
    print(f"{len(result)} results found for query '{query}'.")
    steps = int(1 / reduce)
    result = result[::steps]
    print(len(result))

    with open(f"{out_dir}/{results_name}.tsv", "w") as result_file:
        result_file.write("Localisation\tId. source\tSegment source\tSegment cible\tId. cible\n")
        for item in result:
            result_file.write("\t".join(item) + "\n")

    with open(f"{out_dir}/{results_name}.tex", "w") as result_file_tex:
        table_string = pd.DataFrame([[ligne[0]] + ligne[2:4] for ligne in result]).to_latex(index=False, header=False)
        table_string = table_string.replace("\\\\", "\\\\\hline")
        table_string = table_string.replace("\\toprule", "")
        table_string = table_string.replace("\\bottomrule", "")
        table_string = table_string.replace("lll", "|p{1cm}|p{6.5cm}|p{6.5cm}|")
        table_string = table_string.replace("\midrule", "\hline")
        table_string = table_string.replace("<b>", "\\textbf{")
        table_string = table_string.replace("</b>", "}")
        result_file_tex.write(table_string)
        
    tsv2html(f"{out_dir}/{results_name}.tsv", f"{out_dir}/{results_name}.html", nombre_resultat=len(result))


def tsv2html(tsv_file_name, html_file_name, nombre_resultat):
    df = pd.read_csv(tsv_file_name, sep='\t', header=0)
    old_width = pd.get_option('display.max_colwidth')
    # pd.set_option('display.max_colwidth', -1)
    with open(html_file_name, 'w') as html_file:
        written = df.to_html(index=False)
        html_file.write(written.replace("&lt;", "<").replace("&gt;", ">"))
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

    text = "\n".join([word for index, word in list_of_words.items()])
    for orig, reg in normalisation.items():
        text = text.replace(orig, reg)

    with open(text_file.replace('.txt', '.normalized.txt'), "w") as output_text_file:
        output_text_file.write(text.replace(".", ".\n"))
        output_text_file.write("\n")
    print("Done")


def lemmatisation(fichier, langue, out_dir):
    """
        Lemmatisation du fichier XML et réinjection dans le document xml originel.
        :param temoin: le temoin à lemmatiser
        :param division: la division à traiter
        """
    moteur_transformation = "scripts/saxon9he.jar"
    fichier_sans_extension = fichier.split("/")[-1].replace(".xml", "")
    print(fichier_sans_extension)
    fichier_xsl = "scripts/transformation_pre_lemmatisation.xsl"
    chemin_vers_fichier = fichier
    fichier_sortie_txt = f'{out_dir}/{fichier_sans_extension}.tokenized.txt'
    param_sortie = f"sortie={fichier_sortie_txt}"
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
    print(fichier_sortie_txt)
    tokens = root.xpath(groupe_words, namespaces=ns_decl)
    with open(fichier_sortie_txt, "w") as output_file:
        output_file.write("\n".join([token.text.replace("\n", "") for token in tokens]))
        
    if langue == "es":
        normalize_spelling(fichier_sortie_txt)
        fichier_normalise = f'{out_dir}/{fichier_sans_extension}.tokenized.normalized.txt'
        fichier_lemmatise = f'{out_dir}/{fichier_sans_extension}.lemmatized.txt'
        cmd_sh = ["sh",
                  "scripts/analyze.sh",
                  fichier_normalise,
                  fichier_lemmatise]  # je dois passer par un script externe car un subprocess tourne dans le vide,
        # pas trouvé pourquoi
        subprocess.run(cmd_sh)
        texte_lemmatise = txt_to_liste(fichier_lemmatise)
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
        modele_latin = "scripts/model.tar"
        modele_voice = "scripts/model_voice.tar"
        print(fichier_sortie_txt)
        device = torch.device("cuda") if torch.cuda.is_available() else "cpu"
        cmd = f"pie tag --device {device} {fichier_sortie_txt} " \
              f"<{modele_latin},lemma,pos,Person,Numb,Case><{modele_voice},Mood_Tense_Voice>  --batch_size 2048"
        # subprocess.run(cmd.split())
        fichier_seul = os.path.splitext(fichier_sortie_txt)[0]
        fichier_lemmatise = str(fichier_seul) + "-pie.txt"
        maliste = txt_to_liste(fichier_lemmatise)
        # Nettoyage de la liste
        maliste.pop(0)  # on supprime les titres de colonne
        for index, mot in enumerate(tokens):
            liste_correcte = maliste[index]
            _, cas, mode, temps, voix, nombre, personne, lemme, pos, *autres_arguments = liste_correcte
            assert len(liste_correcte) == 9, "Length error"
            # on nettoie la morphologie pour supprimer les entrées vides
            morph = f"CAS={cas}|NOMB.={nombre}|PERS.={personne}|TEMPS={temps}||MODE={mode}|VOIX={voix}"
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
    Transworde le fichier txt produit par Freeling ou pie en liste de listes pour processage ultérieur.
    :param filename: le nom du fichier txt à transworder
    :return: une liste de listes: pour chaque worde, les différentes analyses
    """
    output_list = []
    fichier = open(filename, 'r')
    for ligne in fichier.readlines():
        if not re.match(r'^\s*$',
                        ligne):  # https://stackoverflow.com/a/3711884 élimination des lignes vides (séparateur de phrase)
            ligne_splittee = re.split(r'[\s+\|]', ligne)
            if ligne_splittee[2] == '_':
                ligne_splittee[2:3] = ["_", "_", "_"]
            output_list.append(ligne_splittee[:-1])
    return output_list


def main(source,
         target,
         query,
         lang,
         lemmatize,
         window,
         out_dir,
         minimal_element,
         filter,
         negative_filter,
         reduce,
         create_json_dict):

    source_as_json = source.replace(".xml", ".json")
    target_as_json = target.replace(".xml", ".json")
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass
    if lemmatize:
        print("Lemmatizing")
        lemmatisation(source, lang, out_dir)
    source_as_tree = ET.parse(source)
    target_as_tree = ET.parse(target)
    if create_json_dict:
        nodes_source = tree_to_dict_of_nodes(source_as_tree, window=window, minimal_element=minimal_element, save_as=source_as_json)
        nodes_target = tree_to_dict_of_nodes(target_as_tree, window=window, minimal_element=minimal_element, save_as=target_as_json)
    else:
        with open(source_as_json, "r") as input_source:
            nodes_source = json.load(input_source)
        with open(target_as_json, "r") as input_target:
            nodes_target = json.load(input_target)


    match_all_nodes(nodes_source,
                    nodes_target,
                    window=window,
                    query=query,
                    out_dir=out_dir,
                    filter=filter,
                    negative_filter=negative_filter,
                    reduce=reduce)
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
                        help="Langue.")
    parser.add_argument("-o", "--out_dir", default="",
                        help="Répertoire de sortie des résultats.")
    parser.add_argument("-me", "--minimal_element", default="cl",
                        help="Élément servant à l'alignement du corpus.")
    parser.add_argument("-f", "--filter", default="",
                        help="Filtre sur la cible.")
    parser.add_argument("-nf", "--negative-filter", default="",
                        help="Filtre négatif sur la cible.")
    parser.add_argument("-r", "--reduce", default="1",
                        help="Ne retenir que la proportion proposée d'exemples.")
    parser.add_argument('--lemmatize', action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument('--create_json', action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()
    source = args.source
    out_dir = args.out_dir
    target = args.target
    query = args.query
    filter = args.filter
    negative_filter = args.negative_filter
    minimal_element = args.minimal_element
    window = int(args.window)
    lang = args.lang
    reduce = float(args.reduce)
    lemmatize = args.lemmatize
    create_json_dict = args.create_json

    main(source=source,
         target=target,
         query=query,
         lang=lang,
         lemmatize=lemmatize,
         window=window,
         out_dir=out_dir,
         minimal_element=minimal_element,
         filter=filter,
         negative_filter=negative_filter,
         reduce=reduce,
         create_json_dict=create_json_dict)

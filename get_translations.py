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
        for form, lemma, pos, morph in list(zip(text, lemmas, pos, morph)):
            annotations.append({"form": form, "lemma": lemma, "pos": pos, "morph": morph})
        current_tokens = " ".join(
            clause.xpath("descendant::node()[self::tei:pc or self::tei:w]/descendant::text()", namespaces=ns_decl))
        dictionnary["idents"][ident] = {"corresp": corresp, "sent": current_tokens, "annotations": annotations,
                              "localisation": localisation}

    with open(save_as, "w") as output_json:
        json.dump(dictionnary, output_json)
    return dictionnary

def process_query(query="[pos='AQ.*'][pos='NC.*']"):
    """
    Implémentation d'un parseur CQL. Pour l'instant, ne gère que des requêtes simples (pas de répétition, *+{})
    """
    all_tokens = re.compile(r"\[([^\]\[]+)\]")
    tokenized = [item for item in re.split(all_tokens, query) if item != ""]
    queries = []
    regex_query = re.compile(r"=")
    for item in tokenized:
        field, expression = re.split(regex_query, item)
        expression = f"^{expression}$"
        expression = re.compile(expression.replace("'", ""))
        queries.append((field, expression))

    return queries

def check_if_match(annotation, query):
    result = re.match(query, annotation)
    if result:
        return True
    else:
        return False


def check_if_its_a_match(segment, queries):
    # On vérifie d'abord que le premier matche
    print("---")
    print(segment)
    first_field, first_regexp = queries[0]
    first_query_matches = [check_if_match(token[first_field], first_regexp) for token in segment[:-len(queries)]]
    if any(first_query_matches):
        print(first_query_matches.index(True))
        print(first_query_matches)
        exit(0)

def check_if_path_matches(sentence, segment, queries, matches, matches_number, match_first=False, debug=False):
    """
    Cette fonction permet de vérifier si une requête est vérifiée dans un texte. Fonction récursive.
    :param segment: le segment sous la forme d'une liste de dictionnaire avec les analyses
    :param queries: la requête sous la forme d'une liste [('pos': re.compile('AQ.*')), (etc)]
    :param matches: le nombre de résultats positifs
    :param matches_number: le nombre de résultats positifs recherchés
    """
    # print(f"Segment: {segment}")
    # print(f"Query: {queries}")
    # On itère sur les différentes requêtes (chacune des unités de la requête)
    if "allegar los castiellos" in sentence:
        debug = True
    else:
        return
        debug = False
    if match_first:
        updated_segment = segment[0:1]
    else:
        updated_segment = segment

    if match_first:
        updated_queries = queries[0:1]
    else:
        updated_queries = queries
    for idx_queries, query in enumerate(updated_queries):
        # print(f"Query: {query}")
        field, regexp = query
        # On itère sur chacun des tokens du segment
        for idx_token, token in enumerate(updated_segment):
            # Si le token matche avec la requête, on relance la fonction avec la requête suivante et le token suivant
            if check_if_match(token[field], regexp):
                matches += 1
                # On vérifie qu'on n'est pas en fin de segment ou en fin de requête
                if idx_token + 1 != len(updated_segment) and idx_queries + 1 != len(updated_queries):
                    # print(queries[idx_queries + 1])
                    # On relance la fonction avec les listes tronquées à la position actuelle, la quantité de matchs et
                    # la taille de la requête
                    if check_if_path_matches(sentence, segment[idx_token + 1:], queries[idx_queries + 1:], matches, matches_number, match_first=True):
                        matches += 1
                        if matches == matches_number:
                            # print(f"It's a match: {token} against {regexp}")
                            print("True 1")
                            return True
                    else:
                        if debug:
                            print("Exit 0")
                        return False
                else:
                    # Si on est en fin de segment ou de requête et qu'on a matché le nombre de requêtes, c'est un match
                    if matches == matches_number:
                        print("True 2")
                        return True
                    else:
                        if check_if_path_matches(sentence, segment[idx_token + 1:], queries[idx_queries + 1:], matches,
                                                 matches_number, match_first=True):
                            matches += 1
                            if matches == matches_number:
                                print("True 3")
                                return True
        if debug:
            print("Exit 3")
        return False


def iter_over_sentence(sentence, pattern):
    if sentence == []:
        return
    match = False
    index = 0
    pattern_index_match = 0
    while not match:
        if pattern_index_match == len(pattern):
            return True
        elif index + 1 == len(sentence):
            return
        current_state = sentence[index]
        field, query = pattern[pattern_index_match]
        if check_if_match(current_state[field], query):
            pattern_index_match += 1
            index += 1
        else:
            pattern_index_match = 0
            index += 1

def match_all_nodes(source_nodes, target_nodes, query, window, out_dir="test_results", filter_by_lemma=False, filter_by_form=False, negative_filter=False, filter="", reduce=1):
    if filter_by_lemma or filter_by_form:
        if negative_filter:
            print(f"Searching for {query} and not {filter}")
            results_name = f"{query}_!{filter}"
        else:
            print(f"Searching for {query} and {filter}")
            results_name = f"{query}_{filter}"
    else:
        print(f"Searching for {query}")
        results_name = f"{query}"
    corresp_sents = {}
    result = []
    processed_query = process_query(query)
    for idx, (ident, node) in enumerate(source_nodes["idents"].items()):
        # if check_if_its_a_match(node['annotations'], processed_query):
        if iter_over_sentence(node['annotations'], processed_query):
        # if check_if_path_matches(node['sent'], node['annotations'], processed_query, matches=0, matches_number=len(processed_query), match_first=False):
            corresponding_sents = [target_nodes["idents"][corr] for corr in node['corresp']]
            # Adapter la fonction de filtre
            first_source_node_id = ident
            index = next(idx for idx, id in enumerate(source_nodes["ids"]) if id == first_source_node_id)
            previous_source_nodes = " ".join([source_nodes["idents"][source_nodes["ids"][rolling_idx]]["sent"] for rolling_idx in range(index - window, index)])
            next_source_nodes = " ".join([source_nodes["idents"][source_nodes["ids"][rolling_idx]]["sent"] for rolling_idx in range(index + 1, index + 1 + window)])

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
            next_target_nodes = " ".join([target_nodes["idents"][target_nodes["ids"][rolling_idx]]["sent"] for rolling_idx in range(last_corresp_index + 1, last_corresp_index + 1 + window)])
            corresponding_sents_as_txt = " | ".join([sent["sent"] for sent in corresponding_sents])
            corresponding_sents_with_context = previous_target_nodes + " <b> " + corresponding_sents_as_txt + " </b> " + next_target_nodes



            corresp_sents[idx] = (node['corresp'], node["sent"], corresponding_sents)
            print("---")
            print(node['localisation'])
            print(f"Segment: {ident}\n{source_nodes_with_context}")
            print(corresponding_sents_with_context)
            print(f"Corresp segments: {', '.join(node['corresp'])}")
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
        
    tsv2html(f"{out_dir}/{results_name}.tsv", f"{out_dir}/{results_name}.html")


def tsv2html(tsv_file_name, html_file_name):
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

    text = "\n".join([form for index, form in list_of_words.items()])
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
        device = torch.device("cuda") if torch.cuda.is_available() else "cpu"
        cmd = f"pie tag --device {device} {fichier_sortie_txt} " \
              f"<{modele_latin},lemma,pos,Person,Numb,Tense,Case,Mood>  --batch_size 2048"
        print(cmd)
        subprocess.run(cmd.split())
        fichier_seul = os.path.splitext(fichier_sortie_txt)[0]
        fichier_lemmatise = str(fichier_seul) + "-pie.txt"
        maliste = txt_to_liste(fichier_lemmatise)
        # Nettoyage de la liste
        maliste.pop(0)  # on supprime les titres de colonne
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


def main(source,
         target,
         query,
         lang,
         lemmatize,
         window,
         out_dir,
         minimal_element,
         filter_by_lemma,
         filter_by_form,
         negative_filter,
         filter,
         reduce):
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass
    if lemmatize:
        lemmatisation(source, lang, out_dir)
    source_as_tree = ET.parse(source)
    target_as_tree = ET.parse(target)
    search_lemma_in_target = filter_by_lemma is True
    # nodes_source = tree_to_dict_of_nodes(source_as_tree, window=window, minimal_element=minimal_element, save_as="/home/mgl/Documents/source.json")
    # nodes_target = tree_to_dict_of_nodes(target_as_tree, window=window, minimal_element=minimal_element, save_as="/home/mgl/Documents/target.json")

    with open("/home/mgl/Documents/source.json", "r") as input_source:
        nodes_source = json.load(input_source)
    with open("/home/mgl/Documents/target.json", "r") as input_target:
        nodes_target = json.load(input_target)


    match_all_nodes(nodes_source,
                    nodes_target,
                    window=window,
                    query=query,
                    out_dir=out_dir,
                    filter_by_lemma=filter_by_lemma,
                    filter_by_form=filter_by_form,
                    negative_filter=negative_filter,
                    filter=filter,
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
    parser.add_argument("-r", "--reduce", default="1",
                        help="Ne retenir que la proportion proposée d'exemples.")
    parser.add_argument('--lemmatize', action=argparse.BooleanOptionalAction)
    parser.add_argument('--filter_by_lemma', action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument('--filter_by_form', action=argparse.BooleanOptionalAction, default=False)
    parser.add_argument('--negative_filter', action=argparse.BooleanOptionalAction, default=False)
    args = parser.parse_args()
    source = args.source
    out_dir = args.out_dir
    target = args.target
    query = args.query
    filter_by_lemma = args.filter_by_lemma
    filter_by_form = args.filter_by_form
    negative_filter = args.negative_filter
    filter = args.filter
    minimal_element = args.minimal_element
    window = int(args.window)
    lang = args.lang
    reduce = float(args.reduce)
    lemmatize = args.lemmatize

    if filter_by_lemma:
        assert filter_by_form != filter_by_lemma, "Options de filtre incompatibles"

    main(source=source,
         target=target,
         query=query,
         lang=lang,
         lemmatize=lemmatize,
         window=window,
         out_dir=out_dir,
         minimal_element=minimal_element,
         filter_by_lemma=filter_by_lemma,
         filter_by_form=filter_by_form,
         negative_filter=negative_filter,
         filter=filter,
         reduce=reduce)

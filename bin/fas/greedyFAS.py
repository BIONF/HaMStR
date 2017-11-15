
from operator import itemgetter
import xml.etree.ElementTree as ET
import logging
import collections
import inspect
import os
import sys
import math
from optparse import OptionParser
from collections import defaultdict
from copy import deepcopy
from multiprocessing import Pool

# version
version = "1.5.1"
tmp = inspect.getfile(inspect.currentframe())
expath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
tmp = ""

parser = OptionParser(description="You are running greedyFAS.py version " + str(version) + ".")
parser.add_option("-q", "--query", dest="query", default=expath+"/annotations_out/query", help="Path to the folder containing the xml-files of the query protein set, default is fas/in/query")
parser.add_option("-s", "--seed", dest="seed", default=expath+"/annotations_out/template", help="Path to the folder containing the xml-files of the template protein set, default is fas/in/template")
parser.add_option("-w", "--weights", dest="weights", default=(0.7, 0.0, 0.3), help="Defines how the three scores MS, CS and PS are weighted (MS, CS, PS), the standart is (0.7, 0.0, 0.3)")
parser.add_option("-j", "--jobname", dest="jobname", default="out", help="Defines the name for the outputfile")
parser.add_option("-r", "--ref_proteome", dest="ref_proteome", default=0, help="Path to a reference proteome which can be used for the weighting of features, by default there is no reference proteome used")
parser.add_option("-a", "--raw_output", dest="raw_output", default=0, help="If set to 1, the FAS score will be printed to STDOUT. If 0, scores will be printed into output file (XML format). If 2, both output variants are conducted.")
parser.add_option("-t", "--priority_threshold", dest="priority_threshold", default=50, help="Change to define the threshold for activating priority mode in the path evaluation.")
parser.add_option("-m", "--max_cardinality", dest="max_cardinality", default=5000, help="Change to define the threshold for the maximal cardinality of feature paths in a graph. If max. cardinality is exceeded the priority mode will be used to for the path evaluation.")
parser.add_option("-l", "--log", dest="loglevel", default="ERROR", help="Change the verbosity of the standard output to the screen. Levels: DEBUG, INFO, WARNING, ERROR, CRITICAL")
parser.add_option("-f", "--efilter", dest="efilter", default="0.001", help="E-value filter for hmm based search methods (feature based/complete sequence).")
parser.add_option("-i", "--inst_efilter", dest="inst_efilter", default="0.01", help="E-value filter for hmm based search methods (instance based/complete sequence).")
parser.add_option("-g", "--weightcorrection", dest="weightcorrection", default="loge", help="Function applied to the frequency of feature types during weighting, options are linear(no function), loge(natural logarithm[Default]), log10(base-10 logarithm), root4(4th root) and root8(8th root).")
parser.add_option("-x", "--weight_constraints", dest="weight_constraints", default=0, help="Apply weight constraints via constraints file, by default there are no constraints.")
parser.add_option("-d", "--featuretypes", dest="featuretypes", default=expath+ "/config/featuretypes", help="inputfile that contains the tools/databases used to predict features")
parser.add_option("-e", "--extendedout", dest="extendedout", default=1, help="0: only the scores and paths will be in the output, 1: a file with the architecure for each protein will be created additionally")
(options, args) = parser.parse_args();

### global vars ###             ### global var looks ###
#naming
template_proteome = {}        #{("protein_id", {("domain_name", [("START", "STOP")])})}
query_proteome = {}           #{("protein_id", {("domain_name", [("START", "STOP")])})}
clan_dict = {}                #{("domain_name", "clan")}
#search_features = {}         #{("F_0", ("domain_name", "POSITION", "Start", "Stop" ))}
a_s_f = {}                    # additional search features 
query_features = {}           #{("F_0", ("domain_name", "POSITION", "Start", "Stop" ))}
a_q_f = {}                    # additional query features
query_protein = {}            #{("domain_name", ["POSITION_1", "POSITION_2"])}
query_clans = {}              #{("clan", "INSTANCES")}
#weights = {}                 #{("domain_name", "weight")}
weights_counter = {}
protein_lengths = {}          #{("protein_id", "LENGTH")}
domain_count = {}             #{("domain", "COUNT")} 
MS_uni = 0
mode = {}                     #{("protein", mode)}
constraints = {}              #{("tool/feature", weight)}
weight_const = 0
weight_correction = 1
input_linearized = []
input_normal = []
e_output = 1
## hidden options
# tab separated table as output file #
taciturn = 1
# two-sided linearization #
entire = 1

### READ OPTIONS ###
score_weights = options.weights
p_path = options.seed
s_path = options.query
ref_proteome = options.ref_proteome
priority_threshold = options.priority_threshold
max_cardinality = options.max_cardinality
loglevel = options.loglevel.upper()
efilter = float(options.efilter)
inst_efilter = float(options.inst_efilter)

if ref_proteome == 0:
    MS_uni = 1
else:
    MS_uni = 0
    if options.weightcorrection == "log10":
        weight_correction = 2
    elif options.weightcorrection == "loge":
        weight_correction = 1
    elif options.weightcorrection == "root4":
        weight_correction = 3
    elif options.weightcorrection == "root8":
        weight_correction = 4
    elif options.weightcorrection == "linear":
        weight_correction = 0
if int(options.raw_output) == 0:
    output = 0
elif int(options.raw_output) == 1:
    output = 1
else:
    output = 2
if int(options.extendedout) == 0:
    e_output = 0
elif int(options.extendedout) == 1:
    e_output = 1
else:
    print(str(options.extendedout) + " is not a valid input for -e [--extendedout]")
    quit() 
if options.ref_proteome != 0 and options.weight_constraints != 0:
    weight_const = 1

### SETUP LOGGING OPTIONS ###

# possible logging configurations
# logging into stdout with time stamp:
#logging.basicConfig(level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')
# logging into file with line number:
logging.basicConfig(filename='testlog.log', filemode='w', level=loglevel, format='%(lineno)s - %(levelname)s - %(message)s')
# logging into stdout with line number:
#logging.basicConfig(level=loglevel, format='%(lineno)s - %(levelname)s - %(message)s')

logging.info('greedyFAS.py started with options: entire='+str(entire)+', priority_threshold='+str(priority_threshold)+', log_level='+str(loglevel))
logging.info('score_weights are set to: '+ str(score_weights[0]) +" "+ str(score_weights[1]) +" "+ str(score_weights[2]))
logging.info('ref_proteome is set to: '+str(ref_proteome))

########## Flow Control ########## <fc>
# flow control functions

# main function
def fc_main():
    global p_path
    global s_path
    global tools
    global ref_proteome
    global template_proteome
    global MS_uni
    global weight_correction
    global input_linearized
    global input_normal

    globalscores = {}

    ## MS_uni set to 0 when no weighting is conducted        
    if MS_uni == 0:
            if os.path.exists(ref_proteome + "/count.ref"):
                ref_reader(ref_proteome + "/count.ref")
            else:
                for ftype in input_linearized:
                    xmlreader(ref_proteome + "/" + ftype + ".xml", 2, ftype,True)
                for ftype in input_normal:
                    xmlreader(ref_proteome + "/" + ftype + ".xml", 2, ftype,True)

                ## keep feature counts for reference proteome (gene set)
                w_count_ref()
                ######################################
                # reusing feature counts for reference protein sets
                #
                # ref_writer writes into directory with reference protein set
                # may cause "permission denied" error
                ######################################
                #ref_writer(ref_proteome + "/count.ref")

                ## clean global protome
                template_proteome = {}
            if weight_correction != 0:
                w_weight_correction(weight_correction)    # use correction function on counts
    for ftype in input_linearized:
        xmlreader(p_path + "/" + ftype + ".xml", 0, ftype,True)
        xmlreader(s_path + "/" + ftype + ".xml", 1, ftype,True)
    for ftype in input_normal:
        xmlreader(p_path + "/" + ftype + ".xml", 0, ftype,False)
        xmlreader(s_path + "/" + ftype + ".xml", 1, ftype,False)
    if MS_uni == 0:
        w_count()
    fc_path()


# finds best paths
def fc_path():
    global template_proteome
    global a_s_f
    global query_features
    global a_q_f
    global query_proteome #named search_proteome before
    global search_graph
    global MS_uni
    global mode
    global priority_threshold
    global max_cardinality
    global weight_const
    global outpath
    global protein_lengths
    global output
    global e_output
    global score_weights
    global weight_const

    logging.info("fc_path")

    #score model M2 - used for iteration/incremental evaluation
    ##############################################################################################################################

    ######## entire mode ############################################
    # both proteins (search and query) will be linearized ##########
    search_features = {}
    weights = {}
    if output == 0 or output == 2:
        if e_output == 1:
            a_out = open(outpath + "_architecure", "w+")
            a_out.write("<?xml version=\"1.0\"?>\n")
            a_out.write("<architectures>\n")
        out = open(outpath, "w+")
        out.write("<?xml version=\"1.0\"?>\n")
        if MS_uni == 0:
            out.write("<out weighting=\"applied\">\n")
        else:
            out.write("<out weighting=\"uniform\">\n")
    for query in query_proteome:
        go_priority = False
        if MS_uni == 0:
            if weight_const == True:
                weights = w_weighting_constraints(query, 0) 
            else:
                weights = w_weighting(query, 0)
        lin_single_set = su_lin_query_protein(query)
        tmp_single_graph = pb_entire_main_graph(lin_single_set)

        #PRIORITY CHECK 1: checking for number of instances - assess complexity of the feature graph
        if int(len(query_features)) > int(priority_threshold):
            go_priority = True
        else:
            #creating all paths
            all_single_paths = pb_graphtraversal(tmp_single_graph, 1)

            #PRIORITY CHECK 2: checking for number of paths in the feature graph
            if (int(len(all_single_paths)) > int(max_cardinality)):
                logging.info("Switched to priority mode due to number of paths.")
                go_priority = True

        if output == 0 or output == 2:
            out.write("\t<query id=\"" + query + "\" length=\"" + str(int(protein_lengths["single_"+query])) + "\">\n")
        elif(taciturn == 1):
            out = open(outpath, "w+")
        su_query_protein(query)
        for protein in template_proteome:

            mode[protein] = 0;
            pathcount = 0
            if MS_uni == 0:
                if weight_const == True:
                    weights = w_weighting_constraints(protein, 1) 
                else:                        
                    weights = w_weighting(protein, 1)
            if go_priority == True:
                all_single_paths = []
                for domain_type in query_proteome[query]:
                    all_single_paths += (pb_entire_graphtraversal_priority(tmp_single_graph, domain_type, protein, 1, search_features, weights))
                    logging.info("domain_type: "+str(domain_type))
                    logging.debug("all_single_paths: "+str(all_single_paths))       
            # max fixture of (search_path, score, query_path)
            max_fixture = ([], (0.0, 0.0, 0.0, 0.0, 0,0, False), [], protein)

            #check for available paths
            #baustelle: calculation
            # SINGLE(query) <--VS-- SET(search)
            #case M2.1: empty(query) 
            if int(len(all_single_paths)) == 0:
                logging.warning("CASE M2.1: No paths (pfam or smart annotations) in query (single).")

                search_protein, search_features, a_s_f = su_search_protein(protein)
                logging.warning("search_features(p): "+str(len(search_features))+" for: "+str(protein))

                # case M2.1.1: empty(query)-empty(search)
                # should be the best fix independent from weight
                if int(len(search_features)) == 0:
                    logging.warning("CASE M2.1.1: empty vs empty.")
                    score_w = sf_entire_calc_score(a_s_f.keys(),a_q_f.keys(), score_weights, weight_const, weights, search_features, a_s_f, query_features, a_q_f, MS_uni)
                    path = a_s_f.keys()
                    query_architecture = a_q_f.keys()
                    mode[protein] = 2
                else:
                    # case M2.1.2: empty(query)-graph(search)
                    logging.warning("CASE M2.1.2: empty vs graph.")
                    tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, [], search_features, weights)
                    path = tmp_path_score[0][0]
                    score_w = tmp_path_score[0][1]
                    query_architecture = tmp_path_score[0][2]
                    mode[protein] = 2

                # set max_fixture according to 
                max_fixture = (path, score_w, query_architecture, protein)

            # handle all paths
            # case M2.2: graph(query)
            for query_path in all_single_paths:
                logging.warning("CASE M2.2: graph.")
                pathcount += 1
                search_protein, search_features, a_s_f = su_search_protein(protein)
                logging.warning("search_features(p): "+str(len(search_features))+" for: "+str(protein))
                logging.debug("query_path No."+str(pathcount)+": "+str(query_path))

                #case M2.2.1: graph(query)-empty(search)
                if int(len(search_features)) == 0:
                    logging.warning("CASE M2.2.1: graph vs empty.")
                    #special case: protein with no pfam or smart domains
                    # get score for a_s_f and query_path directly
                    path = a_s_f.keys()
                    query_path_ad = query_path+a_q_f.keys()
                    score_w = sf_entire_calc_score(path,query_path_ad, score_weights, weight_const, weights, search_features, a_s_f, query_features, a_q_f, MS_uni)

                else:
                    #case M2.2.2 graph(query)-graph(search)
                    # regular traversal of graph based on search_protein
                    logging.warning("CASE M2.2.2: graph vs graph.")
                    tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, query_path, search_features, weights)
                    path = tmp_path_score[0][0]
                    score_w = tmp_path_score[0][1]
                    query_path_ad = tmp_path_score[0][2]
                    mode[protein] = tmp_path_score[1]
                    logging.debug("tmp_path_score "+str(tmp_path_score))#search path, score, mode

                # check for max scoring fixture of path and query_path
                if score_w[5] == True:
                    if score_w[3] >= max_fixture[1][3]:
                        max_fixture = (path, score_w, query_path_ad, protein)
                        logging.warning("max_fixture: "+str(max_fixture))
                else:
                    if (max_fixture[1][5] != True) and (score_w[4] >= max_fixture[1][4]):
                        max_fixture = (path, score_w, query_path_ad, protein)

            logging.info("Found: "+str(pathcount)+" path(s) for query.")
            logging.debug("Path fixture: "+str(max_fixture))
            #quit()

            score = max_fixture[1]
            if output == 0 or output == 2:
                try:
                    if mode[protein] == 2:
                        mode_out = "greedy"
                    elif mode[protein] == 1:
                        mode_out = "priority"
                    else:
                        mode_out = "exhaustive"
                except KeyError: 
                    mode_out = "default"

                out.write("\t\t<template id=\"" + protein + "\" score=\"" + str(score[3]) + "\" MS=\"" + str(score[0]) + "\" PS=\"" + str(score[1]) + "\" CS=\"" + str(score[2]) + "\" LS=\"" + str(score[6]) + "\" length=\"" + str(int(protein_lengths["set_"+str(protein)])) + "\" mode=\"" + mode_out + "\">\n")
                if e_output == 1:
                    a_out.write("\t<template id=\"" + protein + "\" length=\"" + str(int(protein_lengths["set_"+str(protein)])) + "\">\n")
                    a_out.write("\t\t<architecture>\n")

                    for feature in template_proteome[protein]:
                        if MS_uni == 0:
                            a_out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\""+str(template_proteome[protein][feature][1])+"\" weight=\"" + str(weights[feature]) + "\">\n")
                        else:
                            a_out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\""+str(template_proteome[protein][feature][1])+"\" weight=\"" + str(1.0/len(template_proteome[protein])) + "\">\n")
                        for instance in template_proteome[protein][feature][2:]:
                            a_out.write("\t\t\t\t<instance inst_eval=\""+str(instance[0])+"\" start=\"" + str(instance[1]) + "\" end=\"" + str(instance[2]) + "\"/>\n")
                        a_out.write("\t\t\t</feature>\n")
                    a_out.write("\t\t</architecture>\n")
                    a_out.write("\t</template>\n")

            best_template_path = []
            best_query_path = []

            for feature in max_fixture[0]:
                if feature in search_features:
                    logging.debug(str(search_features[feature][0])+" "+str(search_features[feature][1])+" "+str(search_features[feature][2])+" "+str(search_features[feature][3]))
                    best_template_path.append((search_features[feature][0], search_features[feature][1], search_features[feature][2], search_features[feature][3]))
                else:
                    best_template_path.append((a_s_f[feature][0], a_s_f[feature][1], a_s_f[feature][2], a_s_f[feature][3]))

            for feature in max_fixture[2]:
                if feature in query_features:
                    best_query_path.append((query_features[feature][0], query_features[feature][1], query_features[feature][2], query_features[feature][3]))
                else:
                    best_query_path.append((a_q_f[feature][0], a_q_f[feature][1], a_q_f[feature][2], a_q_f[feature][3]))
            if output == 0 or output == 2:
                path_tmp = {}
                path_tmp_query = {}
                scale = 0
                if weight_const == 1:
                    path_tmp2 = []
                    for feature in best_template_path:
                        if feature[0] not in path_tmp2:
                            path_tmp2.append(feature[0])
                    adjusted_weights = w_weight_const_rescale(path_tmp2, weights, search_features, True)
                    weight_tmp = {}
                    for adj_feature in adjusted_weights:
                        weight_tmp[adj_feature] = weights[adj_feature]
                        weights[adj_feature] = adjusted_weights[adj_feature]
                for feature in best_template_path:
                    if feature[0] in path_tmp:
                        path_tmp[feature[0]].append((feature[2],feature[3]))
                    else:
                        path_tmp[feature[0]] = [(feature[2],feature[3])]
                        if MS_uni == 0:
                            scale += weights[feature[0]]
                for feature in best_query_path:
                    if feature[0] in path_tmp_query:
                        path_tmp_query[feature[0]].append((feature[2],feature[3]))
                    else:
                        path_tmp_query[feature[0]] = [(feature[2],feature[3])]
#                            if MS_uni == 0:
#                                scale += weights[feature[0]]
                logging.debug("path_tmp: "+ str(path_tmp))


                ## unweighted case
                if MS_uni == 0:
                    if scale > 0:
                        scale = 1.0/float(scale)
                    else:
                        scale = 1.0
                ## print path
                out.write("\t\t\t<template_path>\n")
                for feature in path_tmp:
                    if MS_uni == 0:
                        out.write("\t\t\t\t<feature type=\"" + feature + "\" corrected_weight=\"" + str(weights[feature]*scale) + "\">\n")
                    else:
                        out.write("\t\t\t\t<feature type=\"" + feature + "\">\n")
                    for tmp_inst in path_tmp[feature]:
                        out.write("\t\t\t\t\t<instance start=\"" + str(tmp_inst[0]) + "\" end=\"" + str(tmp_inst[1]) + "\"/>\n")
                    out.write("\t\t\t\t</feature>\n")
                out.write("\t\t\t</template_path>\n")
                out.write("\t\t\t<query_path>\n")
                for feature in path_tmp_query:
                    out.write("\t\t\t\t<feature type=\"" + feature + "\">\n")
                    for tmp_inst in path_tmp_query[feature]:
                        out.write("\t\t\t\t\t<instance start=\"" + str(tmp_inst[0]) + "\" end=\"" + str(tmp_inst[1]) + "\"/>\n")
                    out.write("\t\t\t\t</feature>\n")
                out.write("\t\t\t</query_path>\n")
                out.write("\t\t</template>\n")
                if weight_const == 1:
                    for adj_feature in adjusted_weights:
                        weights[adj_feature] = weight_tmp[adj_feature]
            if output == 1 or output == 2:
                print (score[3])
                ## hidden ##
                if (taciturn == 1 and output != 2):
                    out.write(protein+"\t"+str(score[3])+"\n")
        if output == 0 or output == 2:
            out.write("\t</query>\n")
    if output == 0 or output == 2:
        out.write("</out>")
        out.close()
        if e_output == 1:
            for query in query_proteome:
                a_out.write("\t<query id=\"" + query + "\" length=\"" + str(int(protein_lengths["single_"+query])) + "\">\n")
                a_out.write("\t\t<architecture>\n")
                for feature in query_proteome[query]:
                    a_out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\""+str(query_proteome[query][feature][1])+ "\">\n")
                    for instance in query_proteome[query][feature][2:]:
                        a_out.write("\t\t\t\t<instance inst_eval=\""+str(instance[0])+"\" start=\"" + str(instance[1]) + "\" end=\"" + str(instance[2]) + "\"/>\n")

                    a_out.write("\t\t\t</feature>\n")
                a_out.write("\t\t</architecture>\n")
                a_out.write("\t</query>\n")
            a_out.write("</architectures>")
            a_out.close()

########## Start Up Functions ########## <su>

# define of 2-dim-dictionary
def tree():
    return collections.defaultdict(tree)

# initializes query protein
def su_query_protein(protein_id):
    global query_proteome
    global query_protein
    global query_clans
    global protein_lengths
    global clan_dict

    query_protein = {}
    query_clans = {}

    for feature in query_proteome[protein_id]:
        query_protein[feature] = []
        clan = clan_dict[feature]
        if clan in query_clans:
            query_clans[clan] += len(query_proteome[protein_id][feature])-2
        else:
            query_clans[clan] = len(query_proteome[protein_id][feature])-2
        for instance in query_proteome[protein_id][feature][2:]:
            position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["single_"+str(protein_id)])
            position = round(position, 8)
            query_protein[feature].append(position)
                

# initializes the search protein
def su_lin_query_protein(protein_id):
    global template_proteome
    global query_proteome
    global query_clans
    global protein_lengths
    global clan_dict
    global query_features
    global a_q_f
    
    logging.info("su_lin_query_protein")
    
    lin_query_protein = []
    query_clans = {}
    query_features = {}
    a_q_f = {}
    tmp = []
    i = 0

    for feature in query_proteome[protein_id]:
        clan = clan_dict[feature]
        if clan in query_clans:
            query_clans[clan] += len(query_proteome[protein_id][feature])-2
        else:
            query_clans[clan] = len(query_proteome[protein_id][feature])-2

        if query_proteome[protein_id][feature][0] == True:
            for instance in query_proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["single_"+str(protein_id)])
                position = round(position, 8)
                key = "F_" + str(i)
                query_features[key] = (feature, position, instance[1], instance[2])
                #{'F_16': ('smart_RRM_1', 0.43473684, 186, 227), 'F_10': ('pfam_RRM_6', 0.39578947, 151, 225), 'F_17': ('smart_RRM_1', 0.62421053, 261, 332), 'F_11': ('pfam_RRM_6', 0.62210526, 261, 330), 'F_14': ('pfam_RRM_5', 0.64210526, 275, 335), 'F_4': ('smart_RRM', 0.90105263, 394, 462), 'F_5': ('pfam_RRM_1', 0.39578947, 151, 225), 'F_6': ('pfam_RRM_1', 0.62315789, 261, 331), 'F_7': ('pfam_RRM_1', 0.90421053, 399, 460), 'F_0': ('pfam_RRM_occluded', 0.40105263, 148, 233), 'F_1': ('pfam_RRM_occluded', 0.62631579, 259, 336), 'F_2': ('smart_RRM', 0.39684211, 150, 227), 'F_3': ('smart_RRM', 0.62421053, 260, 333), 'F_15': ('pfam_RRM_5', 0.91157895, 402, 464), 'F_12': ('pfam_RRM_6', 0.90210526, 398, 459), 'F_19': ('pfam_Transformer', 0.14526316, 6, 132), 'F_8': ('pfam_RRM_7', 0.37894737, 148, 212), 'F_9': ('pfam_RRM_7', 0.60315789, 258, 315), 'F_18': ('smart_RRM_1', 0.88315789, 378, 461), 'F_13': ('pfam_RRM_5', 0.43263158, 183, 228)}
                tmp.append((key, instance[1]))
                i += 1
        else:
            for instance in query_proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["single_"+str(protein_id)])
                position = round(position, 8)
                key = "O_" + str(i)
                a_q_f[key] = (feature, position, instance[1], instance[2])
                i += 1
    
    #sort  instances
    tmp2 = sorted(tmp, key=itemgetter(1))
    for x in tmp2:
        lin_query_protein.append(x[0])

    return lin_query_protein

# initializes the search protein
def su_search_protein(protein_id):
    global template_proteome

    search_features = {}
    search_protein = []
    a_s_f = {}
    tmp = []
    i = 0
    for feature in template_proteome[protein_id]:
        # True if feature is going to be linearized
        if template_proteome[protein_id][feature][0] == True:
            for instance in template_proteome[protein_id][feature][2:]:
                    position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["set_"+str(protein_id)])
                    position = round(position, 8)
                    key = "F_" + str(i)
                    search_features[key] = (feature, position, instance[1], instance[2])
                    tmp.append((key, instance[1]))
                    i += 1
        else:
            for instance in template_proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["set_"+str(protein_id)])
                position = round(position, 8)
                key = "O_" + str(i)
                a_s_f[key] = (feature, position, instance[1], instance[2])
                i += 1

    #sort  instances
    tmp2 = sorted(tmp, key=itemgetter(1))
    for x in tmp2:
        search_protein.append(x[0])
    return search_protein, search_features, a_s_f
                


    
    
# gives path for output
def su_set_path(jobname):
    global expath

    if not os.path.exists(expath + "/out/"):
        os.makedirs(expath + "/out/")
    if os.path.exists(expath + "/out/" + jobname + ".xml"):
        i = 1
        while os.path.exists(expath + "/out/" + jobname + "_" + str(i) + ".xml"):
            i += 1
        jobname = jobname+"_"+str(i)
    return jobname 


########## Pathbuilding Functions ########## <pb>
# Used for the Pfam/Smart domains 

# build graph for a given protein
def pb_entire_main_graph(current_protein):
    global query_features

    region = pb_region_mapper(current_protein, query_features)
    query_graph = pb_region_paths_nongreedy(region)
    return query_graph
    
# creates graph with all paths 
# retrieves best path from graph traversal function
# returns best path, score and mode
# including priority check
def pb_entire_main_nongreedy(search_protein, protein_id, query_path, search_features, weights):
    global max_cardinality

    logging.info("pb_entire_main_nongreedy")
    
    priority_check = False
        
    region = pb_region_mapper(search_protein, search_features)
    search_graph = pb_region_paths_nongreedy(region)
    logging.debug(region)
    logging.debug(search_graph)
    
    if ((int(len(search_features)) >= int(priority_threshold))):
        mode = 1
        priority_check = True
        # checking best path for all domain types
        path_score = pb_entire_priority_mode(protein_id, query_path, search_graph, search_features, weights)
    else:
        #PRIORITY CHECK "2": for every protein in template_proteome
        tmp_path_set_size = len(pb_graphtraversal(search_graph, 1))
        logging.debug("cardinality of tmp graph: "+str(tmp_path_set_size))
        if ((int(tmp_path_set_size) > int(max_cardinality))):
            mode = 1
            priority_check = True
            path_score = pb_entire_priority_mode(protein_id, query_path, search_graph, search_features, weights)
        
    if (priority_check == False):
        mode = 0
        path_score = pb_entire_graphtraversal(search_graph, query_path, mode, search_features, weights)
        #path = path_score[0]
    search_graph == 0
    return (path_score, mode)

        
#priority mode to reduce number of paths in entire mode
def pb_entire_priority_mode(protein, query_path, search_graph, search_features, weights):
    global template_proteome

    logging.info("pb_entire_priority_mode")
    
    # best_path (Path, (MS_score(float), PS_score(float), CS_score(float), final_score(float), path_weight(float), common_feature(bool)), QPath)
    best_path = ("NULL", (0.0, 0.0, 0.0, 0.0, 0.0, False),"NULL")

    for domain_type in template_proteome[protein]:
        if template_proteome[protein][domain_type][0] == True: 
            path_score_w = pb_entire_graphtraversal_priority(search_graph, domain_type, query_path, 0, search_features, weights)
            if path_score_w[1][5] == True:
                if path_score_w[1][3] >= best_path[1][3]:
                    best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
            else:
                #if so far no path with common features found AND the weight is higher
                if ((best_path[1][5] != True) and (path_score_w[1][4] >= best_path[1][4])):
                    best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
                    
    return best_path    

# creates the different junctions of overlap regions for the non-greedy graph
def pb_region_paths_nongreedy(overlap_map):
    graph = {}
    reached_list = {}
    for feature in overlap_map:
        reached = []
        links = []
        for candidate in feature[1]:
            links.append(candidate)
            reached.append(candidate)
            for i in reached_list[candidate]:
                if i in feature[1]:
                    feature[1].remove(i)
                if i not in reached:
                    reached.append(i)
        reached_list[feature[0]] = reached
        graph[feature[0]] = links
    return graph

# creates an overlap map for a given overlap_region for the query_features. That is for each feature a list of features which do not overlap with it.
def pb_region_mapper(overlap_region, features):
    logging.info("pb_region_mapper")
    logging.debug(overlap_region)

    overlap_map = [("START", overlap_region)]
    for i in range(0, len(overlap_region)):
        end = features[overlap_region[i]][3]
        overlap = []
        x = i+1
        while x < len(overlap_region):
            if end <= features[overlap_region[x]][2]:
                overlap.append(overlap_region[x])
            x = x+1
        overlap.append("END")
        overlap_map.append((overlap_region[i], overlap))
    overlap_map.append(("END", []))
    overlap_map.reverse()
    return overlap_map


# return best scoring path 
def pb_entire_graphtraversal(search_graph, query_path, mode, search_features, weights):
        global a_s_f
        global a_q_f
        global score_weights
        global weight_const
        global query_features
        global MS_uni
    
        logging.info("pb_entire_graphtraversal in mode: "+str(mode))
        logging.debug("search graph: "+str(search_graph))
    
        v_stack = ["START"]
        p_stack = []
        paths = []
        if mode == 0:
            best_path = ([], (0.0, 0.0, 0.0, 0.0, 0.0, False, 0.0),[])
        while v_stack:
            vertex = v_stack.pop()
            if len(p_stack) == 0:
                path = []
            else:
                path = p_stack.pop()
            for next in search_graph[vertex]:
                if next == "END":
                    if mode == 1:
                        paths.append(path)
                    else:
                        path_ad = path + a_s_f.keys()
                        query_path_ad = query_path + a_q_f.keys()
                        score_w = sf_entire_calc_score(path_ad, query_path_ad, score_weights, weight_const, weights, search_features, a_s_f, query_features, a_q_f, MS_uni)
                        logging.info("search path "+str(path_ad)+" in pb_entire_graphtraversal")
                        logging.debug("Score info: "+str(score_w)+" for query_path "+str(query_path_ad)+" in pb_entire_graphtraversal")
                        if score_w[5] == True:
                            if score_w[3] >= best_path[1][3]:
                                best_path = (path_ad, score_w, query_path_ad)
                                logging.debug("new best_path by score: "+str(best_path))
                        else:
                            if (best_path[1][5] != True and score_w[4] >= best_path[1][4]):
                                best_path = (path_ad, score_w, query_path_ad)
                                logging.debug("new best_path by weight: "+str(best_path))                                
                else:
                    v_stack.append(next)
                    p_stack.append(path + [next])
        if mode == 1:
            logging.info("return: "+str(paths))
            return paths
        else:
            logging.info("return: "+str(best_path))
            return best_path

# walks a graph to find all paths
# for greedy it returns all possible paths
def pb_graphtraversal(search_graph, mode):

    logging.info("graphtraversal in mode: "+str(mode))
    logging.debug(search_graph)

    v_stack = ["START"]
    p_stack = []
    paths = []
    if mode == 0:
        best_path = ([], (0.0, 0.0, 0.0, 0.0))
    while v_stack:
        vertex = v_stack.pop()
        if len(p_stack) == 0:
            path = []
        else:
            path = p_stack.pop()
        for next in search_graph[vertex]:
            if next == "END":
                if mode == 1:
                    paths.append(path)
                    logging.info("cardinality: "+str(len(paths)))
                    if len(paths) > max_cardinality:
                        return paths
            else:
                v_stack.append(next)
                p_stack.append(path + [next])
    if mode == 1:
        logging.debug("returning paths: "+str(paths))
        logging.info("cardinality: "+str(len(paths)))
        return paths


# walks the graph in priority-mode to reduce number of paths
def pb_entire_graphtraversal_priority(search_graph, priority, query_path, mode, search_features, weights):
    global query_features
    global a_s_f
    global a_q_f
    global score_weights
    global weight_const
    global MS_uni
    global template_proteome


    logging.info("pb_entire_graphtraversal_priority, mode: "+str(mode))
    logging.debug("search_features: "+str(search_features)+"\nsearch_graph: "+str(search_graph) + "\tpriority: "+str(priority) +"\ta_s_f: "+str(a_s_f)+ "\tquery_path: "+str(query_path) + "\ta_q_f: "+ str(a_q_f))
    v_stack = ["START"]
    p_stack = []
    paths = []
    if mode == 0:
        best_path = ([], (0.0, 0.0, 0.0, 0.0, 0.0, False), [])
    else:
        protein = query_path
    while v_stack:
        vertex = v_stack.pop()
        if len(p_stack) == 0:
            path = []
        else:
            path = p_stack.pop()

        p_found = 0
        p_candidates = []
        for next in search_graph[vertex]:
            if next == "END":
                if mode == 1:
                    paths.append(path)
                    p_found = 2
                else:
                    path_ad = path + a_s_f.keys()
                    query_path_ad = query_path + a_q_f.keys()
                    logging.debug("path_ad: "+str(path_ad))
                    logging.debug("query_path_ad: "+str(query_path_ad))

                    score_w = sf_entire_calc_score(path_ad, query_path_ad, score_weights, weight_const, weights, search_features, a_s_f, query_features, a_q_f, MS_uni)
                    if score_w[5] == True:
                        if score_w[3] >= best_path[1][3]:
                            best_path = (path_ad, score_w, query_path_ad)
                    else:
                        if (best_path[1][5] != True) and (score_w[4] >= best_path[1][4]):
                            best_path = (path_ad, score_w, query_path_ad)
                    p_found = 2
            elif mode == 1:
                if query_features[next][0] == priority:
                    p_candidates.append(next)
                    p_found = 1
            elif mode == 0:
                p_candidates.append(next)
                p_found = 1
        if p_found == 1:
            if len(p_candidates) == 1:
                v_stack.append(p_candidates[0])
                p_stack.append(path + p_candidates)
            else:
                if mode == 1:
                    best_priority_bridger = ("NONE",(0.0, 0.0, 0.0, 0.0, 0.0, False))                        
                    for next in p_candidates:
                        ###baustelle: fixed: path elongation without additives
                        score = sf_calc_score(path+[next], mode, protein, score_weights, weight_const, weights, search_features, query_features, query_protein, template_proteome, MS_uni)
                        if score[3] >= best_priority_bridger[1][3]:
                            best_priority_bridger = (next, score)
                    v_stack.append(best_priority_bridger[0])
                    p_stack.append(path + [best_priority_bridger[0]])
                elif mode == 0:
                    #some kind of greedy strategy: if feature type priority appears more than once
                    best_partial_path = ("NONE", (0.0, 0.0, 0.0, 0.0, 0.0, False))            
                    for next in p_candidates:
                        path_ad = path + a_s_f.keys()
                        query_path_ad = query_path + a_q_f.keys()
                        score_w = sf_entire_calc_score(path_ad+[next],query_path_ad, score_weights, weight_const, weights, search_features, a_s_f, query_features, a_q_f, MS_uni)
                        if score_w[5] == True:
                            if score_w[3] >= best_partial_path[1][3]:
                                best_partial_path = (next, score_w)
                        else:
                            if (best_partial_path[1][5] != True) and (score_w[4] >= best_partial_path[1][4]):
                                best_partial_path = (next, score_w)

                    v_stack.append(best_partial_path[0])
                    p_stack.append(path + [best_partial_path[0]])

        elif p_found == 0:
            if mode == 1:
                best_priority_bridger = ("NONE",(0.0, 0.0, 0.0, 0.0, 0.0, False))                        
                for next in search_graph[vertex]:
                    ###baustelle: fixed: path elongation without additives
                    score = sf_calc_score(path+[next], mode, protein, score_weights, weight_const, weights, search_features, query_features, query_protein, template_proteome, MS_uni)
                    if score[3] >= best_priority_bridger[1][3]:
                        best_priority_bridger = (next, score)
                v_stack.append(best_priority_bridger[0])
                p_stack.append(path + [best_priority_bridger[0]])

            elif mode == 0:
                #some kind of greedy strategy: if feature type (p priority) not found
                best_partial_path = ("NONE", (0.0, 0.0, 0.0, 0.0, 0.0, False))            
                for next in search_graph[vertex]:
                    path_ad = path + a_s_f.keys()
                    query_path_ad = query_path + a_q_f.keys()
                    score_w = sf_entire_calc_score(path_ad+[next],query_path_ad, score_weights, weight_const, weights, search_features, a_s_f, query_features, a_q_f, MS_uni)
                    if score_w[5] == True:
                        if score_w[3] >= best_partial_path[1][3]:
                            best_partial_path = (next, score_w)
                    else:
                        if (best_partial_path[1][5] != True) and (score_w[4] >= best_partial_path[1][4]):
                            best_partial_path = (next, score_w)

                v_stack.append(best_partial_path[0])
                p_stack.append(path + [best_partial_path[0]])
################################
    if mode == 1:
        logging.debug("returning: "+str(paths))
        logging.debug("cardinality: "+str(len(paths)))
        return paths
    elif mode == 0:
        logging.debug(best_path)
        return best_path




########## Scoring Functions ########## <sf>
# calculate the fas score sf_calc_score is the main function


# main scoring function
# calculates the final score from MS, CS and PS
def sf_calc_score(path, mode, protein, score_weights, weight_const, weights, search_features, query_features, query_protein, template_proteome, MS_uni):

    if weight_const == 1:
        adjusted_weights = w_weight_const_rescale(path, weights, search_features, False)
        tmp_weight = {}
        for i in adjusted_weights:
            tmp_weight[i] = weights[i]
            weights[i] = adjusted_weights[i]
    if score_weights[1] != 0.0: #baustelle: correct?
        score_CS = round(sf_CS_score(path, mode, search_features), 10)
    else:
        score_CS = 0.0
    tmp = sf_MS_score(path, mode, protein, search_features, weights, query_protein, MS_uni, template_proteome)
    score_MS = round(tmp[0], 10)
    score_PS = round(sf_PS_score(path, tmp[2], mode, protein, search_features, query_features, query_protein, weights, MS_uni, template_proteome), 10)
    final_score = (score_MS * score_weights[0]) + (score_CS * score_weights[1]) + (score_PS * score_weights[2])
    if weight_const == 1:
        for i in adjusted_weights:
            weights[i] = tmp_weight[i]
    return (score_MS, score_PS, score_CS, final_score)

# main scoring function - entire mode
# calculates the final score from MS, CS and PS
def sf_entire_calc_score(path, query_path, score_weights, weight_const, weights, search_features, a_s_f, query_features, a_q_f, MS_uni):

    if weight_const == 1:
        adjusted_weights = w_weight_const_rescale(path, weights, search_features, False)
        tmp_weight = {}
        for i in adjusted_weights:
            tmp_weight[i] = weights[i]
            weights[i] = adjusted_weights[i]
    if score_weights[1] != 0.0: #baustelle
        score_CS = round(sf_entire_cs_score(path, query_path), 10)
    else:
        score_CS = 0.0
    tmp = sf_entire_ms_score(path, query_path, search_features, a_s_f, query_features, a_q_f, weights, MS_uni)
    score_MS = round(tmp[0], 10)
    path_weight = tmp[3]
    common_feature = tmp[4]
    ps_tmp = sf_entire_ps_score(path, tmp[2], query_path, search_features, a_s_f, query_features, a_q_f, weights, MS_uni)
    score_PS = round(ps_tmp[0], 10)
    score_LS = round(ps_tmp[1], 10)
    final_score = (score_MS * score_weights[0]) + (score_CS * score_weights[1]) + (score_PS * score_weights[2])
    if weight_const == 1:
        for i in adjusted_weights:
            weights[i] = tmp_weight[i]
    return (score_MS, score_PS, score_CS, final_score, path_weight, common_feature, score_LS)

# calculates Clan Score 
def sf_CS_score(path, mode, search_features):
    global clan_dict
    global query_clans

    logging.info("sf_CS_score, mode: "+str(mode))
    logging.debug("search_features: "+str(search_features))
    logging.debug("query_features: "+str(query_features))

    counter = 0.0
    path_clans = {}
    score = 0.0

    #counting clans in path
    #path_clans: contains counts for clans in path
    for i in path:
        if mode == 0:
            feature = search_features[i]
        elif mode == 1:
            feature = query_features[i]
        
        if clan_dict[feature[0]] != "NO_CLAN":
            if clan_dict[feature[0]] in path_clans:
                path_clans[clan_dict[feature[0]]] += 1
            else:
                path_clans[clan_dict[feature[0]]] = 1
    for clan in path_clans:
        counter += 1.0
        if clan in query_clans:
            score += float(path_clans[clan] * query_clans[clan]) / float(max(path_clans[clan], query_clans[clan]) * max(path_clans[clan], query_clans[clan]))
            logging.debug("("+str(path_clans[clan])+" * "+str(query_clans[clan])+") / (max("+str(path_clans[clan])+","+str(query_clans[clan])+") * max("+str(path_clans[clan])+","+str(query_clans[clan])+")")
    if counter == 0:
        score = 0
    else:
        score = score / counter
    return score


# calculates Clan Score 
def sf_entire_cs_score(path, query_path):
    global clan_dict

    logging.debug("search_path: "+str(path))
    logging.debug("query_path: "+str(query_path))

    counter = 0.0
    s_clans = {} # clans from path
    q_clans = {} # clans from query_path

    score = 0.0
    for i in path:
        if i in search_features:
            feature = search_features[i]
            if clan_dict[feature[0]] != "NO_CLAN":
                if clan_dict[feature[0]] in s_clans:
                    s_clans[clan_dict[feature[0]]] += 1
                else:
                    s_clans[clan_dict[feature[0]]] = 1
    for i in query_path:
        if i in query_features:
            feature = query_features[i]
            if clan_dict[feature[0]] != "NO_CLAN":
                if clan_dict[feature[0]] in q_clans:
                    q_clans[clan_dict[feature[0]]] += 1
                else:
                    q_clans[clan_dict[feature[0]]] = 1

    for clan in s_clans:
        counter += 1.0
        if clan in q_clans:
            score += float(s_clans[clan] * q_clans[clan]) / float(max(s_clans[clan], q_clans[clan]) * max(s_clans[clan], q_clans[clan]))
    if counter == 0:
        score = 0
    else:
        score = score / counter
    return score


# calculates Multiplicity Score
def sf_MS_score(path, mode, protein, search_features, weights, query_protein, MS_uni, template_proteome):

    logging.info("sf_MS_score, mode: "+str(mode)+", path: "+str(path))
    search_domains = {}
    scale = 0
    scores = []
    final_score = 0
    for i in path:
        if mode == 0:
            feature = search_features[i]
        elif mode == 1:
            feature = query_features[i]
        if feature[0] in search_domains:
            search_domains[feature[0]] += 1
        else:
            search_domains[feature[0]] = 1
            logging.debug("feature in path : "+str(feature[0]))

    logging.debug("weights: "+str(weights))
    for feature in search_domains:
        if MS_uni == 0 and mode == 0:
            scale += weights[feature]
        else:
            scale += 1
        if mode == 0:
            logging.debug("query_protein: "+ str(query_protein))
            if feature in query_protein:
                s_length = len(query_protein[feature])
                p_score =  float(search_domains[feature] * s_length) / float(max(search_domains[feature], s_length) * max(search_domains[feature], s_length))
                scores.append((feature, p_score))
            else:
                scores.append((feature, 0.0))
        elif mode == 1:
            if feature in template_proteome[protein]:
                s_length = len(template_proteome[protein][feature])-2
                p_score =  float(search_domains[feature] * s_length) / float(max(search_domains[feature], s_length) * max(search_domains[feature], s_length))
                scores.append((feature, p_score))
            else:
                scores.append((feature, 0.0))

    if scale > 0:
        scale = 1.0/float(scale)
    for score in scores:
        if MS_uni == 0 and mode == 0:
            final_score += score[1] * scale * weights[score[0]]
        else:
            final_score += score[1] * scale
    return (final_score, search_domains, scale)


# calculates Multiplicity Score - entire mode
def sf_entire_ms_score(path, query_path, search_features, a_s_f, query_features, a_q_f, weights, MS_uni):

    search_domains = {}
    query_domains = {}
    scale = 0
    scores = []
    final_score = 0
    final_weight = 0
    main_features_s = {}
    main_features_q = []
    common_feature = False

    for i in path:
        if i in search_features:
            feature = search_features[i]
            main_features_s[feature[0]] = True
        else:
            feature = a_s_f[i]

        if feature[0] in search_domains:
            search_domains[feature[0]] += 1
        else:
            search_domains[feature[0]] = 1

    for i in query_path:
        if i in query_features:
            feature = query_features[i]
            main_features_q.append(feature[0])
        else:
            feature = a_q_f[i]
        if feature[0] in query_domains:
            query_domains[feature[0]] += 1
        else:
            query_domains[feature[0]] = 1

    #check for common features
    for y in main_features_q:
        try:
            if main_features_s[y] == True:
                common_feature = True
        except KeyError:
            logging.debug("Feature not found: "+str(y))

    for feature in search_domains:
        if MS_uni == 0:
            scale += weights[feature]
        else:
            scale += 1
        if feature in query_domains:
            s_length = query_domains[feature]
            p_score =  float(search_domains[feature] * s_length) / float(max(search_domains[feature], s_length) * max(search_domains[feature], s_length))
            scores.append((feature, p_score))
        else:
            scores.append((feature, 0.0))
    if scale > 0:
        scale = 1.0/float(scale)
    for score in scores:
        if MS_uni == 0:
            final_score += score[1] * scale * weights[score[0]]
            final_weight += weights[score[0]]
        else:
            final_score += score[1] * scale
    logging.debug("Return entire_ms_score: "+ str(final_score)+", "+ str(search_domains)+", "+ str(scale)+", "+ str(final_weight)+", "+ str(common_feature))
    return (final_score, search_domains, scale, final_weight, common_feature)

# calculates Positional Score
def sf_PS_score(path, scale, mode, protein, search_features, query_features, query_protein, weights, MS_uni, template_proteome):

    count = {}
    final_score = 0.0
    scores = {}
    for i in path:
        if mode == 0:
            feature = search_features[i]
            logging.debug("feature: "+str(feature))
            logging.debug("query_protein: "+str(query_protein))
            if feature[0] in query_protein:
                best_match = 0.0
                scores[feature[0]] = 0.0
                for position in query_protein[feature[0]]:
                    match = 1.0 - float(abs(feature[1] - position))
                    if best_match < match:
                        best_match = match
                scores[feature[0]] += best_match
        elif mode == 1:
            feature = query_features[i]
            #position = ((float(instance[0]) + float(instance[1])) / 2.0) / float(protein_lengths[protein_id])
            if feature[0] in template_proteome[protein]:
                best_match = 0.0
                if not feature[0] in scores:
                    scores[feature[0]] = 0.0
                    count[feature[0]] = 0
                for instance in template_proteome[protein][feature[0]][2:]:
                    pos = (float(instance[1]) + float(instance[2])) / 2.0 / float(protein_lengths["set_"+str(protein)])
                    match = 1.0 - float(abs(feature[1]) - pos)
                    logging.debug(str(float(instance[1]))+" + "+str(float(instance[2]))+" / 2.0 / "+str(float(protein_lengths["set_"+str(protein)]))+" = "+str(pos))
                    if best_match < match:
                        best_match = match
                scores[feature[0]] += best_match                    
                count[feature[0]] += 1

    for f_score in scores:
        if MS_uni == 0 and mode == 0:
            final_score += scores[f_score]/count[f_score] * scale * weights[f_score]
        else:
            final_score += scores[f_score]/count[f_score] * scale

    return final_score

# calculates Positional Score - entire mode
def sf_entire_ps_score(path, scale, query_path, search_features, a_s_f, query_features, a_q_f, weights, MS_uni):

    count = {}
    final_score = 0.0
    final_ls_score = 0.0
    scores = {}
    ls_scores = {}
    local_query_protein = {}
    # get current features from query path
    for i in query_path:
        if i in query_features:
            feature = query_features[i]
        else:
            feature = a_q_f[i]
        if feature[0] in local_query_protein:
            length = feature[3]-feature[2]
            local_query_protein[feature[0]].append((feature[1], length))
        else:
            local_query_protein[feature[0]] = [(feature[1], feature[3]-feature[2])]
    logging.debug(str(local_query_protein)+" for query_path "+str(query_path))
            
            
    # compare features in path with features form query path
    for i in path:
        if i in search_features:
            feature = search_features[i]
        else:
            feature = a_s_f[i]
        if feature[0] in local_query_protein:
            best_match = 0.0
            if not feature[0] in scores:
                count[feature[0]] = 0
                scores[feature[0]] = 0.0
                ls_scores[feature[0]] = 0.0
            for instance in local_query_protein[feature[0]]:
                match = 1.0 - float(abs(feature[1] - instance[0]))
                if best_match < match:
                    best_match = match
                    if instance[1] <= (feature[3] - feature[2]):
                        ls = float(instance[1]) / float(feature[3] - feature[2])
                    else:
                        ls = float(feature[3] - feature[2]) / float(instance[1])

            scores[feature[0]] += best_match
            ls_scores[feature[0]] += ls
            count[feature[0]] += 1     
   
    for f_score in scores:
        if MS_uni == 0:
            final_score += scores[f_score]/count[f_score] * scale * weights[f_score]
            final_ls_score += ls_scores[f_score]/count[f_score] * scale * weights[f_score]
        else:
            final_score += scores[f_score]/count[f_score] * scale
            final_ls_score += ls_scores[f_score]/count[f_score] * scale

    return (final_score, final_ls_score)

########## Weighting ########## <w>
# weighting functions

# counts all domains in query_proteome and template_proteome
def w_count():
    global query_proteome
    global domain_count
    global template_proteome

    logging.debug("w_count: counting domains in single protein and protein set.")
    for query in query_proteome:
        for feature in query_proteome[query]:
            if feature not in domain_count:
                domain_count[feature] = 1
    for template in template_proteome:
        protein = template_proteome[template]
        for feature in protein:
            if feature not in domain_count:
                domain_count[feature] = 1

    logging.debug("domain counts: "+str(domain_count))


# counts all domains in a reference pcroteome if one is given
def w_count_ref():
        global domain_count
        global template_proteome
        logging.debug("w_count_ref: counting domains in reference gene set.")

        for i in template_proteome:
            protein = template_proteome[i]
            for feature in protein:
                # counting instances (substract 2 because first and second entry contain assess(bool) and feat_eval(float))
                count = len(protein[feature]) - 2
                if feature in domain_count:
                    domain_count[feature] +=  count
                else:
                    domain_count[feature] = count
        logging.debug("domains counts: "+str(domain_count))



# calculates weights 
def w_weighting(protein, template):
    global domain_count
    global template_proteome
    global query_proteome

    weights = {}
    scaling_factor = 0.0
    sum_of_features = 0.0
    features = []
    if template == 0:
        for feature in query_proteome[protein]:
            features.append(feature)
    else:
        for feature in template_proteome[protein]:
            features.append(feature)
    for feature in features:
        try:
            domain_count[feature]
        except:
            #print "feature not found in ref " + feature
            domain_count[feature] = 1
            sum_of_features += float(domain_count[feature])
        else:
            sum_of_features += float(domain_count[feature])
    for feature in features:
        weights[feature] = round(float(sum_of_features) / float(domain_count[feature]), 8)
        scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 8)
    for feature in features:
        weights[feature] = round(float(weights[feature]) / float(scaling_factor), 8)
    return weights


# calculates weights if constraints are applied (doesn't work with different tools yet)        
def w_weighting_constraints(protein, template):
    global domain_count
    global template_proteome
    global query_proteome
    global constraints

    weights = {}
    tools = {}
    for tool in (input_linearized + input_normal):
        tools[tool] = []
    features = []
    single_constraints = []
    filled = 0.0
    if template == 0:
        for feature in query_proteome[protein]:
            features.append(feature)
    else:
        for feature in template_proteome[protein]:
            features.append(feature)
    for feature in features:
        if feature in constraints:
            filled += constraints[feature]
            weights[feature] = constraints[feature]
            single_constraints.append(feature)
        elif feature.split('_')[0] in constraints:
            tools[feature.split('_')[0]].append(feature)
    for tool in tools:
        if len(tools[tool]) > 0:
            filled += constraints[tool]
            sum_of_features = 0
            scaling_factor = 0.0
            for feature in tools[tool]:
                features.remove(feature)
                try:
                    domain_count[feature]
                except:
                    domain_count[feature] = 1
                    sum_of_features += float(domain_count[feature])
                else:
                    sum_of_features += float(domain_count[feature])
            for feature in tools[tool]:
                weights[feature] = round(float(sum_of_features) / float(domain_count[feature]), 8)
                scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 8)  
            for feature in tools[tool]:
                weights[feature] = round(float(weights[feature]) / float(scaling_factor) * constraints[tool], 8)
    for feature in single_constraints:
        features.remove(feature)
    for feature in features:
        sum_of_features = 0.0
        try:
            domain_count[feature]
        except:
            #print "feature not found in ref " + feature
            domain_count[feature] = 1
            sum_of_features += float(domain_count[feature])
        else:
            sum_of_features += float(domain_count[feature])
    scaling_factor = 0.0
    for feature in features:
        weights[feature] = round(float(sum_of_features) / float(domain_count[feature]), 8)
        scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 8)
    for feature in features:
        weights[feature] = round(float(weights[feature]) / float(scaling_factor) * (1.0-filled), 8)
    return weights

        
              
                        
# weight correction method
def w_weight_correction(method):
    global domain_count

    if method == 1:
        for feature in domain_count:
            domain_count[feature] = int(round(math.log(domain_count[feature]), 0) + 1)
    elif method == 2:
        for feature in domain_count:
            domain_count[feature] = int(round(math.log10(domain_count[feature]), 0) + 1)
    elif method == 3:
        for feature in domain_count:
            domain_count[feature] = int(round(math.pow(domain_count[feature], 0.25), 0))         
    elif method == 4:
        for feature in domain_count:
            domain_count[feature] = int(round(math.pow(domain_count[feature], 0.125), 0))         
            
# rescales weights          
def w_weight_const_rescale(path, weights, search_features, final):
    global constraints
    global input_linearized

    lindict = {}
    for ftype in input_linearized:
        lindict[ftype] = []
    tmp = 0.0
    rescaled_weights = {}
    
    if final:
        for feature in path:
            if feature not in constraints:
                if feature not in lindict[feature.split('_')[0]]:
                    lindict[feature.split('_')[0]].append(feature)
        for ftype in input_linearized:
            if ftype in constraints:
                for feature in lindict[ftype]:
                    tmp += weights[feature]
                if tmp < constraints[ftype] and tmp > 0.0:
                    scale = constraints[ftype] / tmp
                    for feature in lindict[ftype]:
                        rescaled_weights[feature] = weights[feature] * scale 
                tmp = 0.0
    else:
        for feature in path:
            if feature in search_features:
                if search_features[feature][0] not in constraints:
                    if feature not in lindict[search_features[feature][0].split('_')[0]]:
                        lindict[search_features[feature][0].split('_')[0]].append(feature)
        
        for ftype in input_linearized:
            if ftype in constraints:
                for feature in lindict[ftype]:
                    tmp += weights[search_features[feature][0]]
                if tmp < constraints[ftype] and tmp > 0.0:
                    scale = constraints[ftype] / tmp
                    for feature in lindict[ftype]:
                        rescaled_weights[search_features[feature][0]] = weights[search_features[feature][0]] * scale 
                tmp = 0.0
    
    return rescaled_weights


                    
########## Input ########## <OK>

# reads featuretype file
def featuretypes(path):
    global input_linearized
    global input_normal

    ifile = open(path, "r+")
    lines = ifile.readlines()
    mode = "NULL"
    for line in lines:
        tmp = line.rstrip("\n")
        if tmp == "#linearized":
            mode = "lin"
        elif tmp == "#normal":
            mode = "nor"
        elif mode == "NULL":
            print(path + " is not a valid input file")
            quit()
        elif mode == "lin":
            input_linearized.append(tmp)
        elif mode == "nor":
            input_normal.append(tmp)
    ifile.close()

# reads constraints file
## BAUSTELLE: number of annotation tools 
def constraints_in(path):
    global constraints

    cfile = open(path, "r+")
    lines = cfile.readlines()
    i = 1
    if lines [0][0] == "#":
        while lines[i][0] != "#":
            split = (lines[i].rstrip("\n")).split(" ")
            if split[1] != "N":
                constraints[split[0]] = float(split[1])
            i += 1
    else:
        print ("\n" + path + " might be in the wrong format. Please see the sample file in config directory.")
        quit()
    if lines[i][0] == "#":
        i += 1
        while i < len(lines):
            split = (lines[i].rstrip("\n")).split(" ")
            constraints[split[0]] = float(split[1])
            i += 1
    cfile.close()

                                

# reads the xml-files as input
def xmlreader(path, single, tool, assess):
    global query_proteome
    global template_proteome
    global protein_lengths
    global inst_efilter
    global efilter

    if os.path.exists(path):
        xmltree = ET.parse(path) 
        root = xmltree.getroot()
        for protein in root:
            pID = protein.attrib["id"]
            plength = protein.attrib["length"]
            #set up of protein IDs, differentiate between proteins from different files
            if single == 1:
                singleID = "single_"+str(pID)
                protein_lengths[singleID] = float(plength)
            elif single == 0:
                setID = "set_"+str(pID)
                protein_lengths[setID] = float(plength)
            elif single == 2:
                protein_lengths[pID] = float(plength)
            #set up of datastructure to store annotations
            if (single == 0 or single == 2) and not (pID in template_proteome):
                template_proteome[pID]= {}
            elif single == 1 and not (pID in query_proteome):
                query_proteome[pID] = {} 

            for feature in protein:
                if len(feature)>0:
                    ftype = tool + "_" + feature.attrib["type"]
                    feat_eval = 'NULL'

                    # evalue check: family/ftype based
                    if ('evalue' in feature.attrib and float(feature.attrib["evalue"]) > efilter):
                        logging.debug("reject feature type: "+ ftype)
                        #print "reject feature"
                        #skip current ftype and continue with the next one
                        continue
                    else:
                        #keep feature type bases evalue
                        if('evalue' in feature.attrib):
                            feat_eval=float(feature.attrib["evalue"])

                        if ('clan' in feature.attrib):
                            fclan = feature.attrib["clan"]
                        else:
                            fclan = "NO_CLAN"
                        if single == 0 or single == 2:
                            template_proteome[pID][ftype] = []
                            template_proteome[pID][ftype].append(assess)
                            template_proteome[pID][ftype].append(feat_eval)
                        else:
                            query_proteome[pID][ftype] = [] 
                            query_proteome[pID][ftype].append(assess)
                            query_proteome[pID][ftype].append(feat_eval)

                        i = 0
                        # counting appended instances
                        inst_count = 0
                        for instance in feature:
                            inst_eval = 'NULL'
                            # XMLcase 1: feature instance contains evalue information (XML field inst_eval) 
                            if ('inst_eval' in instance.attrib):
                                #print tool + " instance evalue: "+ str(instance.attrib)
                                inst_eval = float(instance.attrib["inst_eval"])
                                start = int(instance.attrib["start"])
                                end = int(instance.attrib["end"])

                                if (inst_eval <= inst_efilter):
                                    logging.debug("append instance: " + ftype + ": " + str(instance.attrib))

                                    if single == 0 or single == 2:
                                        template_proteome[pID][ftype].append((inst_eval, start, end))
                                        inst_count += 1
                                        #print template_proteome[pID][ftype]

                                    else:
                                        query_proteome[pID][ftype].append((inst_eval, start, end))
                                        inst_count += 1
                                else:
                                    logging.debug("reject instance: "+ ftype + ": " + str(instance.attrib))

                            # XMLcase 2: feature instance contains NO evalue information (XML field inst_eval) 
                            else:
                                #NO instance based evalue information --> no instances can be rejected: set inst_count = 1
                                inst_count = 1
                                if (len(instance.attrib) == 2):
                                    #print tool +" no evalue: "+ str(instance.attrib)

                                    start = int(instance.attrib["start"])
                                    end = int(instance.attrib["end"])
                                    if single == 0 or single == 2:
                                        template_proteome[pID][ftype].append((inst_eval, start, end)) 
                                    else:
                                        query_proteome[pID][ftype].append((inst_eval, start, end))

                                else:
                                    if i == 0:
                                        start = int(instance.attrib["start"])
                                        i = 1
                                    else:
                                        end = int(instance.attrib["end"])
                                        if single == 0 or single == 2:
                                            template_proteome[pID][ftype].append((inst_eval, start, end)) 
                                        else:
                                            query_proteome[pID][ftype].append((inst_eval, start, end))
                                        i = 0
                        # any instance appended?
                        if (inst_count < 1):
                            #delete feature type
                            logging.info("Rejecting feature type " + str(ftype) + " due to rejection of all instances. Check for thresholds and E-values (instance based)")
                            if single == 0 or single == 2:
                                template_proteome[pID].pop(ftype)
                            else:
                                query_proteome[pID].pop(ftype)
                        if ftype not in clan_dict:
                            clan_dict[ftype] = fclan
    else:
        print("Error: " + path + " does not exist")
        quit()
    logging.debug("template_proteome: "+str(template_proteome))

# reads count.ref file
def ref_reader(path):
    global domain_count

    ref = open(path, "r+")
    reflines = ref.readlines()
    for line in reflines:
        tmp = line.rstrip('\n')
        tmp = tmp.split('\t')
        domain_count[tmp[0]] = int(tmp[1])

# writes count.ref file
def ref_writer(path):
    global domain_count

    ref = open(path, "w+")
    for domain in domain_count:
        ref.write(domain + "\t" + str(domain_count[domain]) + "\n")


##### start #####
if options.jobname.find("/")!=-1:
    outpath = options.jobname

else:
    outpath = expath + "/out/" + su_set_path(options.jobname)

    outpath = outpath + ".xml"
if weight_const == 1:
    constraints_in(options.weight_constraints)
featuretypes(options.featuretypes)
        
fc_main()

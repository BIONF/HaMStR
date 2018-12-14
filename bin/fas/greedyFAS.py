from operator import itemgetter
import xml.etree.ElementTree as ElTre
import logging
import inspect
import os
import math
from optparse import OptionParser
import time
from functools import partial
# from collections import defaultdict
from copy import deepcopy
import multiprocessing
# import sys

version = "1.6.1"
greedyfas_path = inspect.getfile(inspect.currentframe())
expath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

parser = OptionParser(description="You are running greedyFAS.py version " + str(version) + ".")
parser.add_option("-q", "--query", dest="query", default=expath + "/annotations_out/query",
                  help="Path to the folder containing the xml-files of the query protein set, default is fas/in/query")
parser.add_option("-s", "--seed", dest="seed", default=expath + "/annotations_out/template",
                  help="Path to the folder containing the xml-files of the template protein set, "
                       "default is fas/in/template")
parser.add_option("-w", "--score_weights", dest="score_weights", nargs=3, default=(0.7, 0.0, 0.3),
                  help="Defines how the three scores MS, CS and PS are weighted (MS, CS, PS), "
                       "the default is 0.7, 0.0, 0.3")
parser.add_option("-j", "--jobname", dest="jobname", default="out",
                  help="Defines the name for the outputfile, can also be used to define the output path, if no path is"
                       "given the output will be created in the FAS directory under out/")
parser.add_option("-r", "--ref_proteome", dest="ref_proteome", default=0,
                  help="Path to a reference proteome which can be used for the weighting of features, "
                       "by default there is no reference proteome used")
parser.add_option("-a", "--raw_output", dest="raw_output", default=0,
                  help="If set to 1, the FAS score will be printed to STDOUT. If 0, scores will be printed into output "
                       "file (XML format). If 2, both output variants are conducted.")
parser.add_option("--priority_check", dest="priority_check", default=1,
                  help="Set to 0 to deactivate priority checks prior to the calculation, skipping all priority checks "
                       "and going through all proteins exhaustively until time-limit is reached where it will switch "
                       "to priority mode, NOT RECOMMENDED as this may result in a very long runtime")
parser.add_option("-t", "--priority_threshold", dest="priority_threshold", default=50,
                  help="Change to define the feature number threshold for activating priority mode in the path "
                       "evaluation.")
parser.add_option("-m", "--max_cardinality", dest="max_cardinality", default=5000,
                  help="Change to define the threshold for the maximal cardinality of feature paths in a graph. "
                       "If max. cardinality is exceeded the priority mode will be used to for the path evaluation.")
parser.add_option("-l", "--log", dest="loglevel", default="ERROR",
                  help="Change the verbosity of the standard output to the screen. "
                       "Levels: DEBUG, INFO, WARNING, ERROR, CRITICAL")
parser.add_option("-f", "--efilter", dest="efilter", default="0.001",
                  help="E-value filter for hmm based search methods (feature based/complete sequence).")
parser.add_option("-i", "--inst_efilter", dest="inst_efilter", default="0.01",
                  help="E-value filter for hmm based search methods (instance based/complete sequence).")
parser.add_option("-g", "--weightcorrection", dest="weightcorrection", default="loge",
                  help="Function applied to the frequency of feature types during weighting, "
                       "options are linear(no function), loge(natural logarithm[Default]), log10(base-10 logarithm), "
                       "root4(4th root) and root8(8th root).")
parser.add_option("-x", "--weight_constraints", dest="weight_constraints", default=0,
                  help="Apply weight constraints via constraints file, by default there are no constraints.")
parser.add_option("-d", "--featuretypes", dest="featuretypes", default=expath + "/config/featuretypes",
                  help="inputfile that contains the tools/databases used to predict features")
parser.add_option("-e", "--extendedout", dest="extendedout", default=1,
                  help="0: only the scores and paths will be in the output, 1: "
                       "a file with the architecture for each protein will be created additionally")
parser.add_option("-y", "--feature_info", dest="feature_info", default=0,
                  help="if 1: a file with information on the abundance of all seed and query features in the reference")
parser.add_option("-b", "--bidirectional", dest="bidirectional", default=0,
                  help="if 1: calculate both scoring directions (separate files) on one core, if 2: calculate both "
                       "scoring directions (separate files) parallel on two cores")
parser.add_option("-c", "--max_overlap", dest="max_overlap", default=0,
                  help="maximum size overlap allowed, default is 0")
parser.add_option("--max_overlap_percentage", dest="max_overlap_percentage", default=0.4,
                  help="defines how much percent of a feature the overlap is allowed to cover, default is 0.4 (40%)")
parser.add_option("--classicMS", dest="classicMS", default=1,
                  help="(de)activate classic MS score")
parser.add_option("--timelimit", dest="timelimit", default=7200,
                  help="Sets a maximum time-limit in seconds a calculation between a pair of proteins is allowed to "
                       "take, default is 2 hours after which it will stop, set to 0 to deactivate; As FAS divides this "
                       "time among multiple processes, this limit does not necessarily represent the actual runtime,"
                       "especially if multiple cores are used")
parser.add_option("--cores", dest="cores", default=1,
                  help="Number of cores available for calculation, do not use together with --bidirectional 2, only "
                       "useful when not using priority_mode")

(options, args) = parser.parse_args()

### important vars ###             ###  var looks ###
# naming
# seed_proteome = {}           #{("protein_id", {("domain_name", [("START", "STOP")])})}
# query_proteome = {}          #{("protein_id", {("domain_name", [("START", "STOP")])})}
# clan_dict = {}               #{("domain_name", "clan")}
# search_features = {}         #{("F_0", ("domain_name", "POSITION", "Start", "Stop" ))}
# a_s_f = {}                   # additional seed features [non linearized]
# query_features = {}          #{("F_0", ("domain_name", "POSITION", "Start", "Stop" ))}
# a_q_f = {}                   # additional query features [non linearized]
# query_protein = {}           #{("domain_name", ["POSITION_1", "POSITION_2"])}
# query_clans = {}             #{("clan", "INSTANCES")}
# weights = {}                 #{("domain_name", "weight")}
# protein_lengths = {}         #{("protein_id", "LENGTH")}
# domain_count = {}            #{("domain", "COUNT")}

# hidden options
# tab separated table as output file #
taciturn = 1


# READ OPTIONS #
option_dict = {"p_path": options.seed, "s_path": options.query, "ref_proteome": options.ref_proteome, "weight_const": 0,
               "version": version}
loglevel = options.loglevel.upper()
try:
    option_dict["score_weights"] = []
    test = 0.0
    for weight in options.score_weights:
        option_dict["score_weights"].append(float(weight))
        test += float(weight)
    if test != 1.0:
        print("The sum of all values in -w [--score_weights] should be 1.0")
        quit()
    option_dict["score_weights"] = tuple(option_dict["score_weights"])
except ValueError:
    print(str(options.score_weights) + " is not a valid input for -w [--score_weights]")
    quit()
try:
    option_dict["cores"] = int(options.cores)
    if option_dict["cores"] < 1:
        print(str(options.cores) + " is not a valid input for [--cores], must be 1 or more")
        quit()
except ValueError:
    print(str(options.cores) + " is not a valid input for [--cores], must be 1 or more")
    quit()
try:
    option_dict["classicMS"] = int(options.classicMS)
    if not (option_dict["classicMS"] == 1 or option_dict["classicMS"] == 0):
        print(str(options.classicMS) + " is not a valid input for [--classicMS], must be 0 or 1")
        quit()
except ValueError:
    print(str(options.classicMS) + " is not a valid input for [--classicMS], must be 0 or 1")
    quit()
try:
    option_dict["timelimit"] = int(options.timelimit)
    if not (option_dict["timelimit"] >= 0):
        print(str(options.timelimit) + " is not a valid input for [--timelimit], must be integer higher than 0 or 0 to "
                                       "deactivate")
        quit()
except ValueError:
    print(str(options.timelimit) + " is not a valid input for [--timelimit], must be integer higher than 0 or 0 to "
                                   "deactivate")
    quit()
try:
    option_dict["priority_mode"] = int(options.priority_check)
    if not (option_dict["priority_mode"] == 1 or option_dict["priority_mode"] == 0):
        print(str(options.priority_check) + " is not a valid input for [--priority_check], must be 0 or 1")
        quit()
except ValueError:
    print(str(options.priority_check) + " is not a valid input for [--priority_check], must be 0 or 1")
    quit()
try:
    option_dict["priority_threshold"] = int(options.priority_threshold)
except ValueError:
    print(str(options.priority_threshold) + " is not a valid input for -c [--priority_threshold]")
    quit()
try:
    option_dict["max_cardinality"] = int(options.max_cardinality)
except ValueError:
    print(str(options.max_cardinality) + " is not a valid input for -m [--max_cardinality]")
    quit()
try:
    option_dict["max_overlap"] = int(options.max_overlap)
except ValueError:
    print(str(options.max_overlap) + " is not a valid input for -c [--max_overlap]")
    quit()
try:
    option_dict["efilter"] = float(options.efilter)
except ValueError:
    print(str(options.efilter) + " is not a valid input for -f [--efilter]")
    quit()
try:
    option_dict["inst_efilter"] = float(options.inst_efilter)
except ValueError:
    print(str(options.inst_efilter) + " is not a valid input for -i [--inst_efilter]")
    quit()
if option_dict["ref_proteome"] == 0:
    option_dict["MS_uni"] = 1
    option_dict["weight_correction"] = 0
else:
    option_dict["MS_uni"] = 0
    if options.weightcorrection == "log10":
        option_dict["weight_correction"] = "log10"
    elif options.weightcorrection == "loge":
        option_dict["weight_correction"] = "loge"
    elif options.weightcorrection == "root4":
        option_dict["weight_correction"] = "root4"
    elif options.weightcorrection == "root8":
        option_dict["weight_correction"] = "root8"
    elif options.weightcorrection == "linear":
        option_dict["weight_correction"] = 0
    else:
        print(str(options.weightcorrection) + "is not a valid input for -g [--weightcorrection]")
        quit()
try:
    if int(options.raw_output) == 0:
        option_dict["output"] = 0
    elif int(options.raw_output) == 1:
        option_dict["output"] = 1
    elif int(options.raw_output) == 2:
        option_dict["output"] = 2
    else:
        print(str(options.raw_output) + "is not a valid input for -a [--raw_output]")
        quit()
except ValueError:
    print(str(options.raw_output) + "is not a valid input for -a [--raw_output]")
    quit()
try:
    if int(options.extendedout) == 0:
        option_dict["e_output"] = 0
    elif int(options.extendedout) == 1:
        option_dict["e_output"] = 1
    else:
        print(str(options.extendedout) + " is not a valid input for -e [--extendedout]")
        quit()
except ValueError:
    print(str(options.extendedout) + " is not a valid input for -e [--extendedout]")
    quit()
try:
    if int(options.feature_info) == 0:
        option_dict["feature_info"] = 0
    elif int(options.feature_info) == 1:
        option_dict["feature_info"] = 1
    else:
        print(str(options.feature_info) + " is not a valid input for -y [--feature_info]")
        quit()
except ValueError:
    print(str(options.feature_info) + " is not a valid input for -y [--feature_info]")
    quit()

if options.ref_proteome != 0 and options.weight_constraints != 0:
    option_dict["weight_const"] = 1
elif options.weight_constraints != 0:
    print("[--weight_constraints] only works with a reference proteome")
    quit()
try:
    if int(options.bidirectional) == 2:

        option_dict["bidirectional"] = 2
    elif int(options.bidirectional) == 1:
        option_dict["bidirectional"] = 1
    elif int(options.bidirectional) == 0:
        option_dict["bidirectional"] = 0
    else:
        print(str(options.bidirectional) + " is not a valid input for -b [--bidirectional]")
        quit()
except ValueError:
    print(str(options.bidirectional) + " is not a valid input for -b [--bidirectional]")
    quit()
try:
    if 0.0 <= float(options.max_overlap_percentage) <= 1.0:
        option_dict["max_overlap_percentage"] = float(options.max_overlap_percentage)
    else:
        print("[--max_overlap_percentage] should be between 0.0 and 1.0")
        quit()
except ValueError:
    print(str(options.max_overlap_percentagea) + " is not a valid input for [--max_overlap_percentage]")
    quit()

### SETUP LOGGING OPTIONS ###

# possible logging configurations
# logging into stdout with time stamp:
logging.basicConfig(level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')
# logging into file with line number:
# logging.basicConfig(filename='testlog.log', filemode='w', level=loglevel,
#                    format='%(lineno)s - %(levelname)s - %(message)s')
# logging into stdout with line number:
# logging.basicConfig(level=loglevel, format='%(lineno)s - %(levelname)s - %(message)s')

logging.info(
    'greedyFAS.py started with options: priority_threshold=' + str(option_dict["priority_threshold"]) + ', log_level='
    + str(loglevel))
logging.info(
    'score_weights are set to: ' + str(option_dict["score_weights"][0]) + " " + str(option_dict["score_weights"][1])
    + " " + str(option_dict["score_weights"][2]))
logging.info('ref_proteome is set to: ' + str(option_dict["ref_proteome"]))


########## Flow Control ########## <fc>
# flow control functions

def fc_start(option):
    """Overhead function,
    this function manages the individual functions that read the input files and prepares the data for the main script.
    Function calls: xmlreader(), w_count_ref(), w_count(), fc_main()

    :param option: dictionary that contains the main option variables of FAS
    """
    clan_dict = {}
    seed_proteome = {}
    query_proteome = {}
    protein_lengths = {}
    domain_count = {}
    prot_count = 0

    ## MS_uni set to 0 when no weighting is conducted
    if option["MS_uni"] == 0:
        for ftype in option["input_linearized"]:
            seed_proteome, protein_lengths, clan_dict = xmlreader(option["ref_proteome"] + "/" + ftype + ".xml", 2,
                                                                  ftype, True, seed_proteome, protein_lengths,
                                                                  clan_dict, option)
        for ftype in option["input_normal"]:
            seed_proteome, protein_lengths, clan_dict = xmlreader(option["ref_proteome"] + "/" + ftype + ".xml", 2,
                                                                  ftype, True, seed_proteome, protein_lengths,
                                                                  clan_dict, option)

        prot_count, domain_count = w_count_ref(seed_proteome)
        seed_proteome = {}
        protein_lengths = {}
    for ftype in option["input_linearized"]:
        seed_proteome, protein_lengths, clan_dict = xmlreader(option["p_path"] + "/" + ftype + ".xml", 0, ftype, True,
                                                              seed_proteome, protein_lengths, clan_dict, option)
        query_proteome, protein_lengths, clan_dict = xmlreader(option["s_path"] + "/" + ftype + ".xml", 1, ftype, True,
                                                               query_proteome, protein_lengths, clan_dict, option)
    for ftype in option["input_normal"]:
        seed_proteome, protein_lengths, clan_dict = xmlreader(option["p_path"] + "/" + ftype + ".xml", 0, ftype, False,
                                                              seed_proteome, protein_lengths, clan_dict, option)
        query_proteome, protein_lengths, clan_dict = xmlreader(option["s_path"] + "/" + ftype + ".xml", 1, ftype, False,
                                                               query_proteome, protein_lengths, clan_dict, option)
    if option["MS_uni"] == 0:
        relevant_features, domain_count = w_count(prot_count, domain_count, seed_proteome, query_proteome)
    else:
        relevant_features = {}
    if option["weight_correction"] != 0:
        domain_count = w_weight_correction(option["weight_correction"],
                                           domain_count)  # use correction function on counts

    if option["bidirectional"] == 2:
        p1 = multiprocessing.Process(target=fc_main, args=(relevant_features, prot_count, domain_count, seed_proteome,
                                                           query_proteome, protein_lengths, clan_dict, option))
        p1.start()
        tmp = {}
        for protein in protein_lengths:
            if protein[0:4] == "seed":
                tmp["query_" + protein[5:]] = protein_lengths[protein]
            else:
                tmp["seed_" + protein[6:]] = protein_lengths[protein]
        option["e_output"] = 0
        org_outpath = option["outpath"]
        option["outpath"] += "_reverse"
        option["feature_info"] = 0
        p2 = multiprocessing.Process(target=fc_main, args=(relevant_features, prot_count, domain_count, query_proteome,
                                                           seed_proteome, tmp, clan_dict, option))
        p2.start()
        p1.join()
        p2.join()
        bidirectionout(org_outpath)
    elif option["bidirectional"] == 1:
        fc_main(relevant_features, prot_count, domain_count, seed_proteome, query_proteome, protein_lengths, clan_dict,
                option)
        tmp = {}
        for protein in protein_lengths:
            if protein[0:4] == "seed":
                tmp["query_" + protein[5:]] = protein_lengths[protein]
            else:
                tmp["seed_" + protein[6:]] = protein_lengths[protein]
        option["e_output"] = 0
        org_outpath = option["outpath"]
        option["outpath"] += "_reverse"
        option["feature_info"] = 0
        fc_main(relevant_features, prot_count, domain_count, query_proteome, seed_proteome, tmp, clan_dict, option)
        bidirectionout(org_outpath)
    elif option["bidirectional"] == 0:
        fc_main(relevant_features, prot_count, domain_count, seed_proteome, query_proteome, protein_lengths, clan_dict,
                option)


def fc_main(relevant_features, prot_count, domain_count, seed_proteome, query_proteome, protein_lengths, clan_dict,
            option):
    """Main function,
    manages linearization and scoring, creates the output
    Function calls: w_weighting_constraints(), w_weighting(), su_lin_query_protein(), pb_graphtraversal(),
                    su_query_protein(), su_search_protein(), pb_entire_graphtraversal_priority(),
                    sf_entire_calc_score(), pb_entire_main_nongreedy(), w_weight_const_rescale()

    :param relevant_features: dictionary that contains the rarity of all features from seed and query
    :param prot_count: integer that stores the number of proteins in the reference
    :param domain_count: dictionary that contains the (adjusted) counts of each feature in the reference proteome
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param query_proteome: dictionary that contains the feature architecture of all query proteins
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param clan_dict: dictionary that maps features to clans
    :param option: dictionary that contains the main option variables of FAS
    """
    logging.info("fc_main")

    # score model M2 - used for iteration/incremental evaluation

    # entire mode
    # both proteins (search and query) will be linearized
    mode = {}
    turn = 0
    weights = {}
    adjusted_weights = {}
    weight_tmp = {}
    if option["output"] == 0 or option["output"] == 2:
        if option["e_output"] == 1:
            a_out = open(option["outpath"] + "_architecture.xml", "w+")
            a_out.write("<?xml version=\"1.0\"?>\n")
            a_out.write("<architectures FAS_version=\"" + str(option["version"]) + "\">\n")
        out = open(option["outpath"] + ".xml", "w+")
        out.write("<?xml version=\"1.0\"?>\n")
        settings_out = {"priority": "off", "weighting": "uniform", "constraints": "none", "overlap":
                        str(option["max_overlap"]) + "/" + str(option["max_overlap_percentage"]), "filters":
                        str(option["efilter"]) + "/" + str(option["inst_efilter"]), "lin": "", "norm": "", "ms":
                        "classic"}

        if option["classicMS"] == 0:
            settings_out["ms"] = "new"
        if option["timelimit"] == 0 or option["priority_mode"] == 1:
            settings_out["time"] = "off"
        else:
            settings_out["time"] = str(option["timelimit"]) + "s"
        for tool in option["input_linearized"]:
            settings_out["lin"] = settings_out["lin"] + "/" + tool
        settings_out["lin"] = settings_out["lin"].lstrip("/")
        if len(settings_out["lin"]) == 0:
            settings_out["lin"] = "none"
        for tool in option["input_normal"]:
            settings_out["norm"] = settings_out["norm"] + "/" + tool
        settings_out["norm"] = settings_out["norm"].lstrip("/")
        if len(settings_out["norm"]) == 0:
            settings_out["norm"] = "none"
        if option_dict["weight_const"] == 1:
            settings_out["constraints"] = option["constname"]
        if option["priority_mode"] == 1:
            settings_out["priority"] = str(option["priority_threshold"]) + "/" + str(option["max_cardinality"])
        if option["MS_uni"] == 0 and option["weight_correction"] != 0:
            settings_out["weighting"] = option["ref_proteome"].rstrip("/").split("/")[-1] + "/" + \
                                        option["weight_correction"]
        elif option["MS_uni"] == 0:
            settings_out["weighting"] = option["ref_proteome"].rstrip("/").split("/")[-1] + "/none"
        out.write("<out FAS_version=\"" + str(option["version"]) + "\" weighting=\"" + settings_out["weighting"] +
                  "\" constraints=\"" + settings_out["constraints"] + "\" MS=\"" + settings_out["ms"] + "\" priority=\""
                  + settings_out["priority"] + "\" overlap=\"" + settings_out["overlap"] + "\" efilters=\"" +
                  settings_out["filters"] + "\" timelimit=\"" + settings_out["time"] + "\" scoreweights=\"" +
                  str(option["score_weights"]) + "\" cores=\"" + str(option["cores"]) + "\" linearized=\"" +
                  settings_out["lin"] + "\" normal=\"" + settings_out["norm"] + "\">\n")
    for query in query_proteome:
        go_priority = False
        if option["MS_uni"] == 0:
            if option["weight_const"]:
                weights, domain_count = w_weighting_constraints(query, domain_count, query_proteome, option)
            else:
                weights, domain_count = w_weighting(query, domain_count, query_proteome)
        lin_query_set, query_features, a_q_f, query_clans, clan_dict = su_lin_query_protein(query, query_proteome,
                                                                                            protein_lengths, clan_dict)
        tmp_query_graph = pb_region_paths_nongreedy(pb_region_mapper(lin_query_set, query_features,
                                                                     option["max_overlap"],
                                                                     option["max_overlap_percentage"]))
        # PRIORITY CHECK 1: checking for number of instances - assess complexity of the feature graph
        if int(len(query_features)) > int(option["priority_threshold"]) and option["priority_mode"] == 1:
            go_priority = True
            all_query_paths = "PRIORITY"
        elif option["priority_mode"] == 1:
            # creating all paths
            all_query_paths = pb_graphtraversal(tmp_query_graph, [], [], option)
            # PRIORITY CHECK 2: checking for number of paths in the feature graph
            if int(len(all_query_paths)) > int(option["max_cardinality"]):
                logging.info("Switched to priority mode due to number of paths.")
                go_priority = True
        elif len(tmp_query_graph) == 0:
            all_query_paths = []
        else:
            all_query_paths = "NOPRIORITY"
        if option["output"] == 0 or option["output"] == 2:
            out.write("\t<query id=\"" + query + "\" length=\"" + str(int(protein_lengths["query_" + query])) + "\">\n")
        elif taciturn == 1:
            out = open(option["outpath"], "w+")
        query_protein, query_clans, clan_dict = su_query_protein(query, query_proteome, protein_lengths, clan_dict)
        for protein in seed_proteome:
            go_priority_2 = False
            calcstart = time.time()
            timeover = 0
            mode[protein] = 0
            pathcount = 0
            if option["MS_uni"] == 0:
                if option["weight_const"]:
                    weights, domain_count = w_weighting_constraints(protein, domain_count, seed_proteome, option)
                else:
                    weights, domain_count = w_weighting(protein, domain_count, seed_proteome)
            search_protein, search_features, a_s_f = su_search_protein(protein, seed_proteome, protein_lengths)
            search_protein = tuple(search_protein)
            max_fixture = ([], (0.0, 0.0, 0.0, 0.0, 0, 0, False), [], protein)

            # check for available paths
            # query <--VS-- seed
            # case M2.1: empty(query)

            if int(len(all_query_paths)) == 0:
                logging.warning("CASE M2.1: No paths (pfam or smart annotations) in query.")
                # case M2.1.1: empty(query)-empty(search)
                # should be the best fix independent from weight
                if int(len(search_features)) == 0:
                    logging.warning("CASE M2.1.1: empty vs empty.")
                    path = list(a_s_f.keys())
                    query_architecture = list(a_q_f.keys())
                    score_w = sf_entire_calc_score(path, query_architecture, weights, search_features, a_s_f,
                                                   query_features, a_q_f, clan_dict, option)
                    mode[protein] = 2
                else:
                    # case M2.1.2: empty(query)-graph(search)
                    logging.warning("CASE M2.1.2: empty vs graph.")
                    tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, [], search_features, weights,
                                                              query_features, seed_proteome, a_s_f, a_q_f, clan_dict,
                                                              query_clans, protein_lengths, "OFF", option)
                    path = tmp_path_score[0][0]
                    score_w = tmp_path_score[0][1]
                    query_architecture = tmp_path_score[0][2]
                    mode[protein] = 2

                # set max_fixture according to
                max_fixture = (path, score_w, query_architecture, protein)

            # handle all paths
            # case M2.2: graph(query)
            elif option["priority_mode"] == 0:
                mode[protein] = 0
                stack, jobpaths = pb_create_jobs(tmp_query_graph, option)
                jobpool = multiprocessing.Pool(processes=option["cores"])
                timelimit = 0.0
                if len(stack) > 0:
                    if option["timelimit"] > 0:
                        timelimit = option["timelimit"] / len(stack)
                        func = partial(pb_graph_traversal_sub, search_protein, protein, search_features, weights,
                                       query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                       protein_lengths, timelimit, tmp_query_graph, option)
                    else:
                        func = partial(pb_graph_traversal_sub, search_protein, protein, search_features, weights,
                                       query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                       protein_lengths, 0, tmp_query_graph, option)
                    pool_results = jobpool.map_async(func, stack)
                    jobpool.close()
                    jobpool.join()
                    pool_results.wait()
                    for pool_result in pool_results.get():
                        if pool_result[1] > timelimit and option["timelimit"] > 0:
                            timeover += 1
                            go_priority_2 = True
                        if pool_result[0][1][5] >= max_fixture[1][3]:
                            max_fixture = deepcopy(pool_result[0])
                timecheck = time.time()
                if len(jobpaths) > 0:
                    # case M2.2.1: graph(query)-empty(search)
                    if int(len(search_features)) == 0:
                        for jobpath in jobpaths:
                            logging.warning("CASE M2.2: graph.")
                            pathcount += 1
                            logging.warning("CASE M2.2.1: graph vs empty.")
                            # special case: protein with no pfam or smart domains
                            # get score for a_s_f and query_path directly
                            path = list(a_s_f.keys())
                            query_path_ad = jobpath + list(a_q_f.keys())
                            score_w = sf_entire_calc_score(path, query_path_ad, weights, search_features, a_s_f,
                                                           query_features, a_q_f, clan_dict, option)

                            # check for max scoring fixture of path and query_path
                            if score_w[5]:
                                if score_w[3] >= max_fixture[1][3]:
                                    max_fixture = (path, score_w, query_path_ad, protein)
                                    logging.warning("max_fixture: " + str(max_fixture))
                            else:
                                if (not max_fixture[1][5]) and (score_w[4] >= max_fixture[1][4]):
                                    max_fixture = (path, score_w, query_path_ad, protein)
                    else:
                        # case M2.2.2 graph(query)-graph(search)
                        # regular traversal of graph based on search_protein
                        logging.warning("CASE M2.2.2: graph vs graph.")
                        jobs = []
                        jobcount = (len(jobpaths) / option["cores"], len(jobpaths) % option["cores"])
                        if jobcount[0] == 0:
                            timelimit = (option["timelimit"] - (timecheck - calcstart)) / float(jobcount[1])
                            for i in range(jobcount[1]):
                                jobs.append([jobpaths[-(i + 1)]])
                        else:
                            timelimit = (option["timelimit"] - (timecheck - calcstart)) / float(option["cores"])
                            counter = 0
                            for i in range(option["cores"]):
                                jobs.append(jobpaths[counter: (i + 1) * jobcount[0]])
                                counter = (i + 1) * jobcount[0]
                            for i in range(jobcount[1]):
                                jobs[i].append(jobpaths[-(i + 1)])
                        func = partial(pb_calc_sub, search_protein, protein, search_features, weights, query_features,
                                       seed_proteome, a_s_f, a_q_f, clan_dict, query_clans, protein_lengths, timelimit,
                                       option)
                        pool_results = jobpool.map_async(func, jobs)
                        jobpool.close()
                        jobpool.join()
                        pool_results.wait()
                        for pool_result in pool_results.get():
                            if pool_result[1] > timelimit:
                                timeover += 1
                                go_priority_2 = True
                            if pool_result[0][1][5] >= max_fixture[1][3]:
                                max_fixture = deepcopy(pool_result[0])
            if go_priority or go_priority_2:
                mode[protein] = 1
                all_query_paths = []
                priority_list = []
                for domain_type in query_proteome[query]:
                    if query_proteome[query][domain_type][0]:
                        priority_list.append(domain_type)
                priority_list.append("NONE")
                for domain_type in priority_list:
                    all_query_paths += (
                        pb_entire_graphtraversal_priority(tmp_query_graph, domain_type, protein, 1, search_features,
                                                          weights, query_features, seed_proteome, a_s_f, a_q_f,
                                                          clan_dict, query_clans, protein_lengths, option))
                    logging.info("domain_type: " + str(domain_type))
                    logging.debug("all_query_paths: " + str(all_query_paths))
                    # max fixture of (search_path, score, query_path)
            if option["priority_mode"] == 1 or go_priority_2:
                for query_path in all_query_paths:
                    logging.warning("CASE M2.2: graph.")
                    pathcount += 1

                    # case M2.2.1: graph(query)-empty(search)
                    if int(len(search_features)) == 0:
                        logging.warning("CASE M2.2.1: graph vs empty.")
                        # special case: protein with no pfam or smart domains
                        # get score for a_s_f and query_path directly
                        path = list(a_s_f.keys())
                        query_path_ad = query_path + list(a_q_f.keys())
                        score_w = sf_entire_calc_score(path, query_path_ad, weights, search_features, a_s_f,
                                                       query_features, a_q_f, clan_dict, option)
                    else:
                        # case M2.2.2 graph(query)-graph(search)
                        # regular traversal of graph based on search_protein
                        logging.warning("CASE M2.2.2: graph vs graph.")
                        if option["priority_mode"] == 1:
                            tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, query_path,
                                                                      search_features, weights, query_features,
                                                                      seed_proteome, a_s_f, a_q_f, clan_dict,
                                                                      query_clans, protein_lengths, "OFF", option)
                        else:
                            tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, query_path,
                                                                      search_features, weights, query_features,
                                                                      seed_proteome, a_s_f, a_q_f, clan_dict,
                                                                      query_clans, protein_lengths, "OVER", option)
                        path = tmp_path_score[0][0]
                        score_w = tmp_path_score[0][1]
                        query_path_ad = tmp_path_score[0][2]
                        mode[protein] = tmp_path_score[1]
                        logging.debug("tmp_path_score " + str(tmp_path_score))  # search path, score, mode

                    # check for max scoring fixture of path and query_path
                    if score_w[5]:
                        if score_w[3] >= max_fixture[1][3]:
                            max_fixture = (path, score_w, query_path_ad, protein)
                            logging.warning("max_fixture: " + str(max_fixture))
                    else:
                        if (not max_fixture[1][5]) and (score_w[4] >= max_fixture[1][4]):
                            max_fixture = (path, score_w, query_path_ad, protein)

                logging.info("Found: " + str(pathcount) + " path(s) for query.")
                logging.debug("Path fixture: " + str(max_fixture))

            timecheck = time.time()
            score = max_fixture[1]
            if option["output"] == 0 or option["output"] == 2:
                try:
                    if mode[protein] == 2:
                        mode_out = "greedy"
                    elif mode[protein] == 1:
                        mode_out = "priority"
                    else:
                        mode_out = "exhaustive"
                except KeyError:
                    mode_out = "default"
                if go_priority or go_priority_2:
                    mode_out += "/priority"
                else:
                    mode_out += "/exhaustive"
                runtime = str(round(timecheck - calcstart, 4)) + "s"
                if timeover >= 1:
                    runtime += "/exceeded_timelimit"
                out.write("\t\t<template id=\"" + protein + "\" score=\"" + str(score[3]) + "\" MS=\"" + str(
                    score[0]) + "\" PS=\"" + str(score[1]) + "\" CS=\"" + str(score[2]) + "\" LS=\"" + str(
                    score[6]) + "\" length=\"" + str(int(protein_lengths["seed_" + str(protein)])) + "\" mode=\"" +
                          mode_out + "\" calculationTime=\"" + runtime + "\" >\n")
                if option["e_output"] == 1 and turn == 0:
                    a_out.write("\t<template id=\"" + protein + "\" length=\"" + str(
                        int(protein_lengths["seed_" + str(protein)])) + "\">\n")
                    a_out.write("\t\t<architecture>\n")

                    for feature in seed_proteome[protein]:
                        if option["MS_uni"] == 0:
                            a_out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\"" + str(
                                seed_proteome[protein][feature][1]) + "\" weight=\"" + str(
                                weights[feature]) + "\">\n")
                        else:
                            a_out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\"" + str(
                                seed_proteome[protein][feature][1]) + "\">\n")
                        for instance in seed_proteome[protein][feature][2:]:
                            a_out.write("\t\t\t\t<instance inst_eval=\"" + str(instance[0]) + "\" start=\"" + str(
                                instance[1]) + "\" end=\"" + str(instance[2]) + "\"/>\n")
                        a_out.write("\t\t\t</feature>\n")
                    a_out.write("\t\t</architecture>\n")
                    a_out.write("\t</template>\n")

            best_template_path = []
            best_query_path = []

            for feature in max_fixture[0]:
                if feature in search_features:
                    logging.debug(str(search_features[feature][0]) + " " + str(search_features[feature][1]) + " " + str(
                        search_features[feature][2]) + " " + str(search_features[feature][3]))
                    best_template_path.append((search_features[feature][0], search_features[feature][1],
                                               search_features[feature][2], search_features[feature][3]))
                else:
                    best_template_path.append(
                        (a_s_f[feature][0], a_s_f[feature][1], a_s_f[feature][2], a_s_f[feature][3]))

            for feature in max_fixture[2]:
                if feature in query_features:
                    best_query_path.append((query_features[feature][0], query_features[feature][1],
                                            query_features[feature][2], query_features[feature][3]))
                else:
                    best_query_path.append((a_q_f[feature][0], a_q_f[feature][1], a_q_f[feature][2], a_q_f[feature][3]))
            if option["output"] == 0 or option["output"] == 2:
                path_tmp = {}
                path_tmp_query = {}
                scale = 0
                if option["weight_const"] == 1:
                    path_tmp2 = []
                    for feature in best_template_path:
                        if feature[0] not in path_tmp2:
                            path_tmp2.append(feature[0])
                    adjusted_weights = w_weight_const_rescale(path_tmp2, weights, search_features, True, option)
                    weight_tmp = {}
                    for adj_feature in adjusted_weights:
                        weight_tmp[adj_feature] = weights[adj_feature]
                        weights[adj_feature] = adjusted_weights[adj_feature]
                for feature in best_template_path:
                    if feature[0] in path_tmp:
                        path_tmp[feature[0]].append((feature[2], feature[3]))
                    else:
                        path_tmp[feature[0]] = [(feature[2], feature[3])]
                        if option["MS_uni"] == 0:
                            scale += weights[feature[0]]
                for feature in best_query_path:
                    if feature[0] in path_tmp_query:
                        path_tmp_query[feature[0]].append((feature[2], feature[3]))
                    else:
                        path_tmp_query[feature[0]] = [(feature[2], feature[3])]
                logging.debug("path_tmp: " + str(path_tmp))

                ## unweighted case
                if option["MS_uni"] == 0:
                    if scale > 0:
                        scale = 1.0 / float(scale)
                    else:
                        scale = 1.0
                ## print path
                out.write("\t\t\t<template_path>\n")
                for feature in path_tmp:
                    if option["MS_uni"] == 0:
                        out.write("\t\t\t\t<feature type=\"" + feature + "\" corrected_weight=\"" + str(
                            weights[feature] * scale) + "\">\n")
                    else:
                        out.write("\t\t\t\t<feature type=\"" + feature + "\">\n")
                    for tmp_inst in path_tmp[feature]:
                        out.write("\t\t\t\t\t<instance start=\"" + str(tmp_inst[0]) + "\" end=\"" + str(
                            tmp_inst[1]) + "\"/>\n")
                    out.write("\t\t\t\t</feature>\n")
                out.write("\t\t\t</template_path>\n")
                out.write("\t\t\t<query_path>\n")
                for feature in path_tmp_query:
                    out.write("\t\t\t\t<feature type=\"" + feature + "\">\n")
                    for tmp_inst in path_tmp_query[feature]:
                        out.write("\t\t\t\t\t<instance start=\"" + str(tmp_inst[0]) + "\" end=\"" + str(
                            tmp_inst[1]) + "\"/>\n")
                    out.write("\t\t\t\t</feature>\n")
                out.write("\t\t\t</query_path>\n")
                out.write("\t\t</template>\n")
                if option["weight_const"] == 1:
                    for adj_feature in adjusted_weights:
                        weights[adj_feature] = weight_tmp[adj_feature]
            if option["output"] == 1 or option["output"] == 2:
                print (score[3])
                ## hidden ##
                if taciturn == 1 and option["output"] != 2:
                    out.write(protein + "\t" + str(score[3]) + "\n")
        if option["output"] == 0 or option["output"] == 2:
            out.write("\t</query>\n")
        turn = 1
    if option["output"] == 0 or option["output"] == 2:
        out.write("</out>")
        out.close()
        if option["e_output"] == 1:
            for query in query_proteome:
                a_out.write(
                    "\t<query id=\"" + query + "\" length=\"" + str(int(protein_lengths["query_" + query])) + "\">\n")
                a_out.write("\t\t<architecture>\n")
                for feature in query_proteome[query]:
                    a_out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\"" + str(
                        query_proteome[query][feature][1]) + "\">\n")
                    for instance in query_proteome[query][feature][2:]:
                        a_out.write("\t\t\t\t<instance inst_eval=\"" + str(instance[0]) + "\" start=\"" + str(
                            instance[1]) + "\" end=\"" + str(instance[2]) + "\"/>\n")

                    a_out.write("\t\t\t</feature>\n")
                a_out.write("\t\t</architecture>\n")
                a_out.write("\t</query>\n")
            a_out.write("</architectures>")
            a_out.close()
        if option["feature_info"] == 1:
            f_out = open(option["outpath"] + "_features", "w+")
            f_out.write("#proteins: " + str(prot_count) + "\n")
            tmp = []
            for feature in relevant_features:
                tmp.append((feature, relevant_features[feature]))
            tmp = sorted(tmp, key=lambda x: x[1])
            for feature in tmp:
                f_out.write(feature[0] + "\t" + str(feature[1]) + "\n")
            f_out.close()


########## Start Up Functions ########## <su>

def su_query_protein(protein_id, query_proteome, protein_lengths, clan_dict):
    """Initializes variables for the current query protein

    :param protein_id: String that contains the identifier of the query protein
    :param query_proteome: dictionary that contains the feature architecture of all query proteins
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param clan_dict: dictionary that maps features to clans
    :return: query_protein, query_clans, clan_dict
    """
    query_protein = {}
    query_clans = {}

    for feature in query_proteome[protein_id]:
        query_protein[feature] = []
        clan = clan_dict[feature]
        if clan in query_clans:
            query_clans[clan] += len(query_proteome[protein_id][feature]) - 2
        else:
            query_clans[clan] = len(query_proteome[protein_id][feature]) - 2
        for instance in query_proteome[protein_id][feature][2:]:
            position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(
                protein_lengths["query_" + str(protein_id)])
            position = round(position, 8)
            query_protein[feature].append(position)
    return query_protein, query_clans, clan_dict


def su_lin_query_protein(protein_id, query_proteome, protein_lengths, clan_dict):
    """Initializes variables for the current query protein (linerization version)

    :param protein_id: String that contains the identifier of the query protein
    :param query_proteome: dictionary that contains the feature architecture of all query proteins
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param clan_dict: dictionary that maps features to clans
    :return:lin_query_protein, query_features, a_q_f(additional[not linearized] query features), query_clans, clan_dict
    """
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
            query_clans[clan] += len(query_proteome[protein_id][feature]) - 2
        else:
            query_clans[clan] = len(query_proteome[protein_id][feature]) - 2

        if query_proteome[protein_id][feature][0]:
            for instance in query_proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(
                    protein_lengths["query_" + str(protein_id)])
                position = round(position, 8)
                key = "F_" + str(i)
                query_features[key] = (feature, position, instance[1], instance[2])
                tmp.append((key, instance[1]))
                i += 1
        else:
            for instance in query_proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(
                    protein_lengths["query_" + str(protein_id)])
                position = round(position, 8)
                key = "O_" + str(i)
                a_q_f[key] = (feature, position, instance[1], instance[2])
                i += 1

    # sort  instances
    tmp2 = sorted(tmp, key=itemgetter(1))
    for x in tmp2:
        lin_query_protein.append(x[0])

    return lin_query_protein, query_features, a_q_f, query_clans, clan_dict


def su_search_protein(protein_id, seed_proteome, protein_lengths):
    """Initializes variables for the current seed protein

    :param protein_id: String that contains the identifier of the seed protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :return: search_protein, search_features, a_s_f (additional [not linearized] seed features)
    """
    search_features = {}
    search_protein = []
    a_s_f = {}
    tmp = []
    i = 0
    for feature in seed_proteome[protein_id]:
        # True if feature is going to be linearized
        if seed_proteome[protein_id][feature][0]:
            for instance in seed_proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(
                    protein_lengths["seed_" + str(protein_id)])
                position = round(position, 8)
                key = "F_" + str(i)
                search_features[key] = (feature, position, instance[1], instance[2])
                tmp.append((key, instance[1]))
                i += 1
        else:
            for instance in seed_proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(
                    protein_lengths["seed_" + str(protein_id)])
                position = round(position, 8)
                key = "O_" + str(i)
                a_s_f[key] = (feature, position, instance[1], instance[2])
                i += 1

    # sort  instances
    tmp2 = sorted(tmp, key=itemgetter(1))
    for x in tmp2:
        search_protein.append(x[0])
    return tuple(search_protein), search_features, a_s_f


def su_set_path(jobname, path):
    """Checks default output path if only a name for job is given

    :param jobname: String contains name of output file
    :param path: String contains output path
    :return: jobname
    """
    if not os.path.exists(path + "/out/"):
        os.makedirs(path + "/out/")
    if os.path.exists(path + "/out/" + jobname + ".xml"):
        i = 1
        while os.path.exists(path + "/out/" + jobname + "_" + str(i) + ".xml"):
            i += 1
        jobname = jobname + "_" + str(i)
    return jobname


# Path-building Functions # <pb>
# Used for the Pfam/Smart domains (default)
def pb_create_jobs(graph, option):
    """Does a breadth-first search to divide the feature architecture graph into several sub-problems so that it can be
    worked on by multiple cores. There are more jobs created than there are cores to counter the possibility of a core
    having nothing to do because of different running times for each job.

    :param graph: graph containing all paths through a protein architecture
    :param option: dictionary that contains the main option variables of FAS
    :return: queue (contains the starting points and sub-paths of the jobs), paths (empty unless b-f search reaches end
             node)
    """
    limit = 4 * option["cores"]
    queue = [("START", [])]
    paths = []
    while len(queue) < limit and queue:
        vertex, path = queue.pop(0)
        for next_vertex in graph[vertex]:
            if next_vertex == "END":
                paths.append(path)
            else:
                queue.append((next_vertex, path + [next_vertex]))
    return queue, paths


def pb_graph_traversal_sub(search_protein, protein, search_features, weights, query_features, seed_proteome, a_s_f,
                           a_q_f, clan_dict, query_clans, protein_lengths, timelimit, query_graph, option, stack):
    """Sub function for multi-processing: This function does a depth-first search in the feature graph. The starting
     point in the graph is given in the stack.

    :param search_protein: contains feature architecture of the current seed protein [dictionary]
    :param protein: String that contains the identifier of the seed protein
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param timelimit: the amount of time this job is allowed to take before the calculation is stopped and priority
                      mode is used instead
    :param option: dictionary that contains the main option variables of FAS
    :param stack: contains the starting point in the graph and sub-path
    :return: max_fixture, timecheck - tmp_calcstart (time taken for calculation)
    """
    tmp_calcstart = time.time()
    max_fixture = ([], (0.0, 0.0, 0.0, 0.0, 0, 0, False), [], protein)
    v_stack, p_stack = [stack[0]], [stack[1]]
    timecheck = time.time()
    while v_stack and (timecheck - tmp_calcstart <= timelimit or option["timelimit"] == 0):
        query_path, v_stack, p_stack = pb_graphtraversal(query_graph, v_stack, p_stack, option)
        # case M2.2.1: graph(query)-empty(search)
        timecheck = time.time()
        if int(len(search_features)) == 0:
            logging.warning("CASE M2.2.1: graph vs empty.")
            # special case: protein with no pfam or smart domains
            # get score for a_s_f and query_path directly
            path = list(a_s_f.keys())
            query_path_ad = query_path + list(a_q_f.keys())
            score_w = sf_entire_calc_score(path, query_path_ad, weights, search_features, a_s_f,
                                           query_features, a_q_f, clan_dict, option)
        else:
            # case M2.2.2 graph(query)-graph(search)
            # regular traversal of graph based on search_protein
            logging.warning("CASE M2.2.2: graph vs graph.")
            if option["timelimit"] >= 1:
                tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, query_path, search_features,
                                                          weights, query_features, seed_proteome, a_s_f, a_q_f,
                                                          clan_dict, query_clans, protein_lengths,
                                                          timelimit - (timecheck - tmp_calcstart), option)
            else:
                tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, query_path, search_features,
                                                          weights, query_features, seed_proteome, a_s_f, a_q_f,
                                                          clan_dict, query_clans, protein_lengths, "OFF", option)
            path = tmp_path_score[0][0]
            score_w = tmp_path_score[0][1]
            query_path_ad = tmp_path_score[0][2]
            logging.debug("tmp_path_score " + str(tmp_path_score))  # search path, score, mode

        # check for max scoring fixture of path and query_path
        if score_w[5]:
            if score_w[3] >= max_fixture[1][3]:
                max_fixture = (path, score_w, query_path_ad, protein)
                logging.warning("max_fixture: " + str(max_fixture))
        else:
            if (not max_fixture[1][5]) and (score_w[4] >= max_fixture[1][4]):
                max_fixture = (path, score_w, query_path_ad, protein)
        timecheck = time.time()

    return max_fixture, timecheck - tmp_calcstart


def pb_calc_sub(search_protein, protein, search_features, weights, query_features, seed_proteome, a_s_f, a_q_f,
                clan_dict, query_clans, protein_lengths, timelimit, option, jobpaths):
    """Sub function for multiprocessing [2]: Goes through a list of (query) paths an evaluates them.

    :param search_protein: contains feature architecture of the current seed protein [dictionary]
    :param protein: String that contains the identifier of the seed protein
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param timelimit: the amount of time this job is allowed to take before the calculation is stopped and priority
                      mode is used instead
    :param option: dictionary that contains the main option variables of FAS
    :param jobpaths: a group of (query) paths to be evaluated
    :return: max_fixture, tmpcheck - tmpstart (time taken for calculation)
    """
    max_fixture = ([], (0.0, 0.0, 0.0, 0.0, 0, 0, False), [], protein)
    tmpstart = time.time()
    for jobpath in jobpaths:
        tmpcheck = time.time()
        tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, jobpath, search_features, weights,
                                                  query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                                  protein_lengths, timelimit - (tmpcheck - tmpstart), option)
        path = tmp_path_score[0][0]
        score_w = tmp_path_score[0][1]
        query_path_ad = tmp_path_score[0][2]
        logging.debug("tmp_path_score " + str(tmp_path_score))  # search path, score, mode
        if score_w[5]:
            if score_w[3] >= max_fixture[1][3]:
                max_fixture = (path, score_w, query_path_ad, protein)
                logging.warning("max_fixture: " + str(max_fixture))
        else:
            if (not max_fixture[1][5]) and (score_w[4] >= max_fixture[1][4]):
                max_fixture = (path, score_w, query_path_ad, protein)
    tmpcheck = time.time()
    return max_fixture, tmpcheck - tmpstart


def pb_entire_main_nongreedy(search_protein, protein_id, query_path, search_features, weights, query_features,
                             seed_proteome, a_s_f, a_q_f, clan_dict, query_clans, protein_lengths, tmp_timelimit,
                             option):
    """Main Path-building function,
    creates graph with all paths(for seed/search protein), checks priority mode activation if necessary retrieves best
    path from graph traversal function, returns best path, score and mode
    Function calls: pb_region_mapper(), pb_region_paths_nongreedy(), pb_entire_priority_mode(),
                    pb_entire_graphtraversal()

    :param search_protein: contains feature architecture of the current seed protein [dictionary]
    :param protein_id: String that contains the identifier of the seed protein
    :param query_path: currently evaluated query path
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param option: dictionary that contains the main option variables of FAS
    :return: path_score, mode[priority or extensive]
    """
    logging.info("pb_entire_main_nongreedy")
    path_score = 0
    mode = 0
    priority_check = False

    region = pb_region_mapper(list(search_protein), search_features, option["max_overlap"],
                              option["max_overlap_percentage"])
    search_graph = pb_region_paths_nongreedy(region)
    logging.debug(region)
    logging.debug(search_graph)

    if (int(len(search_features)) >= int(option["priority_threshold"]) and option["priority_mode"]) \
            or tmp_timelimit == "OVER":
        mode = 1
        priority_check = True
        # checking best path for all domain types
        path_score = pb_entire_priority_mode(protein_id, query_path, search_graph, search_features, weights,
                                             query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                             protein_lengths, option)
    elif option["priority_mode"]:
        # PRIORITY CHECK "2": for every protein in seed_proteome
        tmp_path_set_size = len(pb_graphtraversal(search_graph, [], [], option))
        logging.debug("cardinality of tmp graph: " + str(tmp_path_set_size))
        if int(tmp_path_set_size) > int(option["max_cardinality"]):
            mode = 1
            priority_check = True
            path_score = pb_entire_priority_mode(protein_id, query_path, search_graph, search_features, weights,
                                                 query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                                 protein_lengths, option)

    if not priority_check:
        mode = 0
        path_score = pb_entire_graphtraversal(search_graph, query_path, search_features, weights, query_features,
                                              a_s_f, a_q_f, clan_dict, tmp_timelimit, option)
    return path_score, mode


def pb_entire_priority_mode(protein, query_path, search_graph, search_features, weights, query_features, seed_proteome,
                            a_s_f, a_q_f, clan_dict, query_clans, protein_lengths, option):
    """Path-evaluation in priority mode

    :param protein: String that contains the identifier of the seed protein
    :param query_path: currently evaluated query path
    :param search_graph: graph containing all paths through the seed architecture
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param option: dictionary that contains the main option variables of FAS
    :return: best_path
    """
    logging.info("pb_entire_priority_mode")

    # best_path (Path, (MS_score(float), PS_score(float), CS_score(float), final_score(float), path_weight(float),
    # common_feature(bool)), QPath)
    best_path = ("NULL", (0.0, 0.0, 0.0, 0.0, 0.0, False), "NULL")
    priority_list = []
    for domain_type in seed_proteome[protein]:
        if seed_proteome[protein][domain_type][0]:
            priority_list.append(domain_type)
    priority_list.append("NONE")
    for domain_type in priority_list:
        path_score_w = pb_entire_graphtraversal_priority(search_graph, domain_type, query_path, 0, search_features,
                                                         weights, query_features, seed_proteome, a_s_f, a_q_f,
                                                         clan_dict, query_clans, protein_lengths, option)
        if path_score_w[1][5]:
            if path_score_w[1][3] >= best_path[1][3]:
                best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
        else:
            # if so far no path with common features found AND the weight is higher
            if (not best_path[1][5]) and (path_score_w[1][4] >= best_path[1][4]):
                best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
    return best_path


def pb_region_paths_nongreedy(overlap_map):
    """Uses overlap information to build a directional, acyclic graph that contains all linear paths

    :param overlap_map: contains overlap information
    :return: graph
    """
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


def pb_region_mapper(overlap_region, features, max_overlap, max_overlap_percentage):
    """ Finds all (relevant) overlaps in the feature architecture

    :param overlap_region: region that is looked at
    :param features: feature information [start/stop]
    :return: overlap_map
    """
    logging.info("pb_region_mapper")
    logging.debug(overlap_region)

    overlap_map = [("START", overlap_region)]
    for i in range(0, len(overlap_region)):
        end = features[overlap_region[i]][3]
        length_i = features[overlap_region[i]][3] - features[overlap_region[i]][2]
        length_i = length_i * max_overlap_percentage
        overlap = []
        x = i + 1
        while x < len(overlap_region):
            length_x = features[overlap_region[x]][3] - features[overlap_region[x]][2]
            length_x = length_x * max_overlap_percentage
            overlap_size = end - features[overlap_region[x]][2]
            if overlap_size <= max_overlap and overlap_size <= length_i and overlap_size <= length_x:
                overlap.append(overlap_region[x])
            x += 1
        overlap.append("END")
        overlap_map.append((overlap_region[i], overlap))
    overlap_map.append(("END", []))
    overlap_map.reverse()
    return overlap_map


def pb_entire_graphtraversal(search_graph, query_path, search_features, weights, query_features, a_s_f, a_q_f,
                             clan_dict, timelimit, option):
    """Traverses the feature architecture graph (seed) to score all paths

    :param search_graph: graph containing all paths through the seed architecture
    :param query_path: currently evaluated query path
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param option: dictionary that contains the main option variables of FAS
    :return: best_path
    """
    logging.info("pb_entire_graphtraversal")
    logging.debug("search graph: " + str(search_graph))

    calcstart = time.time()
    v_stack = ["START"]
    p_stack = []
    best_path = ([], (0.0, 0.0, 0.0, 0.0, 0.0, False, 0.0), [])
    first = 1
    timecheck = time.time()
    if timelimit == "OFF":
        timelimit_off = True
        timelimit = 0.0
    else:
        timelimit_off = False
    if len(query_path) > 0:
        while v_stack and (timecheck - calcstart < timelimit or first or timelimit_off):
            vertex = v_stack.pop()
            if len(p_stack) == 0:
                path = []
            else:
                path = p_stack.pop()
            for next_vertex in search_graph[vertex]:
                if next_vertex == "END":
                    path_ad = path + list(a_s_f.keys())
                    query_path_ad = query_path + list(a_q_f.keys())
                    score_w = sf_entire_calc_score(path_ad, query_path_ad, weights, search_features, a_s_f,
                                                   query_features, a_q_f, clan_dict, option)
                    logging.info("search path " + str(path_ad) + " in pb_entire_graphtraversal")
                    logging.debug("Score info: " + str(score_w) + " for query_path " + str(
                        query_path_ad) + " in pb_entire_graphtraversal")
                    if score_w[5]:
                        if score_w[3] >= best_path[1][3]:
                            best_path = (path_ad, score_w, query_path_ad)
                            logging.debug("new best_path by score: " + str(best_path))
                    else:
                        if (not best_path[1][5]) and score_w[4] >= best_path[1][4]:
                            best_path = (path_ad, score_w, query_path_ad)
                            logging.debug("new best_path by weight: " + str(best_path))
                    first = 0
                else:
                    v_stack.append(next_vertex)
                    p_stack.append(path + [next_vertex])
            timecheck = time.time()
    else:
        while v_stack and first:
            vertex = v_stack.pop()
            if len(p_stack) == 0:
                path = []
            else:
                path = p_stack.pop()
            for next_vertex in search_graph[vertex]:
                if next_vertex == "END":
                    path_ad = path + list(a_s_f.keys())
                    query_path_ad = query_path + list(a_q_f.keys())
                    score_w = sf_entire_calc_score(path_ad, query_path_ad, weights, search_features, a_s_f,
                                                   query_features, a_q_f, clan_dict, option)
                    logging.info("search path " + str(path_ad) + " in pb_entire_graphtraversal")
                    logging.debug("Score info: " + str(score_w) + " for query_path " + str(
                        query_path_ad) + " in pb_entire_graphtraversal")
                    if score_w[5]:
                        if score_w[3] >= best_path[1][3]:
                            best_path = (path_ad, score_w, query_path_ad)
                            logging.debug("new best_path by score: " + str(best_path))
                    else:
                        if (not best_path[1][5]) and score_w[4] >= best_path[1][4]:
                            best_path = (path_ad, score_w, query_path_ad)
                            logging.debug("new best_path by weight: " + str(best_path))
                    first = 0
                else:
                    v_stack.append(next_vertex)
                    p_stack.append(path + [next_vertex])
    logging.info("return: " + str(best_path))
    return best_path


def pb_graphtraversal(graph, v_stack, p_stack, option):
    """Traverses the feature architecture graph (query) to evaluate number of paths, returns all paths if # of paths is
     smaller than the max_cardinality

    :param graph: graph containing all paths through the (query) architecture
    :param option: dictionary that contains the main option variables of FAS
    :return: paths
    """
    logging.debug(graph)
    paths = []
    if v_stack == [1]:
        v_stack = ["START"]
        p_stack = []
        mode = 1
    elif not v_stack:
        v_stack = ["START"]
        p_stack = []
        mode = 0
    else:
        mode = 1
    while v_stack:
        vertex = v_stack.pop()
        if len(p_stack) == 0:
            path = []
        else:
            path = p_stack.pop()
        for next_vertex in graph[vertex]:
            if next_vertex == "END":
                if mode == 0:
                    paths.append(path)
                    logging.info("cardinality: " + str(len(paths)))
                    if len(paths) > option["max_cardinality"]:
                        return paths
                else:
                    return path, v_stack, p_stack
            else:
                v_stack.append(next_vertex)
                p_stack.append(path + [next_vertex])
    logging.debug("returning paths: " + str(paths))
    logging.info("cardinality: " + str(len(paths)))
    return paths


# NEEDS TO BE UNTANGLED (mode 0 + mode 1, query, seed)
def pb_entire_graphtraversal_priority(search_graph, priority, query_path, mode, search_features, weights,
                                      query_features, seed_proteome, a_s_f, a_q_f, clan_dict, query_clans,
                                      protein_lengths, option):
    """Traverses the feature architecture graph in priority mode

    :param search_graph: graph containing all paths through the seed architecture
    :param priority: current priority feature
    :param query_path: graph containing all paths through the query architecture
    :param mode: query or seed linearization
    :param search_features: all features to be linearized in the seed protein
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param option: dictionary that contains the main option variables of FAS
    :return: paths or best_path
    """
    logging.info("pb_entire_graphtraversal_priority")
    logging.debug(
        "search_features: " + str(search_features) + "\nsearch_graph: " + str(search_graph) + "\tpriority: " + str(
            priority) + "\ta_s_f: " + str(a_s_f) + "\tquery_path: " + str(query_path) + "\ta_q_f: " + str(a_q_f))
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
        for next_vertex in search_graph[vertex]:
            if next_vertex == "END":
                if mode == 1:
                    paths.append(path)
                    p_found = 2
                else:
                    path_ad = path + list(a_s_f.keys())
                    query_path_ad = query_path + list(a_q_f.keys())
                    logging.debug("path_ad: " + str(path_ad))
                    logging.debug("query_path_ad: " + str(query_path_ad))

                    score_w = sf_entire_calc_score(path_ad, query_path_ad, weights, search_features, a_s_f,
                                                   query_features, a_q_f, clan_dict, option)
                    if score_w[5]:
                        if score_w[3] >= best_path[1][3]:
                            best_path = (path_ad, score_w, query_path_ad)
                    else:
                        if (not best_path[1][5]) and (score_w[4] >= best_path[1][4]):
                            best_path = (path_ad, score_w, query_path_ad)
                    p_found = 2
            elif mode == 1:
                if query_features[next_vertex][0] == priority:
                    p_candidates.append(next_vertex)
                    p_found = 1
            elif mode == 0:
                p_candidates.append(next_vertex)
                p_found = 1
        if p_found == 1:
            if len(p_candidates) == 1:
                v_stack.append(p_candidates[0])
                p_stack.append(path + p_candidates)
            else:
                if mode == 1:
                    best_priority_bridger = ("NONE", (0.0, 0.0, 0.0, 0.0, 0.0, False))
                    for next_vertex in p_candidates:
                        # baustelle: fixed: path elongation without additives
                        score = sf_calc_score(path + [next_vertex], protein, weights, search_features,
                                              query_features, seed_proteome, clan_dict, query_clans,
                                              protein_lengths, option)
                        if score[3] >= best_priority_bridger[1][3]:
                            best_priority_bridger = (next_vertex, score)
                    v_stack.append(best_priority_bridger[0])
                    p_stack.append(path + [best_priority_bridger[0]])
                elif mode == 0:
                    # greedy strategy: if feature type priority appears more than once
                    best_partial_path = ("NONE", (0.0, 0.0, 0.0, 0.0, 0.0, False))
                    for next_vertex in p_candidates:
                        path_ad = path + list(a_s_f.keys())
                        query_path_ad = query_path + list(a_q_f.keys())
                        score_w = sf_entire_calc_score(path_ad + [next_vertex], query_path_ad, weights,
                                                       search_features, a_s_f, query_features, a_q_f, clan_dict, option)
                        if score_w[5]:
                            if score_w[3] >= best_partial_path[1][3]:
                                best_partial_path = (next_vertex, score_w)
                        else:
                            if (not best_partial_path[1][5]) and (score_w[4] >= best_partial_path[1][4]):
                                best_partial_path = (next_vertex, score_w)

                    v_stack.append(best_partial_path[0])
                    p_stack.append(path + [best_partial_path[0]])

        elif p_found == 0:
            if mode == 1:
                best_priority_bridger = ("NONE", (0.0, 0.0, 0.0, 0.0, 0.0, False))
                for next_vertex in search_graph[vertex]:
                    score = sf_calc_score(path + [next_vertex], protein, weights, search_features, query_features,
                                          seed_proteome, clan_dict, query_clans, protein_lengths, option)
                    if score[3] >= best_priority_bridger[1][3]:
                        best_priority_bridger = (next_vertex, score)
                v_stack.append(best_priority_bridger[0])
                p_stack.append(path + [best_priority_bridger[0]])

            elif mode == 0:
                # some kind of greedy strategy: if feature type (p priority) not found
                best_partial_path = ("NONE", (0.0, 0.0, 0.0, 0.0, 0.0, False))
                for next_vertex in search_graph[vertex]:
                    path_ad = path + list(a_s_f.keys())
                    query_path_ad = query_path + list(a_s_f.keys())
                    score_w = sf_entire_calc_score(path_ad + [next_vertex], query_path_ad, weights, search_features,
                                                   a_s_f, query_features, a_q_f, clan_dict, option)
                    if score_w[5]:
                        if score_w[3] >= best_partial_path[1][3]:
                            best_partial_path = (next_vertex, score_w)
                    else:
                        if (not best_partial_path[1][5]) and (score_w[4] >= best_partial_path[1][4]):
                            best_partial_path = (next_vertex, score_w)

                v_stack.append(best_partial_path[0])
                p_stack.append(path + [best_partial_path[0]])
    if mode == 1:
        logging.debug("returning: " + str(paths))
        logging.debug("cardinality: " + str(len(paths)))
        return paths
    elif mode == 0:
        logging.debug(best_path)
        return best_path


# Scoring Functions # <sf>

def sf_calc_score(path, protein, weights, search_features, query_features, seed_proteome,
                  clan_dict, query_clans, protein_lengths, option):
    """Used to be main Scoring function, starts the functions for the scores and calculates the complete FAS score,
    only used during priority mode to make greedy decisions

    :param path: Path to score
    :param protein: protein id
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param search_features: all features to be linearized in the seed protein
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: (Pfam)-clans of the current query protein
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param option: dictionary that contains the main option variables of FAS
    :return: score_ms, score_ps, score_cs, final_score
    """
    adjusted_weights = {}
    tmp_weight = {}
    if option["weight_const"] == 1:
        adjusted_weights = w_weight_const_rescale(path, weights, search_features, False, option)
        tmp_weight = {}
        for i in adjusted_weights:
            tmp_weight[i] = weights[i]
            weights[i] = adjusted_weights[i]
    score_cs = round(sf_cs_score(path, clan_dict, query_clans, query_features), 10)
    tmp = sf_ms_score(path, protein, seed_proteome, query_features, option)
    score_ms = round(tmp[0], 10)
    score_ps = round(
        sf_ps_score(path, tmp[2], protein, query_features, seed_proteome, protein_lengths), 10)
    final_score = (score_ms * option["score_weights"][0]) + (score_cs * option["score_weights"][1]) + (
        score_ps * option["score_weights"][2])
    if option["weight_const"] == 1:
        for i in adjusted_weights:
            weights[i] = tmp_weight[i]
    return score_ms, score_ps, score_cs, final_score


def sf_entire_calc_score(path, query_path, weights, search_features, a_s_f, query_features, a_q_f, clan_dict, option):
    """Main Scoring function, starts the functions for the scores and calculates the complete FAS score

    :param path: Path to score
    :param query_path: Path to score (query)
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param search_features: all features to be linearized in the seed protein
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param clan_dict: dictionary that maps features to clans
    :param option: dictionary that contains the main option variables of FAS
    :return: score_ms, score_ps, score_cs, final_score, path_weight, common_feature, score_ls
    """
    adjusted_weights = {}
    tmp_weight = {}
    if option["weight_const"] == 1:
        adjusted_weights = w_weight_const_rescale(path, weights, search_features, False, option)
        tmp_weight = {}
        for i in adjusted_weights:
            tmp_weight[i] = weights[i]
            weights[i] = adjusted_weights[i]
    score_cs = round(sf_entire_cs_score(path, query_path, query_features, clan_dict, search_features), 10)
    tmp = sf_entire_ms_score(path, query_path, search_features, a_s_f, query_features, a_q_f, weights, option)
    score_ms = round(tmp[0], 10)
    path_weight = tmp[3]
    common_feature = tmp[4]
    ps_tmp = sf_entire_ps_score(path, tmp[2], query_path, search_features, a_s_f, query_features, a_q_f, weights,
                                option)
    score_ps = round(ps_tmp[0], 10)
    score_ls = round(ps_tmp[1], 10)
    final_score = (score_ms * option["score_weights"][0]) + (score_cs * option["score_weights"][1]) + (
        score_ps * option["score_weights"][2])
    if option["weight_const"] == 1:
        for i in adjusted_weights:
            weights[i] = tmp_weight[i]
    return score_ms, score_ps, score_cs, final_score, path_weight, common_feature, score_ls


def sf_cs_score(path, clan_dict, query_clans, features):
    """Calculates clan score

    :param path: Path to score
    :param clan_dict: dictionary that maps features to clans
    :param query_clans: Clans in the query
    :param features: features (seed or query) [dict]
    :return: score (PS)
    """

    counter = 0.0
    path_clans = {}
    score = 0.0

    # counting clans in path
    # path_clans: contains counts for clans in path
    for i in path:
        feature = features[i]

        if clan_dict[feature[0]] != "NO_CLAN":
            if clan_dict[feature[0]] in path_clans:
                path_clans[clan_dict[feature[0]]] += 1
            else:
                path_clans[clan_dict[feature[0]]] = 1
    for clan in path_clans:
        counter += 1.0
        if clan in query_clans:
            score += float(path_clans[clan] * query_clans[clan]) / float(
                max(path_clans[clan], query_clans[clan]) * max(path_clans[clan], query_clans[clan]))
            logging.debug("(" + str(path_clans[clan]) + " * " + str(query_clans[clan]) + ") / (max(" + str(
                path_clans[clan]) + "," + str(query_clans[clan]) + ") * max(" + str(path_clans[clan]) + "," + str(
                query_clans[clan]) + ")")
    if counter == 0:
        score = 0
    else:
        score = score / counter
    return score


def sf_entire_cs_score(path, query_path, query_features, clan_dict, search_features):
    """calculates clan score

    :param path: Path to score
    :param query_features: all features to be linearized in the query protein
    :param search_features: all features to be linearized in the seed protein
    :param clan_dict: dictionary that maps features to clans
    :return: score (PS)
    """
    logging.debug("search_path: " + str(path))
    logging.debug("query_path: " + str(query_path))

    counter = 0.0
    s_clans = {}  # clans from path
    q_clans = {}  # clans from query_path

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
            score += float(s_clans[clan] * q_clans[clan]) / float(
                max(s_clans[clan], q_clans[clan]) * max(s_clans[clan], q_clans[clan]))
    if counter == 0:
        score = 0
    else:
        score = score / counter
    return score


def sf_ms_score(path, protein, proteome, features, option):
    """Calculates multiplicity score, only used for priority mode now

    :param path : Path to score
    :param protein : protein id
    :param proteome : seed_proteome is a dictionary that contains the feature architecture of all seed proteins
    :param features : features (seed or query) [dict]
    :param option : specifies the behavior of the MS calculation. Can be switched between classicMS and the new version.
    :return: final_score (MS), search_domains, scale
    """
    domains = {}
    scale = 0
    scores = []
    final_score = 0
    for i in path:
        feature = features[i]
        if feature[0] in domains:
            domains[feature[0]] += 1
        else:
            domains[feature[0]] = 1
    for feature in domains:
        scale += 1
        if feature in proteome[protein]:
            s_length = len(proteome[protein][feature]) - 2
            if option["classicMS"]:
                p_score = float(domains[feature] * s_length) / \
                          float(max(domains[feature], s_length) * max(domains[feature], s_length))
            else:
                p_score = min(float(domains[feature] * s_length) / float(domains[feature] * domains[feature]), 1.0)
            scores.append((feature, p_score))
        else:
            scores.append((feature, 0.0))
    if scale > 0:
        scale = 1.0 / float(scale)
    for score in scores:
        final_score += score[1] * scale
    return final_score, domains, scale


def sf_entire_ms_score(path, query_path, search_features, a_s_f, query_features, a_q_f, weights, option):
    """Calculates multiplicity score

    :param path: Path to score
    :param query_path: Path to score (query)
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param search_features: all features to be linearized in the seed protein
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param option: dictionary that contains the main option variables of FAS
    :return: final_score (MS), search_domains, scale, final_weight, common_feature
    """
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

    # check for common features
    for y in main_features_q:
        try:
            if main_features_s[y]:
                common_feature = True
        except KeyError:
            logging.debug("Feature not found: " + str(y))

    for feature in search_domains:
        if option["MS_uni"] == 0:
            scale += weights[feature]
        else:
            scale += 1
        if feature in query_domains:
            s_length = query_domains[feature]
            if option["classicMS"]:
                p_score = float(search_domains[feature] * s_length) / \
                          float(max(search_domains[feature], s_length) * max(search_domains[feature], s_length))
            else:
                p_score = min(float(search_domains[feature] * s_length) /
                              float(search_domains[feature] * search_domains[feature]), 1.0)
            scores.append((feature, p_score))
        else:
            scores.append((feature, 0.0))
    if scale > 0:
        scale = 1.0 / float(scale)
    for score in scores:
        if option["MS_uni"] == 0:
            final_score += score[1] * scale * weights[score[0]]
            final_weight += weights[score[0]]
        else:
            final_score += score[1] * scale
    logging.debug(
        "Return entire_ms_score: " + str(final_score) + ", " + str(search_domains) + ", " + str(scale) + ", " + str(
            final_weight) + ", " + str(common_feature))
    return final_score, search_domains, scale, final_weight, common_feature


def sf_ps_score(path, scale, protein, features, seed_proteome, protein_lengths):
    """Calculates positional score

    :param path: Path to score
    :param scale: contains scaling factor for each weight or number of feature types for uniform weighting
    :param protein: protein id
    :param features: features (seed or query) [dict]
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :return: final_score (PS)
    """
    count = {}
    final_score = 0.0
    scores = {}
    for i in path:
        feature = features[i]
        if feature[0] in seed_proteome[protein]:
            best_match = 0.0
            if not feature[0] in scores:
                scores[feature[0]] = 0.0
                count[feature[0]] = 0
            for instance in seed_proteome[protein][feature[0]][2:]:
                pos = (float(instance[1]) + float(instance[2])) / 2.0 / float(
                    protein_lengths["seed_" + str(protein)])
                match = 1.0 - float(abs(feature[1]) - pos)
                logging.debug(str(float(instance[1])) + " + " + str(float(instance[2])) + " / 2.0 / " + str(
                    float(protein_lengths["seed_" + str(protein)])) + " = " + str(pos))
                if best_match < match:
                    best_match = match
            scores[feature[0]] += best_match
            count[feature[0]] += 1

    for f_score in scores:
        final_score += scores[f_score] / count[f_score] * scale
    return final_score


def sf_entire_ps_score(path, scale, query_path, search_features, a_s_f, query_features, a_q_f, weights, option):
    """Calculates positional score

    :param path: Path to score
    :param scale: contains scaling factor for each weight or number of feature types for uniform weighting
    :param query_path: Path to score (query)
    :param weights: feature weights
    :param query_features: all features to be linearized in the query protein
    :param search_features: all features to be linearized in the seed protein
    :param a_s_f: additional seed features [non linearized]
    :param a_q_f: additional query features [non linearized]
    :param option: dictionary that contains the main option variables of FAS
    :return: final_score (PS), final_ls_score
    """
    count = {}
    final_score = 0.0
    final_ls_score = 0.0
    scores = {}
    ls_scores = {}
    local_query_protein = {}
    ls = 0
    # get current features from query path
    for i in query_path:
        if i in query_features:
            feature = query_features[i]
        else:
            feature = a_q_f[i]
        if feature[0] in local_query_protein:
            length = feature[3] - feature[2]
            local_query_protein[feature[0]].append((feature[1], length))
        else:
            local_query_protein[feature[0]] = [(feature[1], feature[3] - feature[2])]
    logging.debug(str(local_query_protein) + " for query_path " + str(query_path))

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
        if option["MS_uni"] == 0:
            final_score += scores[f_score] / count[f_score] * scale * weights[f_score]
            final_ls_score += ls_scores[f_score] / count[f_score] * scale * weights[f_score]
        else:
            final_score += scores[f_score] / count[f_score] * scale
            final_ls_score += ls_scores[f_score] / count[f_score] * scale

    return final_score, final_ls_score


########## Weighting ########## <w>
# weighting functions


def w_count(prot_count, domain_count, seed_proteome, query_proteome):
    """Goes through all feature in query_proteome and seed_proteome and groups them by rarity (according to the
    reference), also adds features that are not in the reference proteome with a count of 1

    :param prot_count: number of proteins in reference
    :param domain_count: count for each feature in reference
    :param seed_proteome: dictionary that contains the feature architecture of all seed proteins
    :param query_proteome: dictionary that contains the feature architecture of all query proteins
    :return: relevant_feature, domain_count
    """
    logging.debug("w_count: processing query and search domains.")
    relevant_features = {}
    for query in query_proteome:
        for feature in query_proteome[query]:
            if feature not in domain_count:
                domain_count[feature] = 1
                relevant_features[feature] = 1
            elif feature not in relevant_features:
                if domain_count[feature] == 1:
                    relevant_features[feature] = 1
                elif float(domain_count[feature]) / prot_count <= 0.1:
                    relevant_features[feature] = domain_count[feature]
                elif float(domain_count[feature]) / prot_count <= 0.5:
                    relevant_features[feature] = domain_count[feature]
                else:
                    relevant_features[feature] = domain_count[feature]
    for template in seed_proteome:
        protein = seed_proteome[template]
        for feature in protein:
            if feature not in domain_count:
                domain_count[feature] = 1
                relevant_features[feature] = 1
            elif feature not in relevant_features:
                if float(domain_count[feature]) / prot_count <= 0.1:
                    relevant_features[feature] = domain_count[feature]
                if float(domain_count[feature]) / prot_count <= 0.5:
                    relevant_features[feature] = domain_count[feature]
                else:
                    relevant_features[feature] = domain_count[feature]

    logging.debug("domain counts: " + str(domain_count))
    return relevant_features, domain_count


def w_count_ref(proteome):
    """Counts all features in reference proteome

    :param proteome: dictionary that contains the feature architecture of all ref proteins
    :return: prot_count, domain_count
    """
    logging.debug("w_count_ref: counting domains in reference gene set.")

    domain_count = {}
    prot_count = 0
    for i in proteome:
        prot_count += 1
        protein = proteome[i]
        for feature in protein:
            # counting instances (subtracts 2 because first and second entry contain assess(bool) and feat_eval(float))
            count = len(protein[feature]) - 2
            if feature in domain_count:
                domain_count[feature] += count
            else:
                domain_count[feature] = count
    logging.debug("domains counts: " + str(domain_count))
    return float(prot_count), domain_count


def w_weighting(protein, domain_count, proteome):
    """Calculates weights

    :param protein: protein id
    :param domain_count: feature counts of the reference
    :param proteome: seed or query proteome
    :return: weights, domain_count
    """
    weights = {}
    features = []
    scaling_factor = 0.0
    sum_of_features = 0.0
    for feature in proteome[protein]:
        features.append(feature)
    for feature in features:
        try:
            domain_count[feature]
        except KeyError:
            # print "feature not found in ref " + feature
            domain_count[feature] = 1
            sum_of_features += float(domain_count[feature])
        else:
            sum_of_features += float(domain_count[feature])
    for feature in features:
        weights[feature] = round(float(sum_of_features) / float(domain_count[feature]), 8)
        scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 8)
    for feature in features:
        weights[feature] = round(float(weights[feature]) / float(scaling_factor), 8)
    return weights, domain_count


def w_weighting_constraints(protein, domain_count, proteome, option):
    """Calculates weights, while fulfilling the constraints

    :param protein: protein id
    :param domain_count: feature counts of the reference
    :param proteome: seed or query proteome
    :param option: dictionary that contains the main option variables of FAS
    :return: weights, domain_count
    """
    weights = {}
    tools = {}
    for tool in (option["input_linearized"] + option["input_normal"]):
        tools[tool] = []
    features = []
    single_constraints = []
    filled = 0.0
    for feature in proteome[protein]:
        features.append(feature)
    for feature in features:
        if feature in option["constraints"]:
            filled += option["constraints"][feature]
            weights[feature] = option["constraints"][feature]
            single_constraints.append(feature)
        elif feature.split('_')[0] in option["constraints"]:
            tools[feature.split('_')[0]].append(feature)
    for tool in tools:
        if len(tools[tool]) > 0:
            filled += option["constraints"][tool]
            sum_of_features = 0
            scaling_factor = 0.0
            for feature in tools[tool]:
                features.remove(feature)
                try:
                    domain_count[feature]
                except KeyError:
                    domain_count[feature] = 1
                    sum_of_features += float(domain_count[feature])
                else:
                    sum_of_features += float(domain_count[feature])
            for feature in tools[tool]:
                weights[feature] = round(float(sum_of_features) / float(domain_count[feature]), 8)
                scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 8)
            for feature in tools[tool]:
                weights[feature] = round(float(weights[feature]) / float(scaling_factor) * option["constraints"][tool],
                                         8)
    for feature in single_constraints:
        features.remove(feature)
    for feature in features:
        sum_of_features = 0.0
        try:
            domain_count[feature]
        except KeyError:
            # print "feature not found in ref " + feature
            domain_count[feature] = 1
            sum_of_features += float(domain_count[feature])
        else:
            sum_of_features += float(domain_count[feature])
    scaling_factor = 0.0
    for feature in features:
        weights[feature] = round(float(sum_of_features) / float(domain_count[feature]), 8)
        scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 8)
    for feature in features:
        weights[feature] = round(float(weights[feature]) / float(scaling_factor) * (1.0 - filled), 8)
    return weights, domain_count


def w_weight_correction(method, domain_count):
    """Rescales counts by applying one of the 4 functions (loge, log10, root4, root8)

    :param method: rescaling-function
    :param domain_count: feature counts of the reference
    :return: domain_count
    """
    if method == "loge":
        for feature in domain_count:
            domain_count[feature] = int(round(math.log(domain_count[feature]), 0) + 1)
    elif method == "log10":
        for feature in domain_count:
            domain_count[feature] = int(round(math.log10(domain_count[feature]), 0) + 1)
    elif method == "root4":
        for feature in domain_count:
            domain_count[feature] = int(round(math.pow(domain_count[feature], 0.25), 0))
    elif method == "root8":
        for feature in domain_count:
            domain_count[feature] = int(round(math.pow(domain_count[feature], 0.125), 0))
    return domain_count


def w_weight_const_rescale(path, weights, search_features, final, option):
    """Rescales weights to 1 after linearization

    :param path: Path to rescale
    :param weights: feature weights
    :param search_features: all features to be linearized in the seed protein
    :param final: final calculation? (1/0)
    :param option: dictionary that contains the main option variables of FAS
    :return: rescaled_weights
    """
    lindict = {}
    for ftype in option["input_linearized"]:
        lindict[ftype] = []
    if final:
        for ftype in option["input_normal"]:
            lindict[ftype] = []
    tmp = 0.0
    rescaled_weights = {}
    if final:
        for feature in path:
            if feature not in option["constraints"]:
                if feature not in lindict[feature.split('_')[0]]:
                    lindict[feature.split('_')[0]].append(feature)
        for ftype in option["input_linearized"]:
            if ftype in option["constraints"]:
                for feature in lindict[ftype]:
                    tmp += weights[feature]
                if option["constraints"][ftype] > tmp > 0.0:
                    scale = option["constraints"][ftype] / tmp
                    for feature in lindict[ftype]:
                        rescaled_weights[feature] = weights[feature] * scale
                tmp = 0.0
    else:
        for feature in path:
            if feature in search_features:
                if search_features[feature][0] not in option["constraints"]:
                    if feature not in lindict[search_features[feature][0].split('_')[0]]:
                        lindict[search_features[feature][0].split('_')[0]].append(feature)

        for ftype in option["input_linearized"]:
            if ftype in option["constraints"]:
                for feature in lindict[ftype]:
                    tmp += weights[search_features[feature][0]]
                if option["constraints"][ftype] > tmp > 0.0:
                    scale = option["constraints"][ftype] / tmp
                    for feature in lindict[ftype]:
                        rescaled_weights[search_features[feature][0]] = weights[search_features[feature][0]] * scale
                tmp = 0.0

    return rescaled_weights


# Input # <OK>

def featuretypes(path, option):
    """
    Input function,
    reads the tools/featuretypes input-file and stores information in option

    :param path: path to input-file
    :param option: dictionary that contains the main option variables of FAS
    :return: option
    """
    option["input_linearized"] = []
    option["input_normal"] = []
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
            option["input_linearized"].append(tmp)
        elif mode == "nor":
            option["input_normal"].append(tmp)
    ifile.close()
    return option


def constraints_in(path):
    """
    Input function,
    reads the constraints file

    :param path: path to input-file
    :return: constraints
    """
    constraints = {}
    cfile = open(path, "r+")
    lines = cfile.readlines()
    i = 1
    if lines[0][0] == "#":
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
    return constraints


def xmlreader(path, mode, tool, assess, proteome, protein_lengths, clan_dict, option):
    """
    Input function,
    read input-files for seed, query and reference

    :param path: path to input file
    :param mode: 0 (seed), 1 (query), 2 (reference)
    :param tool: tool/database name (pfam, seg, ...)
    :param assess: e-value assess
    :param proteome: proteome var to write in [dict]
    :param protein_lengths: dictionary that contains the length of each protein in the seed and query
    :param clan_dict: dictionary that maps features to clans
    :param option: dictionary that contains the main option variables of FAS
    :return: proteome, protein_lengths, clan_dict
    """
    start = 0
    if os.path.exists(path):
        xmltree = ElTre.parse(path)
        root = xmltree.getroot()
        for protein in root:
            p_id = protein.attrib["id"]
            plength = protein.attrib["length"]
            # set up of protein IDs, differentiate between proteins from different files
            if mode == 1:
                protein_lengths["query_" + str(p_id)] = float(plength)
            elif mode == 0:
                protein_lengths["seed_" + str(p_id)] = float(plength)
            elif mode == 2:
                protein_lengths[p_id] = float(plength)
            # set up of datastructure to store annotations
            if not (p_id in proteome):
                proteome[p_id] = {}

            for feature in protein:
                if len(feature) > 0:
                    ftype = tool + "_" + feature.attrib["type"]
                    feat_eval = 'NULL'

                    # evalue check: family/ftype based
                    if 'evalue' in feature.attrib and float(feature.attrib["evalue"]) > option["efilter"]:
                        logging.debug("reject feature type: " + ftype)
                        # print "reject feature"
                        # skip current ftype and continue with the next one
                        continue
                    else:
                        # keep feature type bases evalue
                        if 'evalue' in feature.attrib:
                            feat_eval = float(feature.attrib["evalue"])

                        if 'clan' in feature.attrib:
                            fclan = feature.attrib["clan"]
                        else:
                            fclan = "NO_CLAN"
                        proteome[p_id][ftype] = []
                        proteome[p_id][ftype].append(assess)
                        proteome[p_id][ftype].append(feat_eval)

                        i = 0
                        # counting appended instances
                        inst_count = 0
                        for instance in feature:
                            inst_eval = 'NULL'
                            # XMLcase 1: feature instance contains evalue information (XML field inst_eval)
                            if 'inst_eval' in instance.attrib:
                                # print tool + " instance evalue: "+ str(instance.attrib)
                                inst_eval = float(instance.attrib["inst_eval"])
                                start = int(instance.attrib["start"])
                                end = int(instance.attrib["end"])

                                if inst_eval <= option["inst_efilter"]:
                                    logging.debug("append instance: " + ftype + ": " + str(instance.attrib))

                                    proteome[p_id][ftype].append((inst_eval, start, end))
                                    inst_count += 1

                                else:
                                    logging.debug("reject instance: " + ftype + ": " + str(instance.attrib))

                            # XMLcase 2: feature instance contains NO evalue information (XML field inst_eval)
                            else:
                                # NO instance based evalue information --> no instances can be rejected:
                                # set inst_count = 1
                                inst_count = 1
                                if len(instance.attrib) == 2:
                                    start = int(instance.attrib["start"])
                                    end = int(instance.attrib["end"])
                                    proteome[p_id][ftype].append((inst_eval, start, end))

                                else:
                                    if i == 0:
                                        start = int(instance.attrib["start"])
                                        i = 1
                                    else:
                                        end = int(instance.attrib["end"])
                                        proteome[p_id][ftype].append((inst_eval, start, end))
                                        i = 0
                        # any instance appended?
                        if inst_count < 1:
                            # delete feature type
                            logging.info("Rejecting feature type " + str(ftype) + " due to rejection of all instances. "
                                                                                  "Check for thresholds and E-values "
                                                                                  "(instance based)")
                            proteome[p_id].pop(ftype)
                        if ftype not in clan_dict:
                            clan_dict[ftype] = fclan
    else:
        print("Error: " + path + " does not exist")
        quit()
    logging.debug("proteome: " + str(proteome))
    return proteome, protein_lengths, clan_dict


def bidirectionout(outpath):
    """Output function for bidirectional mode: This function summarizes the scores of the both scoring directions into a
    table (csv format).

    :param outpath: path to the out directory
    :return: no returns
    """
    outdict = {}
    forwardtree = ElTre.parse(outpath + ".xml")
    reversetree = ElTre.parse(outpath + "_reverse.xml")
    forwardroot = forwardtree.getroot()
    reverseroot = reversetree.getroot()
    for query in forwardroot:
        query_id = query.attrib["id"]
        for seed in query:
            seed_id, forward_score, seed_mode = seed.attrib["id"], seed.attrib["score"], seed.attrib["mode"]
            for path in seed:
                if path.tag == "template_path":
                    forward_s_path = {}
                    for feature in path:
                        forward_s_path[feature.attrib["type"]] = []
                        for instance in feature:
                            forward_s_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                           instance.attrib["end"]))
                if path.tag == "query_path":
                    forward_q_path = {}
                    for feature in path:
                        forward_q_path[feature.attrib["type"]] = []
                        for instance in feature:
                            forward_q_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                           instance.attrib["end"]))
            for node in reverseroot:
                if node.attrib["id"] == seed_id:
                    for child in node:
                        if child.attrib["id"] == query_id:
                            reverse_score, query_mode = child.attrib["score"], child.attrib["mode"]
                            for path in child:
                                if path.tag == "query_path":
                                    reverse_s_path = {}
                                    for feature in path:
                                        reverse_s_path[feature.attrib["type"]] = []
                                        for instance in feature:
                                            reverse_s_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                                           instance.attrib["end"]))
                                if path.tag == "template_path":
                                    reverse_q_path = {}
                                    for feature in path:
                                        reverse_q_path[feature.attrib["type"]] = []
                                        for instance in feature:
                                            reverse_q_path[feature.attrib["type"]].append((instance.attrib["start"],
                                                                                           instance.attrib["end"]))
                            consistence = "yes"
                            for feature in forward_q_path:
                                if feature in reverse_q_path:
                                    for instance in forward_q_path[feature]:
                                        if instance not in reverse_q_path[feature]:
                                            consistence = "no"
                                else:
                                    consistence = "no"
                            for feature in forward_s_path:
                                if feature in reverse_s_path:
                                    for instance in forward_s_path[feature]:
                                        if instance not in reverse_s_path[feature]:
                                            consistence = "no"
                                else:
                                    consistence = "no"
                            outdict[(seed_id, query_id)] = (forward_score, reverse_score, consistence,
                                                            seed_mode + "/" + query_mode)
    out = open(outpath + "_table.csv", "w")
    out.write("seedID,queryID,forward,reverse,Path_consistency,linearization_mode\n")
    for pair in outdict:
        out.write(pair[0] + "," + pair[1] + "," + outdict[pair][0] + "," + outdict[pair][1] + "," +
                  outdict[pair][2] + "," + outdict[pair][3] + "\n")
    out.close()


# start #
if options.jobname.find("/") != -1:
    option_dict["outpath"] = options.jobname

else:
    option_dict["outpath"] = expath + "/out/" + su_set_path(options.jobname, expath)

if option_dict["weight_const"] == 1:
    option_dict["constraints"] = constraints_in(options.weight_constraints)
    option_dict["constname"] = options.weight_constraints.rstrip("/").split("/")[-1]

option_dict = featuretypes(options.featuretypes, option_dict)

fc_start(option_dict)

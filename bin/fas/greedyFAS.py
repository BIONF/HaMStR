
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
version = 1
tmp = inspect.getfile(inspect.currentframe())
expath = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
tmp = ""

parser = OptionParser(description="You are running greedyFAS.py version " + str(version) + ".")
parser.add_option("-s", "--singlepath", dest="singlepath", default=expath+"/annotations_out/single", help="Path to the folder containing the xml-files of the single protein, default is fas/in/single")
parser.add_option("-p", "--proteomepath", dest="proteomepath", default=expath+"/annotations_out/proteome", help="Path to the folder containing the xml-files of the protein set, default is fas/in/proteome")
parser.add_option("-o", "--one_vs_all", dest="one_vs_all", default=0, help="If this value is 1, the single protein architecture will be the base of the search, default is 0")
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
(options, args) = parser.parse_args();

### global vars ###             ### global var looks ###
proteome = {}			#{("protein_id", {("domain_name", [("START", "STOP")])})}
single_protein = ()		#("protein_id", {("domain_name", [("START", "STOP")])})
clan_dict = {}			#{("domain_name", "clan")}
search_features = {}		#{("F_0", ("domain_name", "POSITION", "Start", "Stop" ))}
a_s_f = {}                      # additional search features 
query_features = {}             #{("F_0", ("domain_name", "POSITION", "Start", "Stop" ))}
a_q_f = {}                      # additional query features
query_protein = {}		#{("domain_name", ["POSITION_1", "POSITION_2"])}
query_clans = {}		#{("clan", "INSTANCES")}		
weights = {}			#{("domain_name", "weight")}
weights_counter = {}
protein_lengths = {}		#{("protein_id", "LENGTH")}
domain_count = {}		#{("domain", "COUNT")} 
best_paths = {}			#{("protein", ["path"])}
best_entire_paths = 0
best_entire_fixtures = {}
search_graph = 0
query_graph = 0
MS_uni = 0
mode = {}			#{("protein", mode)}
constraints = {}                #{("tool/feature", weight)}
weight_const = 0
weight_correction = 1
## hidden options
# tab separated table as output file #
taciturn = 1
# two-sided linearization #
entire = 1

### READ OPTIONS ###
if int(options.one_vs_all) == 1:
	one_against_all = 1
else:	
	one_against_all = 0
score_weights = options.weights
p_path = options.proteomepath
s_path = options.singlepath
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
        
if options.weight_constraints != 0:
        weight_const = 1

### SETUP LOGGING OPTIONS ###

# possible logging configurations
# logging into stdout with time stamp:
#logging.basicConfig(level=loglevel, format='%(asctime)s - %(levelname)s - %(message)s')
# logging into file with line number:
logging.basicConfig(filename='testlog.log', filemode='w', level=loglevel, format='%(lineno)s - %(levelname)s - %(message)s')
# logging into stdout with line number:
#logging.basicConfig(level=loglevel, format='%(lineno)s - %(levelname)s - %(message)s')

logging.info('greedyFAS.py started with options: entire='+str(entire)+', one_vs_all='+str(one_against_all)+', priority_threshold='+str(priority_threshold)+', log_level='+str(loglevel))
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
	global proteome
	global MS_uni
        global best_entire_paths
        global weight_correction
        
        best_entire_paths = tree()
	globalscores = {}
                
        ## MS_uni set to 0 when no weighting is conducted        
	if MS_uni == 0:
		if os.path.exists(ref_proteome + "/count.ref"):
			ref_reader(ref_proteome + "/count.ref")
		else:
			xmlreader(ref_proteome + "/pfam.xml", 2, "pfam",True)
			xmlreader(ref_proteome + "/smart.xml", 2, "smart",True)
			xmlreader(ref_proteome + "/cast.xml", 2, "cast",True)
			xmlreader(ref_proteome + "/coils.xml", 2, "coils",True)
			xmlreader(ref_proteome + "/seg.xml", 2, "seg",True)
			xmlreader(ref_proteome + "/signalp.xml", 2, "signalp",True)
			xmlreader(ref_proteome + "/tmhmm.xml", 2, "tmhmm",True)
                        ## keep feature counts for reference proteome (gene set)
			w_count_ref()
			######################################
			# reusing feature counts for reference protein sets
			#
			# ref_writer writes into directory with reference protein set
			# may cause "permission denied" error
			######################################
			#ref_writer(ref_proteome + "/count.ref")
                        
                        ## clean gloabl protome
			proteome = {}
		if weight_correction != 0:
                        w_weight_correction(weight_correction)    # use correction function on counts
	xmlreader(p_path + "/pfam.xml", 0, "pfam", True)
	xmlreader(s_path + "/pfam.xml", 1, "pfam", True)
	xmlreader(p_path + "/smart.xml", 0, "smart", True)
	xmlreader(s_path + "/smart.xml", 1, "smart", True)
	xmlreader(p_path + "/cast.xml", 0, "cast",False)
	xmlreader(s_path + "/cast.xml", 1, "cast",False)
	xmlreader(p_path + "/coils.xml", 0, "coils",False)
	xmlreader(s_path + "/coils.xml", 1, "coils",False)
	xmlreader(p_path + "/seg.xml", 0, "seg",False)
	xmlreader(s_path + "/seg.xml", 1, "seg",False)
	xmlreader(p_path + "/signalp.xml", 0, "signalp",False)
	xmlreader(s_path + "/signalp.xml", 1, "signalp",False)
	xmlreader(p_path + "/tmhmm.xml", 0, "tmhmm",False)
	xmlreader(s_path + "/tmhmm.xml", 1, "tmhmm",False)
        if MS_uni == 0:	
		w_count()
	fc_path()
        fc_normal()


# finds best paths
def fc_path():
	global one_against_all
	global proteome
	global best_paths
        global best_entire_paths
        global best_entire_fixtures
	global search_features
        global a_s_f
        global query_features
        global a_q_f
	global single_protein
	global search_graph
        global query_graph
	global MS_uni
	global mode
	global priority_threshold
        global max_cardinality
        global weight_const
        
        logging.info("fc_path")
        
        #score model M2 - used for iteration/incremental evaluation
        ##############################################################################################################################
	
        if one_against_all == 0:
            ######## entire mode ############################################
            # both proteins (search and query) will be linearized ##########
            if MS_uni == 0:
                if weight_const == True:
                    w_weighting_constraints_counter("none", entire) 
                else:
                    w_weighting_counterpart("none", entire)
            # single features kept as query features (dep. on parameter -o)
            # protein (from proteom) kept as search features
            if entire == 1:
                go_priority = False
                lin_single_set = su_lin_query_protein(single_protein[0], 1)
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
                
                for protein in proteome:

                    mode[protein] = 0;
                    pathcount = 0
                    if MS_uni == 0:
                        if weight_const == True:
                            w_weighting_constraints(protein, entire) 
                        else:                        
                            w_weighting(protein, entire)
                    if go_priority == True:
                        all_single_paths = []
                        for domain_type in single_protein[1]:
                            all_single_paths += (pb_entire_graphtraversal_priority(tmp_single_graph, domain_type, protein, 1))
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
                        
                        search_features = {}    # set in su_search_protein
                        a_s_f = {}              # set in su_search_protein
                        search_protein = su_search_protein(protein, 1)
                        logging.warning("search_features(p): "+str(len(search_features))+" for: "+str(protein))
                        
                        # case M2.1.1: empty(query)-empty(search)
                        # should be the best fix independent from weight
                        if int(len(search_features)) == 0:
                            logging.warning("CASE M2.1.1: empty vs empty.")
                            score_w = sf_entire_calc_score(a_s_f.keys(),a_q_f.keys())
                            path = a_s_f.keys()
                            query_architecture = a_q_f.keys()
                            mode[protein] = 2
                        else:
                            # case M2.1.2: empty(query)-graph(search)
                            logging.warning("CASE M2.1.2: empty vs graph.")
                            tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, [])
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
                        search_features = {}
                        a_s_f = {}
                        search_protein = su_search_protein(protein, 1)
                        logging.warning("search_features(p): "+str(len(search_features))+" for: "+str(protein))
                        logging.debug("query_path No."+str(pathcount)+": "+str(query_path))
                        
                        #case M2.2.1: graph(query)-empty(search)
                        if int(len(search_features)) == 0:
                            logging.warning("CASE M2.2.1: graph vs empty.")
                            #special case: protein with no pfam or smart domains
                            # get score for a_s_f and query_path directly
                            path = a_s_f.keys()
                            query_path_ad = query_path+a_q_f.keys()
                            score_w = sf_entire_calc_score(path,query_path_ad)
                            
                        else:
                            #case M2.2.2 graph(query)-graph(search)
                            # regular traversal of graph based on search_protein
                            logging.warning("CASE M2.2.2: graph vs graph.")
                            tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, query_path)
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
                    
                    best_entire_paths[protein]["search_protein"] = []
                    best_entire_paths[protein]["query_protein"] = []
                    best_entire_fixtures[protein] = max_fixture

                    for feature in max_fixture[0]:
                        if feature in search_features:
                            logging.debug(str(search_features[feature][0])+" "+str(search_features[feature][1])+" "+str(search_features[feature][2])+" "+str(search_features[feature][3]))
                            
                            best_entire_paths[protein]["search_protein"].append((search_features[feature][0], search_features[feature][1], search_features[feature][2], search_features[feature][3]))
                        else:
                            best_entire_paths[protein]["search_protein"].append((a_s_f[feature][0], a_s_f[feature][1], a_s_f[feature][2], a_s_f[feature][3]))
                            
                    for feature in max_fixture[2]:
                        if feature in query_features:
                            best_entire_paths[protein]["query_protein"].append((query_features[feature][0], query_features[feature][1], query_features[feature][2], query_features[feature][3]))
                        else:
                            best_entire_paths[protein]["query_protein"].append((a_q_f[feature][0], a_q_f[feature][1], a_q_f[feature][2], a_q_f[feature][3]))

                    logging.warning("Best_entire_searchpaths: "+str(best_entire_paths[protein]["search_protein"])+"\n\tProtein: "+protein)
                    logging.warning("Best_entire_querypaths: "+str(best_entire_paths[protein]["query_protein"])+"\n\tProtein: "+protein)
                   
            
            ################################################################
            ## classic mode ## only one protein will be linearized #########
            elif entire == 0:
		# erroneous
                logging.info("False predictions in mode entire=0. Not recommended. Please use mode entire=1. Quitting.")
                quit()
                                
        #score model M1 - used for final evaluation
	elif one_against_all == 1:
                search_features = {}    # set in su_search_protein
                a_s_f = {}              # set in su_search_protein
		search_protein = su_search_protein(single_protein[0], 1)                
                logging.info("search_protein(s): "+str(search_protein))
                                
                ####### entire mode ############################################
                # both proteins (search and query) will be linearized ##########
                if entire == 1:
                    if MS_uni == 0:
                        if weight_const == True:
                            w_weighting_constraints("none", entire) 
                        else:
                            logging.info("Weighting without constraints.")
                            w_weighting("none", entire)

                    #### SINGLE(search) --VS--> SET(query) ##
                    #set up a var for search graph complexity
                    tmp_path_set_size = 0
                    traversed = False
                    for protein in proteome:
                        logging.debug("protein:"+str(protein))
                        if MS_uni == 0:
                            if weight_const == True:
                                w_weighting_constraints_counter(protein, entire)
                            else:
                                w_weighting_counterpart(protein, entire)
                        # query_protein setup
                        mode[protein] = 0
                        lin_query_protein_set = su_lin_query_protein(protein, 1)
                        tmp_protein_graph = pb_entire_main_graph(lin_query_protein_set)
                        priority_check = False
                        
                        # check for size of query_protein, priority to reduce number of paths
                        if (int(len(query_features)) > int(priority_threshold)): 
                            logging.info("checking number of instances for query: "+str(protein))
                            priority_check = True
                            
                        else:
                            logging.info("checking cardinality of paths for query: "+str(protein))
                            #PRIORITY CHECK "2":
                            tmp_path_set_protein = pb_graphtraversal(tmp_protein_graph, 1)
                            if(int(len(tmp_path_set_protein)) > int(max_cardinality)):
                                priority_check = True

                                    
                        if priority_check == False:
                            #exhausitive case (priority_check: false)
                            #reuse of path_set from priority check above
                            all_protein_paths = tmp_path_set_protein #pb_graphtraversal(tmp_protein_graph, 1)
                            logging.debug("all_protein_paths:" +str(all_protein_paths))
                            logging.info("No. of all_protein_paths: "+ str(len(all_protein_paths)))
                        else:
                            #priority case (priority_check: true)
                            logging.info("priority graph traversal for query: "+str(protein))
                            all_protein_paths = []
                            for domain_type in proteome[protein]:
                                if proteome[protein][domain_type][0] == True:
                                    all_protein_paths += (pb_entire_graphtraversal_priority(tmp_protein_graph, domain_type, "none", 1))
                                    logging.info("domain_type: "+str(domain_type))
                                    logging.debug("all_protein_paths: "+str(all_protein_paths))
                        
                        pathcount = 0

                        # max fixture of (search_path, score, query_path, protenID)
                        max_fixture = ([], (0.0, 0.0, 0.0, 0.0, 0.0, False), [], protein)
                        
                        # check for available paths
                        # case M1.1: empty(query)
                        if len(all_protein_paths) == 0:
                            logging.warning("No paths (pfam or smart annotations) in query protein (set): "+str(protein))
                            
                            #case M1.1.1: empty(search)-empty(query)
                            # should be the best fix independent from weight
                            if int(len(search_features)) == 0:
                                score_w = sf_entire_calc_score(a_s_f.keys(),a_q_f.keys())
                                path = a_s_f.keys()
                                query_architecture = a_q_f.keys()
                                mode[protein] = 2
                            else:
                                # case M1.1.2: graph(search)-empty(query)   
                                if search_graph == 0:
                                    tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, [])
                                    path = tmp_path_score[0][0]
                                    score_w = tmp_path_score[0][1]
                                    query_architecture = tmp_path_score[0][2]
                                    mode[protein] = tmp_path_score[1]
                                    tmp_path_score = 0

                                elif (int(len(search_features)) > int(priority_threshold)):
                                    mode[protein] = 1
                                    tmp_path_score = pb_entire_priority_mode(single_protein[0], [])
                                    path = tmp_path_score[0]
                                    score_w = tmp_path_score[1]
                                    query_architecture = tmp_path_score[2]
                                    logging.debug("path, score_w, query_architecture, mode: "+str(path)+", "+str(score_w)+", "+str(query_architecture)+", "+str(mode))

                                else:
                                    mode[protein] = 0
                                    path_score = pb_entire_graphtraversal(search_graph, [], 0)
                                    path = path_score[0]
                                    score_w = path_score[1]
                                    query_architecture = path_score[2]
                                    logging.debug("path, score_w, mode: "+str(path)+", "+str(score_w)+", "+str(mode) )
                                
                            # set max_fixture according to the previous settings
                            max_fixture = (path, score_w, query_architecture, protein)
                            
                        ## handle each query_path ##
                        ## case 2: graph(query) ##
                        
                                                
                        for query_path in all_protein_paths:
                            logging.info("\n\tProtein and Query_path: "+str(protein)+": " +str(query_path))
                            logging.debug("search_protein graph: "+str(search_graph))
                            logging.debug("query_protein graph: "+str(query_graph))
                            
                            pathcount += 1
                                                        
                            # case M1.2.1: empty(search)-graph(query)
                            if int(len(search_features)) == 0:
                                query_path_ad = query_path+a_q_f.keys()
                                path = a_s_f.keys()
                                score_w = sf_entire_calc_score(path,query_path_ad)
                            else:
                                # caser M1.2.2 graph(search)-graph(query)
                                if search_graph == 0:
                                    tmp_path_score = pb_entire_main_nongreedy(search_protein, protein, query_path)
                                    path = tmp_path_score[0][0]
                                    score_w = tmp_path_score[0][1]
                                    query_path_ad = tmp_path_score[0][2]
                                    mode[protein] = tmp_path_score[1]
                                    tmp_path_score = 0
                                    logging.debug("path, score_w, mode: "+str(path)+", "+str(score_w)+", "+str(mode) )

                                    # check for max scoring fixture of path and query_path
                                    if score_w[5] == True:
                                        if score_w[3] >= max_fixture[1][3]:
                                            max_fixture = (path, score_w, query_path_ad, protein)
                                            logging.warning("max_fixture: "+str(max_fixture))
                                    else:
                                        if (max_fixture[1][5] != True) and (score_w[4] >= max_fixture[1][4]):
                                            max_fixture = (path, score_w, query_path_ad, protein)

                                else:
                                    if (int(len(search_features)) > int(priority_threshold)):
                                        # priotity check 1 (quantity)
                                        priority_check = True
                                        mode[protein] = 1

                                    else:
                                        if (traversed == False):
                                            tmp_path_set_size = int(len(pb_graphtraversal(search_graph, 1)))
                                            traversed = True
                                        if(int(tmp_path_set_size) > int(max_cardinality)):
                                            # priority check 2 (complexity)
                                            priority_check=True
                                            mode[protein]=1

                                    if priority_check == False:
                                        # exhaustive case
                                        mode[protein] = 0
                                        path_score = pb_entire_graphtraversal(search_graph, query_path, 0)
                                        path = path_score[0]
                                        score_w = path_score[1]
                                        query_path_ad = path_score[2]
                                        logging.debug("path, score_w, mode: "+str(path)+", "+str(score_w)+", "+str(mode) )
                                    else:
                                        # priority case
                                        logging.info("go priority mode for "+str(protein))

                                        path_score = pb_entire_priority_mode(single_protein[0], query_path)
                                        path = path_score[0]
                                        score_w = path_score[1]
                                        query_path_ad = path_score[2]
                                        logging.debug("path, score_w, query_path_ad, mode: "+str(path)+", "+str(score_w)+", "+str(query_path_ad)+", "+str(mode))

                                        #check for max scoring fixture of path and query_path
                                        if score_w[5] == True:
                                            if score_w[3] >= max_fixture[1][3]:
                                                max_fixture = (path, score_w, query_path_ad, protein)
                                                logging.warning("max_fixture: "+str(max_fixture))
                                        else:
                                            if max_fixture[1][5] != True and score_w[4] >= max_fixture[1][4]:
                                                max_fixture = (path, score_w, query_path_ad, protein)                                    
                                
                            #check for max scoring fixture of path and query_path
                            if score_w[5] == True:
                                if score_w[3] >= max_fixture[1][3]:
                                    max_fixture = (path, score_w, query_path_ad, protein)
                                    logging.warning("max_fixture: "+str(max_fixture))
                            else:
                                if (max_fixture[1][5] != True) and (score_w[4] >= max_fixture[1][4]):
                                    max_fixture = (path, score_w, query_path_ad, protein)
                                
                            
                        logging.info("Found: "+str(pathcount)+" path(s) for query.")
                        logging.debug("Path max_fixture: "+str(max_fixture))
                        
                        best_entire_paths[protein]["search_protein"] = []
                        best_entire_paths[protein]["query_protein"] = []
                        best_entire_fixtures[protein] = max_fixture
                        
                        for feature in max_fixture[0]:
                            if feature in search_features:
                                logging.debug(str(search_features[feature][0])+" "+str(search_features[feature][1])+" "+str(search_features[feature][2])+" "+str(search_features[feature][3]))
                                best_entire_paths[protein]["search_protein"].append((search_features[feature][0], search_features[feature][1], search_features[feature][2], search_features[feature][3]))
                            else:
                                logging.debug(str(a_s_f[feature][0])+" "+str(a_s_f[feature][1])+" "+str(a_s_f[feature][2])+" "+str(a_s_f[feature][3]))
                                best_entire_paths[protein]["search_protein"].append((a_s_f[feature][0], a_s_f[feature][1], a_s_f[feature][2], a_s_f[feature][3]))
                                
                        for feature in max_fixture[2]:
                            if feature in query_features:
                                best_entire_paths[protein]["query_protein"].append((query_features[feature][0], query_features[feature][1], query_features[feature][2], query_features[feature][3]))
                            else:
                                best_entire_paths[protein]["query_protein"].append((a_q_f[feature][0], a_q_f[feature][1], a_q_f[feature][2], a_q_f[feature][3]))

                        logging.warning("Best_entire_searchpaths: "+str(best_entire_paths[protein]["search_protein"]))
                        logging.warning("Best_entire_querypaths: "+str(best_entire_paths[protein]["query_protein"]))
                        
                ################################################################
                ## classic mode ## only one protein will be linearized #########
                elif entire == 0:
                    # erroneous
                    logging.info("False predictions in mode entire=0. Not recommended. Please use mode entire=1. Quitting.")
                    quit()


# fas score calculation after best paths for pfam/smart where created
def fc_normal():
	global one_against_all
	global proteome
	global best_paths
        global best_entire_paths
        global best_entire_fixtures#Baustelle: TBD
	global search_features
        global query_features
	global outpath
	global single_protein
	global protein_lengths
	global MS_uni
	global weights
	global output
	global mode
        global weight_const
        
            
        # output for single protein (xml)
	if output == 0 or output == 2:
            if one_against_all == 1:
                tmp = "single-->set"
            else:
                tmp = "set-->single"
            out = open(outpath, "w+")
            out.write("<?xml version=\"1.0\"?>\n")
            if MS_uni == 0:
                out.write("<out direction=\"" + tmp + "\" weighting=\"applied\">\n")
            else:
                out.write("<out direction=\"" + tmp + "\" weighting=\"uniform\">\n")
            out.write("\t<single_protein id=\"" + single_protein[0] + "\" length=\"" + str(int(protein_lengths["single_"+str(single_protein[0])])) + "\">\n")

        elif(taciturn == 1):
            out = open(outpath, "w+")
	if one_against_all == 1:
            #w_weighting("none", entire)
            if output == 0 or output == 2:
                for feature in single_protein[1]:

                    if MS_uni == 0:
                        out.write("\t\t<feature type=\"" + feature + "\" evalue=\""+str(single_protein[1][feature][1])+"\" weight=\"" + str(weights[feature]) + "\">\n")
                    else:
                        out.write("\t\t<feature type=\"" + feature + "\" evalue=\""+str(single_protein[1][feature][1])+"\" weight=\"" + str(1.0/len(single_protein[1])) + "\">\n")

                    for instance in single_protein[1][feature][2:]:
                        out.write("\t\t\t<instance inst_eval=\""+str(instance[0])+"\" start=\"" + str(instance[1]) + "\" end=\"" + str(instance[2]) + "\"/>\n")

                    out.write("\t\t</feature>\n")

	else:
            if output == 0 or output == 2:

                for feature in single_protein[1]:

                    if one_against_all == 0:
                        if MS_uni == 0:
                            out.write("\t\t<feature type=\"" + feature + "\" evalue=\""+str(single_protein[1][feature][1])+"\" weight=\"" + str(weights_counter[feature]) + "\">\n") #hiernochcounterweight
                        else:
                            out.write("\t\t<feature type=\"" + feature + "\" evalue=\""+str(single_protein[1][feature][1])+"\" weight=\"" + str(1.0/len(single_protein[1])) + "\">\n")

                    for instance in single_protein[1][feature][2:]:
                        out.write("\t\t\t<instance inst_eval=\""+str(instance[0])+"\" start=\"" + str(instance[1]) + "\" end=\"" + str(instance[2]) + "\"/>\n")

                    out.write("\t\t</feature>\n")

	if output == 0 or output == 2:
            out.write("\t</single_protein>\n")
	if one_against_all == 0:
            su_query_protein(single_protein[0])
        
        # selecting paths for each protein         
	for protein in proteome:
                logging.debug("protein: "+str(protein))

                # setup for standard and entire mode
                if entire == 0:
                    # erroneous
                    logging.info("False predictions in mode entire=0. Not recommended. Quitting.")
                    quit()
                        
                elif entire == 1:
                    logging.warning("Best_entire_searchpaths: "+str(best_entire_paths[protein]["search_protein"])+"\n\tProtein: "+protein)
                    logging.warning("Best_entire_querypaths: "+str(best_entire_paths[protein]["query_protein"])+"\n\tProtein: "+protein)
                    logging.warning("Best_entire_fixtures: "+str(best_entire_fixtures[protein])+"\n\tProtein: "+protein)

		if entire == 0:
                    # erroneous
                    logging.info("False predictions in mode entire=0. Not recommended. Quitting.")
                    quit()
                elif entire == 1:
                    score = best_entire_fixtures[protein][1]
                    
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
                    out.write("\t<set_protein id=\"" + protein + "\" score=\"" + str(score[3]) + "\" MS=\"" + str(score[0]) + "\" PS=\"" + str(score[1]) + "\" CS=\"" + str(score[2]) + "\" length=\"" + str(int(protein_lengths["set_"+str(protein)])) + "\" mode=\"" + mode_out + "\">\n")
                    out.write("\t\t<architecture>\n")
                    
                    for feature in proteome[protein]:

                        if one_against_all == 0 and MS_uni == 0:
                            out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\""+str(proteome[protein][feature][1])+"\" weight=\"" + str(weights[feature]) + "\">\n")
                        elif one_against_all == 1 and MS_uni == 0:
                            out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\""+str(proteome[protein][feature][1])+"\" weight=\"" + str(weights_counter[feature]) + "\">\n")
                        else:
                            out.write("\t\t\t<feature type=\"" + feature + "\" evalue=\""+str(proteome[protein][feature][1])+"\" weight=\"" + str(1.0/len(proteome[protein])) + "\">\n")
                        for instance in proteome[protein][feature][2:]:
                            out.write("\t\t\t\t<instance inst_eval=\""+str(instance[0])+"\" start=\"" + str(instance[1]) + "\" end=\"" + str(instance[2]) + "\"/>\n")
                        out.write("\t\t\t</feature>\n")

                    out.write("\t\t</architecture>\n")
                    
                    ### handling feature paths    
                    out.write("\t\t<path>\n")

                    # rescaling
                    if entire == 1:
                        logging.debug("best_entire_paths[protein][search_protein]: "+ str(best_entire_paths[protein]["search_protein"]))
                        logging.debug("best_entire_paths[protein][query_protein]: "+ str(best_entire_paths[protein]["query_protein"]))
                    elif entire == 0:
                        # erroneous
                        logging.info("False predictions in mode entire=0. Not recommended. Quitting.")
                        quit()
                    
                    path_tmp = {}
                    scale = 0
                    if entire == 0:
                        # erroneous
                        logging.info("False predictions in mode entire=0. Not recommended. Quitting.")
                        quit()
                    elif entire == 1:
                        if weight_const == 1:
                            path_tmp2 = []
                            for feature in best_entire_paths[protein]["search_protein"]:
                                if feature[0] not in path_tmp2:
                                    path_tmp2.append(feature[0])
                            adjusted_weights = w_weight_const_rescale_final(path_tmp2)
                            weight_tmp = {}
                            for adj_feature in adjusted_weights:
                                weight_tmp[adj_feature] = weights[adj_feature]
                                weights[adj_feature] = adjusted_weights[adj_feature]
                        for feature in best_entire_paths[protein]["search_protein"]:
                            if feature[0] in path_tmp:
                                path_tmp[feature[0]].append((feature[2],feature[3]))
                            else:
                                path_tmp[feature[0]] = [(feature[2],feature[3])]
                                if MS_uni == 0:				
                                    scale += weights[feature[0]]
                        logging.debug("path_tmp: "+ str(path_tmp))


                    ## unweighted case
                    if MS_uni == 0:	
                        if scale > 0:
                            scale = 1.0/float(scale)
                        else:
                            scale = 1.0
                    ## print path
                    for feature in path_tmp:
                        if MS_uni == 0:
                            out.write("\t\t\t<feature type=\"" + feature + "\" corrected_weight=\"" + str(weights[feature]*scale) + "\">\n")
                        else:
                            out.write("\t\t\t<feature type=\"" + feature + "\">\n")
                        for tmp_inst in path_tmp[feature]:
                            out.write("\t\t\t\t<instance start=\"" + str(tmp_inst[0]) + "\" end=\"" + str(tmp_inst[1]) + "\"/>\n")
                        out.write("\t\t\t</feature>\n")
                    out.write("\t\t</path>\n")
                    out.write("\t</set_protein>\n")
                    if weight_const == 1:
                         for adj_feature in adjusted_weights:
                            weights[adj_feature] = weight_tmp[adj_feature]
		if output == 1 or output == 2:
                    print score[3]
                    ## hidden ##
                    if (taciturn == 1 and output != 2):
                        out.write(protein+"\t"+str(score[3])+"\n")
	if output == 0 or output == 2:
		out.write("</out>")
		out.close()


########## Start Up Functions ########## <su>

# define of 2-dim-dictionary
def tree():
    return collections.defaultdict(tree)

# initializes query protein
def su_query_protein(protein_id):
    global one_against_all
    global proteome
    global single_protein
    global query_protein
    global query_clans
    global protein_lengths
    global clan_dict

    query_protein = {}
    query_clans = {}	
    if one_against_all == 0:
        for feature in single_protein[1]:
            query_protein[feature] = []
            clan = clan_dict[feature]
            if clan in query_clans:
                query_clans[clan] += len(single_protein[1][feature])-2
            else:
                query_clans[clan] = len(single_protein[1][feature])-2
            for instance in single_protein[1][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["single_"+str(protein_id)])
                position = round(position, 8)
                query_protein[feature].append(position)
                
    else:
        for feature in proteome[protein_id]:
            query_protein[feature] = []
            clan = clan_dict[feature]
            if clan in query_clans:
                query_clans[clan] += len(proteome[protein_id][feature])-2
            else:
                query_clans[clan] = len(proteome[protein_id][feature])-2
            for instance in proteome[protein_id][feature][2:]:
                position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["set_"+str(protein_id)])
                position = round(position, 8)
                #features and relative positions
                query_protein[feature].append(position)
                #{'pfam_RRM_occluded': [0.40105263, 0.62631579], 'smart_RRM': [0.39684211, 0.62421053, 0.90105263], 'pfam_RRM_1': [0.39578947, 0.62315789, 0.90421053], 'pfam_RRM_7': [0.37894737, 0.60315789], 'pfam_RRM_6': [0.39578947, 0.62210526, 0.90210526], 'pfam_RRM_5': [0.43263158, 0.64210526, 0.91157895], 'smart_RRM_1': [0.43473684, 0.62421053, 0.88315789], 'pfam_Transformer': [0.14526316]}

# initializes the search protein
def su_lin_query_protein(protein_id, pfam):
    global one_against_all
    global proteome
    global single_protein
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
    if one_against_all == 0:
	for feature in single_protein[1]:
            clan = clan_dict[feature]
            if clan in query_clans:
                query_clans[clan] += len(single_protein[1][feature])-2
            else:
                query_clans[clan] = len(single_protein[1][feature])-2

            if single_protein[1][feature][0] == True:
                for instance in single_protein[1][feature][2:]:
                    position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["single_"+str(protein_id)])
                    position = round(position, 8)
                    key = "F_" + str(i)
                    query_features[key] = (feature, position, instance[1], instance[2])
                    #{'F_16': ('smart_RRM_1', 0.43473684, 186, 227), 'F_10': ('pfam_RRM_6', 0.39578947, 151, 225), 'F_17': ('smart_RRM_1', 0.62421053, 261, 332), 'F_11': ('pfam_RRM_6', 0.62210526, 261, 330), 'F_14': ('pfam_RRM_5', 0.64210526, 275, 335), 'F_4': ('smart_RRM', 0.90105263, 394, 462), 'F_5': ('pfam_RRM_1', 0.39578947, 151, 225), 'F_6': ('pfam_RRM_1', 0.62315789, 261, 331), 'F_7': ('pfam_RRM_1', 0.90421053, 399, 460), 'F_0': ('pfam_RRM_occluded', 0.40105263, 148, 233), 'F_1': ('pfam_RRM_occluded', 0.62631579, 259, 336), 'F_2': ('smart_RRM', 0.39684211, 150, 227), 'F_3': ('smart_RRM', 0.62421053, 260, 333), 'F_15': ('pfam_RRM_5', 0.91157895, 402, 464), 'F_12': ('pfam_RRM_6', 0.90210526, 398, 459), 'F_19': ('pfam_Transformer', 0.14526316, 6, 132), 'F_8': ('pfam_RRM_7', 0.37894737, 148, 212), 'F_9': ('pfam_RRM_7', 0.60315789, 258, 315), 'F_18': ('smart_RRM_1', 0.88315789, 378, 461), 'F_13': ('pfam_RRM_5', 0.43263158, 183, 228)}
                    tmp.append((key, instance[1]))
                    i += 1	
            else:
                for instance in single_protein[1][feature][2:]:
                    position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["single_"+str(protein_id)])
                    position = round(position, 8)
                    key = "O_" + str(i)					
                    a_q_f[key] = (feature, position, instance[1], instance[2])
                    i += 1	
    elif one_against_all == 1:
        for feature in proteome[protein_id]:
            clan = clan_dict[feature]
            if clan in query_clans:
                query_clans[clan] += len(proteome[protein_id][feature])-2
            else:
                query_clans[clan] = len(proteome[protein_id][feature])-2
            if proteome[protein_id][feature][0] == True:
                for instance in proteome[protein_id][feature][2:]:				
                    position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["set_"+str(protein_id)])
                    position = round(position, 8)
                    key = "F_" + str(i)
                    query_features[key] = (feature, position, instance[1], instance[2])
                    tmp.append((key, instance[1]))
                    i += 1	
            else:
                for instance in proteome[protein_id][feature][2:]:
                    position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["set_"+str(protein_id)])
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
def su_search_protein(protein_id, pfam):
	global one_against_all
	global proteome
	global single_protein
	global search_features
        global a_s_f

	search_protein = []
	tmp = []
	i = 0
	if one_against_all == 0:
            for feature in proteome[protein_id]:
                # True if feature is going to be linearized
                if proteome[protein_id][feature][0] == True:
                    for instance in proteome[protein_id][feature][2:]:				
                        position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["set_"+str(protein_id)])
                        position = round(position, 8)
                        key = "F_" + str(i)
                        search_features[key] = (feature, position, instance[1], instance[2])
                        tmp.append((key, instance[1]))
                        i += 1
                else:
                    for instance in proteome[protein_id][feature][2:]:				
                        position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["set_"+str(protein_id)])
                        position = round(position, 8)
                        key = "O_" + str(i)
                        a_s_f[key] = (feature, position, instance[1], instance[2])
                        i += 1
                                			
	else:
            for feature in single_protein[1]:
                if single_protein[1][feature][0] == True:
                    for instance in single_protein[1][feature][2:]:
                        position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["single_"+str(protein_id)])
                        position = round(position, 8)
                        key = "F_" + str(i)
                        search_features[key] = (feature, position, instance[1], instance[2])
                        #{'F_16': ('smart_RRM_1', 0.43473684, 186, 227), 'F_10': ('pfam_RRM_6', 0.39578947, 151, 225), 'F_17': ('smart_RRM_1', 0.62421053, 261, 332), 'F_11': ('pfam_RRM_6', 0.62210526, 261, 330), 'F_14': ('pfam_RRM_5', 0.64210526, 275, 335), 'F_4': ('smart_RRM', 0.90105263, 394, 462), 'F_5': ('pfam_RRM_1', 0.39578947, 151, 225), 'F_6': ('pfam_RRM_1', 0.62315789, 261, 331), 'F_7': ('pfam_RRM_1', 0.90421053, 399, 460), 'F_0': ('pfam_RRM_occluded', 0.40105263, 148, 233), 'F_1': ('pfam_RRM_occluded', 0.62631579, 259, 336), 'F_2': ('smart_RRM', 0.39684211, 150, 227), 'F_3': ('smart_RRM', 0.62421053, 260, 333), 'F_15': ('pfam_RRM_5', 0.91157895, 402, 464), 'F_12': ('pfam_RRM_6', 0.90210526, 398, 459), 'F_19': ('pfam_Transformer', 0.14526316, 6, 132), 'F_8': ('pfam_RRM_7', 0.37894737, 148, 212), 'F_9': ('pfam_RRM_7', 0.60315789, 258, 315), 'F_18': ('smart_RRM_1', 0.88315789, 378, 461), 'F_13': ('pfam_RRM_5', 0.43263158, 183, 228)}
                        tmp.append((key, instance[1]))
                        i += 1	
                else:
                    for instance in single_protein[1][feature][2:]:
                        position = ((float(instance[1]) + float(instance[2])) / 2.0) / float(protein_lengths["single_"+str(protein_id)])
                        position = round(position, 8)
                        key = "O_" + str(i)
                        a_s_f[key] = (feature, position, instance[1], instance[2])
                        i += 1	
        #sort  instances
	tmp2 = sorted(tmp, key=itemgetter(1))
	for x in tmp2:
		search_protein.append(x[0])
                
	return search_protein



# adds best pfam/smart path to the other features
def su_add_path(search_protein, protein_id):
	global search_features
	global best_paths

	i = 0
	for feature in best_paths[protein_id]:
		key = "P_" + str(i)
		search_features[key] = (feature[1], feature[2])
		search_protein.append(key)
		i += 1
	return search_protein

# adds best pfam/smart path to other features in entire mode
def su_entire_add_searchpath(search_protein, protein_id):
    global search_features
    global best_entire_paths
    logging.info("su_entire_add_searchpath")
    logging.debug("search_protein: "+str(search_protein)+" in su_entire_add_searchpath")
    logging.debug("entire_paths: "+str(best_entire_paths[protein_id]["search_protein"])+" in su_entire_add_searchpath")
    
    i = 0
    for feature in best_entire_paths[protein_id]["search_protein"]:
        key = "P_" + str(i)
        search_features[key] = (feature[1], feature[2])
        search_protein.append(key)
        i += 1
    
    return search_protein

# adds best pfam/smart path to other features for the query protein in entire mode
def su_entire_add_querypath(query_protein, protein_id):
    global query_features
    global best_entire_paths
    i = 0
    for feature in best_entire_paths[protein_id]["query_protein"]:
        key = "P_" + str(i)
        query_features[key] = (feature[1], feature[2])
        query_protein.append(key)
        i += 1
        
    return query_protein

# gives path for output
def su_set_path(jobname):
	global expath

	if os.path.exists(expath + "/out/" + jobname + ".xml"):
		i = 1
		while os.path.exists(expath + "/out/" + jobname + "_" + str(i) + ".xml"):
			i += 1
		jobname = jobname+"_"+str(i)
	return jobname
			 


########## Pathbuilding Functions ########## <pb>
# Used for the Pfam/Smart domains 

# Checks if best path can be found using Greedy Strategy
def pb_greedycheck(search_protein, protein):
	global single_protein
	global proteome
	global one_against_all
	global search_features
        
        logging.info("function pb_greedycheck")

	multi_features = {}
	greedy = True
	i = 0
	if one_against_all == 0:
		for feature in proteome[protein]:
			if len(proteome[protein][feature]) > 3:
                                # substract first 2 entry from length feature -> (True, Evalue, instance1, instance2)
				multi_features[feature] = len(proteome[protein][feature]) - 2
	else:
		for feature in single_protein[1]:
			if len(single_protein[1][feature]) > 3:
                                logging.info("Number of instances: "+str(feature)+" "+str(len(single_protein[1][feature])-2))
                                multi_features[feature] = len(single_protein[1][feature]) - 2
	## check for general feasibility of greedy strategy
        ## are multi_features contained in any overlap region?
        while i<len(search_protein) and len(multi_features) > 0:
                logging.debug("greedycheck while loop 1")
                    
                ## Baustelle ##
		if i+1 == len(search_protein):
                        logging.debug("breaks in greedyCheck")
			break
		elif search_features[search_protein[i]][3] >= search_features[search_protein[i+1]][2]:
			if search_features[search_protein[i]][3] < search_features[search_protein[i+1]][3]:
				end = search_features[search_protein[i+1]][3]
			else:
				end = search_features[search_protein[i]][3];
                        if search_features[search_protein[i]][0] in multi_features:
				if len(multi_features) == 1 and multi_features[search_features[search_protein[i]][0]] <= 1:
                                        logging.debug("breaks in greedyCheck")
					break
				else:
					greedy = False
					break
                        if search_features[search_protein[i+1]][0] in multi_features:
				if len(multi_features) == 1 and multi_features[search_features[search_protein[i+1]][0]] <= 1:
                                        logging.debug("breaks in greedyCheck")
					break
				else:
					greedy = False
					break
			i += 2
			while (i<len(search_protein)) and (search_features[search_protein[i]][2] <= end):
                                logging.debug("greedyCheck while loop 2")
				if search_features[search_protein[i]][0] in multi_features:
					if len(multi_features) == 1 and multi_features[search_features[search_protein[i]][0]] <= 1:
                                                logging.debug("breaks in greedyCheck")
						break
					else:
						greedy = False
                                                logging.debug("breaks in greedyCheck")
                                                break
                                if end < search_features[search_protein[i]][3]:
					end = search_features[search_protein[i]][3]
				i += 1
		else:
			if search_features[search_protein[i]][0] in multi_features:
                                logging.info("reduce multifeature "+ str(search_features[search_protein[i]][0]))
                                    
				multi_features[search_features[search_protein[i]][0]] = multi_features[search_features[search_protein[i]][0]] - 1
				if multi_features[search_features[search_protein[i]][0]] == 0:
					del multi_features[search_features[search_protein[i]][0]]
                                        logging.info("deletion of: "+str(multi_features[search_features[search_protein[i]][0]]))                                        
			i += 1
        ## further check for size of largest overlap region
        logging.info("check size of overlap regions")
        overlap_control = True
        while i<len(search_protein) and overlap_control:
		current = search_protein[i]
		if i+1 == len(search_protein):
			break
		elif search_features[current][3] >= search_features[search_protein[i+1]][2]:
			if search_features[current][3] < search_features[search_protein[i+1]][3]:
				end = search_features[search_protein[i+1]][3]
			else:
				end = search_features[current][3];
   
			overlap_region_check = [current, search_protein[i+1]]
			i += 2
			while (i<len(search_protein)) and (search_features[search_protein[i]][2] <= end):
				overlap_region_check.append(search_protein[i])
                                if int(len(overlap_region_check)) >= int(priority_threshold):
                                    logging.debug("overlap control: " + str(len(overlap_region_check)))
                                    logging.info("breaking --> no greedy")
                                    greedy = False
                                    overlap_control = False
                                    break
				if end < search_features[search_protein[i]][3]:
					end = search_features[search_protein[i]][3]
				i += 1
                        logging.debug("overlap region check: " + str(overlap_region_check) + " of size " + str(len(overlap_region_check)))
                else:
                    i += 1
                        
        logging.info("greedy set to: " + str(greedy) + " in function pb_greedycheck")
        return greedy

# build graph for a given protein
def pb_entire_main_graph(current_protein):
    global query_graph
    region = pb_query_region_mapper(current_protein)
    query_graph = pb_region_paths_nongreedy(region)
    return query_graph
    
# creates graph with all paths 
# retrieves best path from graph traversal function
# returns best path, score and mode
# including priority check
def pb_entire_main_nongreedy(search_protein, protein_id, query_path):
    global search_graph
    global max_cardinality
    logging.info("pb_entire_main_nongreedy")
    
    priority_check = False
        
    region = pb_search_region_mapper(search_protein)
    search_graph = pb_region_paths_nongreedy(region)
    logging.debug(region)
    logging.debug(search_graph)
    
    if ((int(len(search_features)) >= int(priority_threshold))):
        mode = 1
        priority_check = True
        # checking best path for all domain types
        path_score = pb_entire_priority_mode(protein_id, query_path)
    else:
        #PRIORITY CHECK "2": for every protein in proteome
        tmp_path_set_size = len(pb_graphtraversal(search_graph, 1))
        logging.debug("cardinality of tmp graph: "+str(tmp_path_set_size))
        if ((int(tmp_path_set_size) > int(max_cardinality))):
            mode = 1
            priority_check = True
            path_score = pb_entire_priority_mode(protein_id, query_path)
        
    if (priority_check == False):
        mode = 0
        path_score = pb_entire_graphtraversal(search_graph, query_path, mode)
        #path = path_score[0]
    if one_against_all == 0:
        search_graph == 0
    return (path_score, mode)

# non-Greedy version 
# creates a graph which contains all paths
def pb_main_nongreedy(current_protein, protein_id):
	global search_features
	global search_graph
	global one_against_all
	global priority_threshold
        
        logging.info("pb_main_nongreedy")

	region = pb_search_region_mapper(current_protein)
	search_graph = pb_region_paths_nongreedy(region)
                
	if (int(len(search_features)) >= int(priority_threshold)):
		mode = 1
		path = pb_priority_mode(protein_id)
	else:
		mode = 0
		path = pb_graphtraversal(search_graph, 0)
	if one_against_all == 0:
		search_graph == 0
	return (path, mode)

#priority mode to reduce paths in entire mode
def pb_entire_priority_mode(protein, query_path):
    global search_graph
    global one_against_all
    global single_protein
    global proteome
    
    logging.info("pb_entire_priority_mode")
    
    # best_path (Path, (MS_score(float), PS_score(float), CS_score(float), final_score(float), path_weight(float), common_feature(bool)), QPath)
    best_path = ("NULL", (0.0, 0.0, 0.0, 0.0, 0.0, False),"NULL")
    if one_against_all == 0:
        for domain_type in proteome[protein]:
            if proteome[protein][domain_type][0] == True: 
                path_score_w = pb_entire_graphtraversal_priority(search_graph, domain_type, query_path, 0)
                if path_score_w[1][5] == True:
                    if path_score_w[1][3] >= best_path[1][3]:
                        best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
                else:
                    #if so far no path with common features found AND the weight is higher
                    if ((best_path[1][5] != True) and (path_score_w[1][4] >= best_path[1][4])):
                        best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
                    
    else:
        for domain_type in single_protein[1]:
            if single_protein[1][domain_type][0] == True:
                path_score_w = pb_entire_graphtraversal_priority(search_graph, domain_type, query_path, 0)
                if path_score_w[1][5] == True:
                    if path_score_w[1][3] >= best_path[1][3]:
                        best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
                else:
                    if ((best_path[1][5] != True) and (path_score_w[1][4] >= best_path[1][4])):
                        best_path = (path_score_w[0], path_score_w[1], path_score_w[2])
    return best_path    

# priority mode to reduce calculated paths
def pb_priority_mode(protein):
	global search_graph
	global one_against_all
	global single_protein
	global proteome
	
	best_path = ("NULL", (0.0, 0.0, 0.0, 0.0))
	if one_against_all == 0:
		for domain_type in proteome[protein]:
			path = pb_graphtraversal_priority(search_graph, domain_type)
			if path[1][3] >= best_path[1][3]:
				best_path = (path[0], path[1])
	else:
		for domain_type in single_protein[1]:
			path = pb_graphtraversal_priority(search_graph, domain_type)
			if path[1][3] >= best_path[1][3]:
				best_path = (path[0], path[1])
	return best_path[0]
	
	

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



# Greedy version
# directly creates the best path 
# every path of an overlap region is processed immediatly
def pb_main(search_protein):
	global search_features
        
	i = 0
	path = []
	while i<len(search_protein):
		current = search_protein[i]
		if i+1 == len(search_protein):
			path.append(search_protein[i])
			break
		elif search_features[current][3] >= search_features[search_protein[i+1]][2]:
                            
			if search_features[current][3] < search_features[search_protein[i+1]][3]:
				end = search_features[search_protein[i+1]][3]
			else:
				end = search_features[current][3];

			overlap_region = [current, search_protein[i+1]]
			i += 2
			while (i<len(search_protein)) and (search_features[search_protein[i]][2] <= end):
				overlap_region.append(search_protein[i])
				if end < search_features[search_protein[i]][3]:
					end = search_features[search_protein[i]][3]
				i += 1
                        logging.debug("--> overlap region: "+str(overlap_region)+" size: "+str(int(len(overlap_region))))
			overlap_map = pb_search_region_mapper(overlap_region)
			overlap_paths = pb_region_paths(overlap_map)
                        logging.debug("--> overlap_paths: "+str(overlap_paths))
                    	best_path = ([], (0.0, 0.0, 0.0, 0.0))
			for overlap_path in overlap_paths:
				scored_path = path + overlap_path
                                score = sf_calc_score(scored_path, 0, "none")
				if score[3] >= best_path[1][3]:
					best_path = (scored_path, score)
			path = best_path[0]
		else:
                        
                        logging.debug("added to path"+str(search_protein[i]))
			path.append(search_protein[i])
        		i += 1
        logging.debug("pb_main returns: "+str(path))
	return (path)

# creates an overlap map for a given overlap_region for the query_features. That is for each feature a list of features which do not overlap with it.
def pb_query_region_mapper(overlap_region):
    global query_features
    logging.info("pb_query_region_mapper")
    logging.debug(overlap_region)

    overlap_map = [("START", overlap_region)]
    for i in range(0, len(overlap_region)):
        end = query_features[overlap_region[i]][3]
        overlap = []
	x = i+1
	while x < len(overlap_region):
            if end <= query_features[overlap_region[x]][2]:
                overlap.append(overlap_region[x])
            x = x+1
        overlap.append("END")
        overlap_map.append((overlap_region[i], overlap))
    overlap_map.append(("END", []))
    overlap_map.reverse()
    return overlap_map

# creates an overlap map for a given overlap_region for the search_features. That is for each feature a list of features which do not overlap with it.
def pb_search_region_mapper(overlap_region):
    global search_features
    logging.info("pb_search_region_mapper")
    logging.debug(overlap_region)

    overlap_map = [("START", overlap_region)]
    
    for i in range(0, len(overlap_region)):
        end = search_features[overlap_region[i]][3]
        overlap = []
	x = i+1
	while x < len(overlap_region):
            if end <= search_features[overlap_region[x]][2]:
                overlap.append(overlap_region[x])
            x = x+1
        overlap.append("END")
        overlap_map.append((overlap_region[i], overlap))
    overlap_map.append(("END", []))
    overlap_map.reverse()
    return overlap_map

# creates the different paths for an overlap region for the Greedy Strategy
def pb_region_paths(overlap_map):
    
        logging.info("pb_region_paths")
	search_graph = {}
	reached_list = {}
	for feature in overlap_map:
            reached = len(feature[1])
            rcount=0
            links = []
            for candidate in feature[1]:
                    links.append(candidate)
                    if (-reached_list[candidate]+rcount) < 0:
                        del feature[1][(-reached_list[candidate]+rcount):]
                    if rcount < reached_list[candidate]:			
                        rcount = reached_list[candidate]
            reached_list[feature[0]] = reached
            search_graph[feature[0]] = links
	paths = pb_graphtraversal(search_graph, 1)
	return paths

# return best scoring path 
def pb_entire_graphtraversal(search_graph, query_path, mode):
        global a_s_f
        global a_q_f
    
        logging.info("pb_entire_graphtraversal in mode: "+str(mode))
        logging.debug("search graph: "+str(search_graph))
    
    	v_stack = ["START"]
	p_stack = []
	paths = []
	if mode == 0:
            best_path = ([], (0.0, 0.0, 0.0, 0.0, 0.0, False),[])
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
                        score_w = sf_entire_calc_score(path_ad, query_path_ad)
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
# for non-greedy it returns the best scoring
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
                        score = sf_calc_score(path, mode, "none")
                        logging.info("search path "+str(path)+" in pb_graphtraversal")
                        logging.info("Score: "+str(score)+" in pb_graphtraversal")
                        if score[3] >= best_path[1][3]:
                            best_path = (path, score)
                else:
                    v_stack.append(next)
                    p_stack.append(path + [next])
	if mode == 1:
            logging.debug("returning paths: "+str(paths))
            logging.info("cardinality: "+str(len(paths)))
            return paths
	else:
            logging.info("returning "+str(best_path[0]))
            return best_path[0]
#
# walks the graph in priority-mode to reduce number of paths
def pb_entire_graphtraversal_priority(search_graph, priority, query_path, mode):
	global search_features
        global query_features
        global a_s_f
        global a_q_f
        
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
            if len(search_graph[vertex]) < 2: ## unsafe threshold ##
                for next in search_graph[vertex]:
                    if next == "END":
                        if mode == 1:
                            paths.append(path)
                        else:
                            path_ad = path + a_s_f.keys()
                            query_path_ad = query_path + a_q_f.keys()
                            score_w = sf_entire_calc_score(path_ad, query_path_ad)
                            # common feature (score_w[5] is true if two paths have at least one common feature (main feature)
                            if score_w[5] == True:
                                if score_w[3] >= best_path[1][3]:
                                    best_path = (path_ad, score_w, query_path_ad)
                            else:
                                if (best_path[1][5] != True) and (score_w[4] >= best_path[1][4]):
                                    best_path = (path_ad, score_w, query_path_ad)
                                    #print best_path
                    else:
                        v_stack.append(next)
                        p_stack.append(path + [next])
            else:
                p_found = 0
                for next in search_graph[vertex]:
                    if next == "END":
                        if mode == 1:
                            paths.append(path)
                        else:
                            path_ad = path + a_s_f.keys()
                            query_path_ad = query_path + a_q_f.keys()
                            logging.debug("path_ad: "+str(path_ad))
                            logging.debug("query_path_ad: "+str(query_path_ad))
                            
                            score_w = sf_entire_calc_score(path_ad, query_path_ad)
                            if score_w[5] == True:
                                if score_w[3] >= best_path[1][3]:
                                    best_path = (path_ad, score_w, query_path_ad)
                            else:
                                if (best_path[1][5] != True) and (score_w[4] >= best_path[1][4]):
                                    best_path = (path_ad, score_w, query_path_ad)
                            p_found = 1
                    elif mode == 1:
                        if query_features[next][0] == priority:
                            v_stack.append(next)
                            p_stack.append(path + [next])
                            p_found = 1
                    elif mode == 0:
                        if search_features[next][0] == priority:
                            v_stack.append(next)
                            p_stack.append(path + [next])
                            p_found = 1
                                                        
                if p_found == 0:
                    if mode == 1:
                        best_priority_bridger = ("NONE",(0.0, 0.0, 0.0, 0.0, 0.0, False))                        
                        for next in search_graph[vertex]:
                            ###baustelle: fixed: path elongation without additives
                            score = sf_calc_score(path+[next], mode, protein)
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
                            score_w = sf_entire_calc_score(path_ad+[next],query_path_ad)
                            if score_w[5] == True:
                                if score_w[3] >= best_partial_path[1][3]:
                                    best_partial_path = (next, score_w)
                            else:
                                if (best_partial_path[1][5] != True) and (score_w[4] >= best_partial_path[1][4]):
                                    best_partial_path = (next, score_w)
                                    
                        v_stack.append(best_partial_path[0])
                        p_stack.append(path + [best_partial_path[0]])
        if mode == 1:
            logging.debug("returning: "+str(paths))
            logging.debug("cardinality: "+str(len(paths)))
            return paths
        elif mode == 0:
            logging.debug(best_path)
            return best_path

# walks the graph in priority-mode to reduce number of paths
def pb_graphtraversal_priority(search_graph, priority):
	global search_features
        
        logging.info("pb_graphtraversal_priority")
        logging.debug(search_graph)
        
        v_stack = ["START"]
	p_stack = []
	best_path = ([], (0.0, 0.0, 0.0, 0.0))
	while v_stack:
		vertex = v_stack.pop()
		if len(p_stack) == 0:
			path = []
		else:
			path = p_stack.pop()
                if len(search_graph[vertex]) < 2: ## unsafe threshold ##
			for next in search_graph[vertex]:
				if next == "END":
                                        score = sf_calc_score(path, 0, "none")
					if score[3] >= best_path[1][3]:
						best_path = (path, score)
				else:
					v_stack.append(next)
					p_stack.append(path + [next])
		else:
			p_found = 0
			for next in search_graph[vertex]:
				if next == "END":
					score = sf_calc_score(path, 0, "none")
					if score[3] >= best_path[1][3]:
						best_path = (path, score)
					p_found = 1
				elif search_features[next][0] == priority:
                                        v_stack.append(next)
					p_stack.append(path + [next])
					p_found = 1
			if p_found == 0:
                                #some kind of greedy strategy: if feature type (p priority) not found
				best_partial_path = ("NONE", (0.0, 0.0, 0.0, 0.0))
				for next in search_graph[vertex]:
					score = sf_calc_score(path+[next], 0, "none")
					if score[3] >= best_partial_path[1][3]:
                                                best_partial_path = (next, score)
                                v_stack.append(best_partial_path[0])
                                p_stack.append(path + [best_partial_path[0]])
	return best_path

	

########## Scoring Functions ########## <sf>
# calculate the fas score sf_calc_score is the main function


# main scoring function
# calculates the final score from MS, CS and PS
def sf_calc_score(path, mode, protein):
    global score_weights
    global weight_const
    global weights

    if weight_const == 1:
        adjusted_weights = w_weight_const_rescale(path)
        tmp_weight = {}
        for i in adjusted_weights:
            tmp_weight[i] = weights[i]
            weights[i] = adjusted_weights[i]
    if score_weights[1] != 0.0: #baustelle: correct?
        score_CS = round(sf_CS_score(path, mode), 10)
    else:
        score_CS = 0.0
    tmp = sf_MS_score(path, mode, protein)
    score_MS = round(tmp[0], 10)
    score_PS = round(sf_PS_score(path, tmp[1], tmp[2], mode, protein), 10)
    final_score = (score_MS * score_weights[0]) + (score_CS * score_weights[1]) + (score_PS * score_weights[2])
    if weight_const == 1:
        for i in adjusted_weights:
            weights[i] = tmp_weight[i]
    return (score_MS, score_PS, score_CS, final_score)

# main scoring function - entire mode
# calculates the final score from MS, CS and PS
def sf_entire_calc_score(path, query_path):
    global score_weights
    global weight_const
    global weights

    if weight_const == 1:
        adjusted_weights = w_weight_const_rescale(path)
        tmp_weight = {}
        for i in adjusted_weights:
            tmp_weight[i] = weights[i]
            weights[i] = adjusted_weights[i]
    if score_weights[1] != 0.0: #baustelle
        score_CS = round(sf_entire_cs_score(path, query_path), 10)
    else:
        score_CS = 0.0
    tmp = sf_entire_ms_score(path, query_path)
    score_MS = round(tmp[0], 10)
    path_weight = tmp[3]
    common_feature = tmp[4]
    score_PS = round(sf_entire_ps_score(path, tmp[1], tmp[2], query_path), 10)
    final_score = (score_MS * score_weights[0]) + (score_CS * score_weights[1]) + (score_PS * score_weights[2])
    if weight_const == 1:
        for i in adjusted_weights:
            weights[i] = tmp_weight[i]
    return (score_MS, score_PS, score_CS, final_score, path_weight, common_feature)

# calculates Clan Score 
def sf_CS_score(path, mode):
    global clan_dict
    global query_clans
    global search_features

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
def sf_MS_score(path, mode, protein):
	global search_features
	global weights
        global single_protein
	global query_protein
	global MS_uni
        global proteome

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
                if one_against_all == 1:
                    logging.debug("single_protein[1]: "+ str(single_protein[1]))
                    if feature in single_protein[1]:
                        logging.debug("feature count: "+str(feature)+"=> "+str(len(single_protein[1][feature])-2))
                        s_length = len(single_protein[1][feature])-2
                        p_score =  float(search_domains[feature] * s_length) / float(max(search_domains[feature], s_length) * max(search_domains[feature], s_length))
                        scores.append((feature,p_score))
                    else:
                        scores.append((feature, 0.0))
                elif one_against_all == 0:
                    if feature in proteome[protein]:
                        s_length = len(proteome[protein][feature])-2
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
def sf_entire_ms_score(path, query_path):
	global search_features
        global a_s_f
        global query_features
        global a_q_f
	global weights
	global MS_uni

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
def sf_PS_score(path, search_domains, scale, mode, protein):
    global search_features
    global query_features
    global query_protein
    global single_protein
    global weights
    global MS_uni
    global proteome

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
            logging.debug("feature: "+str(feature))
            logging.debug("single_protein: "+str(single_protein))
            #position = ((float(instance[0]) + float(instance[1])) / 2.0) / float(protein_lengths[protein_id])
            if one_against_all == 1:
                if feature[0] in single_protein[1]:
                    best_match = 0.0
                    scores[feature[0]] = 0.0
                    for instance in single_protein[1][feature[0]][2:]:
                        pos = (float(instance[1]) + float(instance[2])) / 2.0 / float(protein_lengths["single_"+str(single_protein[0])])
                        match = 1.0 - float(abs(feature[1]) - pos)
                        logging.debug(str(float(instance[1]))+" + "+str(float(instance[2]))+" / 2.0 / "+str(float(protein_lengths["single_"+str(single_protein[0])]))+" = "+str(pos))
                        if best_match < match:
                            best_match = match
                    scores[feature[0]] += best_match
            elif one_against_all == 0:
                if feature[0] in proteome[protein]:
                    best_match = 0.0
                    scores[feature[0]] = 0.0
                    for instance in proteome[protein][feature[0]][2:]:
                        pos = (float(instance[1]) + float(instance[2])) / 2.0 / float(protein_lengths["set_"+str(protein)])
                        match = 1.0 - float(abs(feature[1]) - pos)
                        logging.debug(str(float(instance[1]))+" + "+str(float(instance[2]))+" / 2.0 / "+str(float(protein_lengths["set_"+str(protein)]))+" = "+str(pos))
                        if best_match < match:
                            best_match = match
                    scores[feature[0]] += best_match                    
        
    for f_score in scores:
        if MS_uni == 0 and mode == 0:
            final_score += scores[f_score] * scale * weights[f_score]
        else:
            final_score += scores[f_score] * scale
    return final_score

# calculates Positional Score - entire mode
def sf_entire_ps_score(path, search_domains, scale, query_path):
    global search_features
    global a_s_f
    global query_features
    global a_q_f
    global weights
    global MS_uni

    final_score = 0.0
    scores = {}
    local_query_protein = {}

    # get current features from query path
    for i in query_path:
        if i in query_features:
            feature = query_features[i]
        else:
            feature = a_q_f[i]

        if feature[0] in local_query_protein:
            local_query_protein[feature[0]].append(feature[1])
        else:
            local_query_protein[feature[0]] = [feature[1]]

    logging.debug(str(local_query_protein)+" for query_path "+str(query_path))
            
            
    # compare features in path with features form query path
    for i in path:
        if i in search_features:
            feature = search_features[i]
        else:
            feature = a_s_f[i]
        if feature[0] in local_query_protein:
            best_match = 0.0
            scores[feature[0]] = 0.0
            for position in local_query_protein[feature[0]]:
                match = 1.0 - float(abs(feature[1] - position))
                if best_match < match:
                    best_match = match
            scores[feature[0]] += best_match
    for f_score in scores:
        if MS_uni == 0:
            final_score += scores[f_score] * scale * weights[f_score]
        else:
            final_score += scores[f_score] * scale
    return final_score

########## Weighting ########## <w>
# weighting functions

# counts all domains in single_protein and proteome
def w_count():
	global single_protein
	global domain_count
	global proteome
	global one_against_all
        
        logging.debug("w_count: counting domains in single protein and protein set.")

	for feature in single_protein[1]:
            if feature not in domain_count:
		domain_count[feature] = 1
	for i in proteome:
            protein = proteome[i]
            for feature in protein:
                if feature not in domain_count and one_against_all == 0:
                    domain_count[feature] = 1
        
        logging.debug("domain counts: "+str(domain_count))



# counts all domains in a reference proteome if one is given
def w_count_ref():
	global domain_count
	global proteome
        logging.debug("w_count_ref: counting domains in reference gene set.")

	for i in proteome:
		protein = proteome[i]
		for feature in protein:
                    # counting instances (substract 2 because first and second entry contain assess(bool) and feat_eval(float))
                    count = len(protein[feature]) - 2
                    if feature in domain_count:
                        domain_count[feature] +=  count
                    else:
                        domain_count[feature] = count
        logging.debug("domains counts: "+str(domain_count))



# calculates weights 
def w_weighting(protein, mode):
	global domain_count
	global weights
	global proteome
	global single_protein

	scaling_factor = 0.0
	sum_of_features = 0.0
	features = []
	if protein == "none":
            for feature in single_protein[1]:
                features.append(feature)
	else:
            for feature in proteome[protein]:
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

# calculates weights 
def w_weighting_counterpart(protein, mode):
	global domain_count
	global weights_counter
	global proteome
	global single_protein

	scaling_factor = 0.0
	sum_of_features = 0.0
	features = []
	if protein == "none":
            for feature in single_protein[1]:
                features.append(feature)
	else:
            for feature in proteome[protein]:
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
            weights_counter[feature] = round(float(sum_of_features) / float(domain_count[feature]), 8)
            scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 8)
	for feature in features:
            weights_counter[feature] = round(float(weights_counter[feature]) / float(scaling_factor), 8)

# calculates weights if constraints are applied            
def w_weighting_constraints(protein, mode):
        global domain_count
        global weights
        global proteome
        global single_protein
        global constraints
        
        tools = {}
        tools["pfam"] = []
        tools["smart"] = []
        tools["cast"] = []
        tools["coils"] = []
        tools["signalp"] = []
        tools["seg"] = []
        tools["tmhmm"] = []
        features = []
        single_constraints = []
        filled = 0.0
	if protein == "none":
            for feature in single_protein[1]:
                features.append(feature)
	else:
            for feature in proteome[protein]:
                features.append(feature)
	for feature in features:
            if feature in constraints:
                filled += constraints[feature]
                weights[feature] = constraints[feature]
                single_constraints.append(feature)
            elif feature[0:4] == "pfam" and "pfam" in constraints:
                tools["pfam"].append(feature)
            elif feature[0:5] == "smart" and "smart" in constraints:
                tools["smart"].append(feature)
            elif feature[0:4] == "cast" and "cast" in constraints:
                tools["cast"].append(feature)
            elif feature[0:5] == "coils" and "coils" in constraints:
                tools["coils"].append(feature)
            elif feature[0:7] == "signalp" and "signalp" in constraints:
                tools["signalp"].append(feature)
            elif feature[0:3] == "seg" and "seg" in constraints:
                tools["seg"].append(feature)
            elif feature[0:5] == "tmhmm" and "tmhmm" in constraints:
                tools["tmhmm"].append(feature)
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

def w_weighting_constraints_counter(protein, mode):
        global domain_count
        global weights_counter
        global proteome
        global single_protein
        global constraints
        
        tools = {}
        tools["pfam"] = []
        tools["smart"] = []
        tools["cast"] = []
        tools["coils"] = []
        tools["signalp"] = []
        tools["seg"] = []
        tools["tmhmm"] = []
        features = []
        single_constraints = []
        filled = 0.0
	if protein == "none":
            for feature in single_protein[1]:
                features.append(feature)
	else:
            for feature in proteome[protein]:
                features.append(feature)
	for feature in features:
            if feature in constraints:
                filled += constraints[feature]
                weights_counter[feature] = constraints[feature]
                single_constraints.append(feature)
            elif feature[0:4] == "pfam" and "pfam" in constraints:
                tools["pfam"].append(feature)
            elif feature[0:5] == "smart" and "smart" in constraints:
                tools["smart"].append(feature)
            elif feature[0:4] == "cast" and "cast" in constraints:
                tools["cast"].append(feature)
            elif feature[0:5] == "coils" and "coils" in constraints:
                tools["coils"].append(feature)
            elif feature[0:7] == "signalp" and "signalp" in constraints:
                tools["signalp"].append(feature)
            elif feature[0:3] == "seg" and "seg" in constraints:
                tools["seg"].append(feature)
            elif feature[0:5] == "tmhmm" and "tmhmm" in constraints:
                tools["tmhmm"].append(feature)
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
                    weights_counter[feature] = round(float(sum_of_features) / float(domain_count[feature]), 8)
                    scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 8)  
                for feature in tools[tool]:
                    weights_counter[feature] = round(float(weights_counter[feature]) / float(scaling_factor) * constraints[tool], 8)
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
            weights_counter[feature] = round(float(sum_of_features) / float(domain_count[feature]), 8)
            scaling_factor += round(float(sum_of_features) / float(domain_count[feature]), 8)
	for feature in features:
            weights_counter[feature] = round(float(weights_counter[feature]) / float(scaling_factor) * (1.0-filled), 8)
              
                        
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
            
def w_weight_const_rescale(path):
        global weights
        global search_features
        global constraints

        pfam = []
        smart = []
        tmp = 0.0
        rescaled_weights = {}
        for feature in path:
            if feature in search_features:
                if search_features[feature][0] not in constraints:
                    if search_features[feature][0][0:4] == "pfam" and search_features[feature][0] not in pfam:
                        pfam.append(search_features[feature][0])
                    elif search_features[feature][0][0:5] == "smart" and search_features[feature][0] not in smart:
                        smart.append(search_features[feature][0])
        if "pfam" in constraints:
            for feature in pfam:
                tmp += weights[feature]
            if tmp < constraints["pfam"] and tmp > 0.0:
                scale = constraints["pfam"] / tmp
                for feature in pfam:
                    rescaled_weights[feature] = weights[feature] * scale 
        tmp = 0.0
        if "smart" in constraints:
            for feature in smart:
                tmp += weights[feature]
            if tmp < constraints["smart"] and tmp > 0.0:
                scale = constraints["smart"] / tmp
                for feature in pfam:
                    rescaled_weights[feature] = weights[feature] * scale 
                    
        return rescaled_weights

def w_weight_const_rescale_final(path):
        global weights
        global constraints

        pfam = []
        smart = []
        tmp = 0.0
        rescaled_weights = {}
        for feature in path:
            if feature not in constraints:
                if feature[0:4] == "pfam" and feature not in pfam:
                    pfam.append(feature)
                elif feature[0:5] == "smart" and feature[0] not in smart:
                    smart.append(feature)
        if "pfam" in constraints:
            for feature in pfam:
                tmp += weights[feature]
            if tmp < constraints["pfam"] and tmp > 0.0:
                scale = constraints["pfam"] / tmp
                for feature in pfam:
                    rescaled_weights[feature] = weights[feature] * scale 
        tmp = 0.0
        if "smart" in constraints:
            for feature in smart:
                tmp += weights[feature]
            if tmp < constraints["smart"] and tmp > 0.0:
                scale = constraints["smart"] / tmp
                for feature in pfam:
                    rescaled_weights[feature] = weights[feature] * scale 
                    
        return rescaled_weights


########## Input ########## <OK>

# reads constraints file
def constraints_in(path):
        global constraints

        cfile = open(path, "r+")
        lines = cfile.readlines()
        for i in range(7):
            split = (lines[i].rstrip("\n")).split(" ")
            if split[1] != "N":
                constraints[split[0]] = float(split[1])
        if len(lines) > 7:
            for i in range(8, len(lines)):
                split = (lines[i].rstrip("\n")).split(" ")
                constraints[split[0]] = float(split[1])

                                

# reads the xml-files as input
def xmlreader(path, single, tool, assess):
	global single_protein
	global proteome
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
			if (single == 0 or single == 2) and not (pID in proteome):
				proteome[pID]= {}
			elif single == 1 and single_protein == ():
				single_protein = (pID, {}) 

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
                                        proteome[pID][ftype] = []
                                        proteome[pID][ftype].append(assess)
                                        proteome[pID][ftype].append(feat_eval)
                                    else:
                                        single_protein[1][ftype] = [] 
                                        single_protein[1][ftype].append(assess)
                                        single_protein[1][ftype].append(feat_eval)
                                    
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
                                                    proteome[pID][ftype].append((inst_eval, start, end))
                                                    inst_count += 1
                                                    #print proteome[pID][ftype]

                                                else:
                                                    single_protein[1][ftype].append((inst_eval, start, end))
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
                                                    proteome[pID][ftype].append((inst_eval, start, end)) 
                                                else:
                                                    single_protein[1][ftype].append((inst_eval, start, end))

                                            else:
                                                if i == 0:
                                                    start = int(instance.attrib["start"])
                                                    i = 1
                                                else:
                                                    end = int(instance.attrib["end"])
                                                    if single == 0 or single == 2:
                                                            proteome[pID][ftype].append((inst_eval, start, end)) 
                                                    else:
                                                            single_protein[1][ftype].append((inst_eval, start, end))
                                                    i = 0
                                    # any instance appended?
                                    if (inst_count < 1):
                                        #delete feature type
                                        logging.info("Rejecting feature type " + str(ftype) + " due to rejection of all instances. Check for thresholds and E-values (instance based)")
                                        if single == 0 or single == 2:
                                            proteome[pID].pop(ftype)
                                        else:
                                            single_protein[1].pop(ftype)
                                    if ftype not in clan_dict:
                                        clan_dict[ftype] = fclan
        logging.debug("Proteome: "+str(proteome))

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

        
fc_main()

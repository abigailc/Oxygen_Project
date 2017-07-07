#!/usr/bin/python

#abigailc@Actaeon Jan 4 2017; edit july 6 2017


#this is going to be an outside module that can take input of (gene tree, species tree) and run ranger dtl, return some result.

#THEORETICALLY if a progress file already exists, we can now go straight to doing ranger, subsampling, and analysis.
##might not work because information that we are not tracking includes 1) loss sequences 2) rooting sequences.
##this should be easy enough to re-work out, but maybe should have been tracking the whole time...
##a final thing: the difference between short_clades (singleton) vs (exterior to anything else) (if that is even tracked?????) vs ERRORED.

# def Get_All_Unseen_Tips(???):
     #this should return all tips that were not inlcuded in any of the trees examined (big/small)
     #will include: external loners
######### PERSONAL_SETTINGS #########
ssh_inst = "ssh -l abigailc -i ~/.ssh/id_rsa eofe4.mit.edu"
clus_head = "abigailc@eofe4.mit.edu:/home/abigailc/"

import os
import time
import re

def Master_Ranger_SS(list_of_big_clades, projectname, list_of_short_clades, original_fasta):
    #list of big clades is the subtree object for each
    print("Beginning Master Ranger SS in Run_Ranger.py")
    base_dir = os.getcwd()
    os.chdir(projectname)

    os.system("mkdir Ranger_DTL")
 
    os.chdir("Ranger_DTL")
    os.system("mkdir ranger_input; mkdir ranger_output; mkdir biparts")
    ###############BOOTSTRAP FUNCTION###################
    print("num_big_clades: "+str(len(list_of_big_clades)))
    good, error = Check_Boots_Exist(list_of_big_clades, base_dir)
    print("gene_tree_bootstraps_exist: "+str(len(good)))
    good, error = Check_Species_Exist(good, error, base_dir)
    print("species_besttrees_exist: "+str(len(good)))
    ranger_names, i = Make_Ranger_In_For_Each_Bootstrap(good)
    print("ranger_input_files:")
    print(ranger_names)
    #check for already existing bipartitions files.. if they exist, skip.
    exist_already = 0
    exist_list = []
    for item in list_of_big_clades:
        bipartition_table = item.prefix+"biparts_table.txt"
        exists = os.path.isfile("biparts/"+bipartition_table)
        if exists is True:
            exist_already+=1
            ranger_names.remove(item.prefix+"RangerIn.txt")
            exist_list.append(item.prefix+"RangerOut.txt")
    print(" Bipart files already exist for: "+str(exist_already))

    if ranger_names == []:
        ranger_outputs = []
        print("all bipartitions tables already existed")
    else:
    #these names should be JUST THE BASE and not include numbers in front of them.
    #    ranger_names = Make_Ranger_Input_Files(list_of_big_clades)
        ranger_outputs = RunRangerOnInputs_boots(ranger_names, projectname, i)
    #ranger_outputs includes... just the base, or also the number??
#    print(ranger_outputs)
    for item in exist_list:
        ranger_outputs.append(item)
    print("rangers should exist")
    #just attaches things to the subtree object
    Add_Ranger_To_Subtree_Boots(ranger_outputs, list_of_big_clades)
    ######################PARSING FUNCTIONS#############
    for item in list_of_big_clades:
        rangerfile = item.rangerout_file
        rangerfiles = item.rangerout_list
        bipartition_table = item.prefix+"biparts_table.txt"
        print(item.prefix)
        if os.path.isfile(bipartition_table) is True:
            bipart_file = bipartition_table
            transfer_recips_list = []
        else:
            bipart_file, empty_list = Parse_Ranger_Clade_BOTHWAYS(rangerfiles, bipartition_table)
#            bipart_file, transfer_recips_list = Parse_Ranger_Clade_M_not_N(rangerfiles, bipartition_table)
        item.bipartitions_file = bipart_file
#        item.transfer_recipient_tips = transfer_recips_list
    
    print("finished making bipartitions tables.")




    #clean up the utter shitshow that is left on the cluster in clus_head/Ranger
    print("moving files on the cluster... don't want to leave a mess! ;)")
    os.system(ssh_inst+" \'cd Ranger; mkdir "+projectname+"; mv *"+projectname+"* "+projectname+"/\'")
    print("reorganizing files locally...")

    os.system("mv *RangerIn.txt ranger_input/")

    os.system("mv *RangerOut.txt ranger_output/")
    os.system("mv *biparts_table.txt biparts/")

    #go back to MakeSpeciesTrees
    os.chdir("..")
    os.chdir("..")

    #next step : compile and send a list of "tips to ignore (transfer recipients)" to Choose Seqs Large ST.
    #this function was removed because EVERY TIP was part of a deeper transfer somewhere. could maybe do something to avoid really shallow transferS??

    ######################SUBSAMPLING FUNCTIONS########

#    Add_Deepest_Split_To_Subtree(list_of_big_clades)
#    list_of_big_clades, list_of_short_clades = Check_if_Big_Clade_Is_Sub(list_of_big_clades, list_of_short_clades)
    
 #   seqs_to_keep = Choose_Seqs_Large_ST(list_of_big_clades)
    return "None yet"


def Master_Ranger_Single(subtree_ob, projectname):
    #i still want this to do it with bootsraps.
    os.system("mkdir "+projectname+"/SS_Ranger_DTL")
    #print(subtree_ob.besttree_str)
    print("trying to move to a better place.")
    os.system("cp "+subtree_ob.besttree_str+" /"+projectname+"/SS_Ranger_DTL")
    os.chdir(projectname)
    os.chdir("SS_Ranger_DTL")
    ###############see if species tree and gene tree exist###################
    list_of_big_clades = [subtree_ob]
    ###############BOOTSTRAP FUNCTION###################
    good, error = Check_Boots_Exist(list_of_big_clades)
    print("gene_tree_bootstraps_exist: "+str(len(good)))
    good, error = Check_Species_Exist(good, error)
    print("species_besttrees_exist: "+str(len(good)))
    #make the ranger input files
    #this is messing up hard. printing "besttree" instead of the tree itself!
    #fix a mistake


    a = os.path.isfile(subtree_ob.besttree_str)
    if a is False:
        c = os.path.isfile(subtree_ob.besttree_str)
        if c is True:
            subtree_ob.besttree_str = subtree_ob.besttree_str
        else:
            b = raw_input("give me the full path to and name of the species besttree because this shit isn't working")
            subtree_ob.besttree_str = b
    with open(subtree_ob.besttree_str) as old:
        for line in old:
            if line != "":
                besttree_str_no_really = line
                break
    subtree_ob.besttree_str = besttree_str_no_really
    ranger_names, i = Make_Ranger_In_For_Each_Bootstrap(good)
    
    print("ranger_input_files:")
    #print(ranger_names)
    #these names should be JUST THE BASE and not include numbers in front of them.
    #    ranger_names = Make_Ranger_Input_Files(list_of_big_clades)
    #make ranger outputs
    ranger_outputs = RunRangerOnInputs_boots(ranger_names, projectname, i)
    #print(ranger_outputs)
    print("rangers should exist")
    #just attaches things to the subtree object
    Add_Ranger_To_Subtree_Boots(ranger_outputs, list_of_big_clades)
    ######################PARSING FUNCTIONS#############
    for item in list_of_big_clades:
        rangerfile = item.rangerout_file
        rangerfiles = item.rangerout_list
        bipartition_table = item.prefix+"biparts_table.txt"
        #print(item.prefix)
        bipart_file, transfer_recips_list = Parse_Ranger_Clade_BOTHWAYS(rangerfiles, bipartition_table)
        item.bipartitions_file = bipart_file
        item.transfer_recipient_tips = transfer_recips_list
    #parsing

    
    print("finished making bipartitions tables for single species and gene tree bootstraps.")
    print(bipart_file)
    print("you should run add numbers to nodes to visualize it?")
    

#obsolete -- USED FOR EASYMODE
def RunRangerOnInputs(ranger_names, projectname):
    list_ranger_outs = []
    for item in ranger_names:
        output_name = item[:-6]+"Out.txt"
        list_ranger_outs.append(output_name)
    move_ran_cluster(ranger_names)
    run_ranger_on_cluster(ranger_names, list_ranger_outs, projectname)
    move_home_ranger_cluster(list_ranger_outs)
        
    #        os.system("rangeru -i "+item+" -o "+output_name)
        #someday, this should be parallelized or run on the cluster.
        #would need to generate a script, do scp there and back verification, and get rangeru installed on cluster.
    return list_ranger_outs





#obsolete
def move_ran_cluster(ranger_names):
    clus_path = "Ranger/"
    for item in ranger_names:
         os.system("scp "+item+" "+clus_head+clus_path)

#obsolete
def run_ranger_on_cluster(ranger_in_names, ranger_out_names, projectname):
    indexname = "Ranger_"+projectname+"_Corr.txt"
    scriptfile = "Ranger_"+projectname+"_Sc.sh"
    a = gen_ranger_script(str(len(ranger_in_names)), indexname, scriptfile)
    b = gen_ranger_corr(ranger_in_names, ranger_out_names, indexname)
    move_ran_cluster([a, b])
    os.system(ssh_inst+" 'cd ~/Ranger/;echo $PWD;sbatch "+a+"'")

#obsolete
def gen_ranger_corr(ran_in, ran_out, indexname):
    with open(indexname, "w") as corr:
        for i in range(len(ran_in)):
            corr.write(str(i)+" "+ran_in[i]+" "+ran_out[i]+"\n")
    return indexname

#obsolete
def gen_ranger_script(n, indexname, scriptfile):
    a =  """#!/bin/bash                                                                   
#SBATCH -p sched_mit_g4nier                                                               
#SBATCH -t 0-10:00:00    
#SBATCH -J Ranger   
#SBATCH -o Ranger.out                                                                     
#SBATCH --array=1-"""+n+"""
. /etc/profile.d/modules.sh
module add engaging/openmpi/1.8.8
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
THE_OUTPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $3}' )
THEROOT=${THE_INPUT_FILE%%.*}
OPT=OptRoot.txt
OPTOUT=$THEROOT$OPT
RAN=$PostOpt_RangerIn.txt
RANIN=$THEROOT$RAN
./OptRoot.linux -i $THE_INPUT_FILE -o $OPTOUT
mpirun python OptRootParser.py -i $OPTOUT -o $RANIN
./Ranger-DTL.linux -i $RANIN -o $THE_OUTPUT_FILE   
exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile

#obsolete
def move_home_ranger_cluster(list_ranger_outs):
    return "done"
    movehome = []
    direct = os.getcwd()
    clus_path = "/Ranger"
    iteration = 0
    for i in list_ranger_outs:
        movehome.append(i)
    finished = False
    while finished is not True:
        #try and move each home.
        for filename in movehome:
            os.system("scp "+clus_head[:-1]+clus_path+"/"+filename+" "+direct)
        finished = "yes"
        for item in list_ranger_outs:
            #see if everything got moved home.
            exists = os.path.isfile(item)
            if exists is True:
                if item in movehome:
                    movehome.remove(item)
            else:
                finished = False
                print(item+" not found")
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            #wait five minutes and then try again.
            
            iteration += 1
            print("CURRENTLY SKIPPING FOR SPEED")
            return("done")
            if iteration < 10:
                print("checking.... some ranger outputs do not exist yet. sleeping for 5 minutes.")
                time.sleep(300)
            elif iteration > 10:
                a = raw_input("type resume to resume, skip to skip to next step")
                if a == "skip":
                    return "done"
            else:
                print("checking.... some ranger outputs do not exist yet. sleeping for 15 minutes.")
                time.sleep(900)
#obsolete
def Make_Ranger_Input_Files(list_of_big_clades):
    list_of_bc_ranger_names = []
    for subtree in list_of_big_clades:
        new_name = subtree.ret_prefix()+"RangerIn.txt"
        list_of_bc_ranger_names.append(new_name)
        with open (new_name,"w") as new:
            #FOR NOW COMPLETELY REMOVING _ BECAUSE THEY RUIN EVERYTHING. SO DO # and ~
            genet = subtree.gene_tree_species_tips
            genet = genet.replace("_", "XX0XX")
            spect = subtree.species_tree
            spect = spect.replace("_", "XX0XX")
            new.write("[&R]"+spect+"\n")
            new.write("[&U]"+genet)
           
    return list_of_bc_ranger_names

#########################SUBSAMPLING FUNCTIONS ######################

def Add_Deepest_Split_To_Subtree(list_of_big_clades):
    for item in list_of_big_clades:
        first, last = DeepestSplit(item.gene_tree_species_tips)
        item.first_list = first
        item.last_list = last

def Add_Ranger_To_Subtree_Boots(ranger_outputs, list_of_big_clades):
    #add each rangeroutput name to subtree for each.
    if len(ranger_outputs) == len(list_of_big_clades):
        pass
    else:
        print("error, rangerouts do not match bigclades in length")
    i = 0
    for largest in list_of_big_clades:
        list_of_rangers = []
        largest.rangerout_file = ranger_outputs[i]
        for num in range(100):
            num += 1
            #this deals with if you couldn't run all 100 bootstraps for any reason. 
            a=os.path.isfile(str(num)+ranger_outputs[i])
            if a is False:
                b=os.path.isfile("ranger_output/"+str(num)+ranger_outputs[i])
                if b is False:
                #warning message.
                    print("total of :"+str(num-1)+" rangerfiles found for "+largest.prefix)
                    break
                if b is True:
                    list_of_rangers.append("ranger_output/"+str(num)+ranger_outputs[i])
            if a is True:
                list_of_rangers.append(str(num)+ranger_outputs[i])
        largest.rangerout_list = list_of_rangers
        i += 1

def Add_Ranger_To_Subtree(ranger_outputs, list_of_big_clades):
    #add each rangeroutput name to subtree for each.
    #open and parse each ranger out - result should be (???? for transfers) AND a list of recips to ignore....
    #recipe ignore list needs to be saved to the subtree.
    if len(ranger_outputs) == len(list_of_big_clades):
        pass
    else:
        print("error, rangerouts do not match bigclades in length")
    i = 0
    for largest in list_of_big_clades:
        largest.rangerout_file = ranger_outputs[i]
        i += 1
        print("You need to code something that parses the output of rangerDTL and identifies all transfer recipient tips, so they can be excluded by tip-selector.")
        largest.trans_recips = ParseRangerOut(largest.rangerout_file)


#####BOOTSTRAP SPECIFIC THINGS######

#####BOOTSTRAP SPECIFIC THINGS######

#####BOOTSTRAP SPECIFIC THINGS######


def Make_Ranger_In_For_Each_Bootstrap(list_of_big_clades):
    ranger_in_base_names = []
    
    if len(list_of_big_clades) == 0:
        print("Something has gone wrong, or there isn't enough data. No clades remain to be run through rangerDTL!!")
        raise SystemExit
    print("making ranger input files for each bootstrap from a total number of clades: "+str(len(list_of_big_clades)))
    for subtree in list_of_big_clades:
        subtree.bootstrap_trees = Separate_Bootstraps(subtree.boot_gene_file)
        subtree.boots_gene_with_species_names = []
        subtree.ranger_in_list = []
        i = 0
        ranger_in_base_names.append(subtree.ret_prefix()+"RangerIn.txt")
        for boot in subtree.bootstrap_trees:
            i+=1
            new_name = str(i)+subtree.ret_prefix()+"RangerIn.txt"
            subtree.ranger_in_list.append(new_name)
            with open (new_name,"w") as new:
                genet = Make_Species_Tips(boot, subtree)
                spect = subtree.besttree_str
                spect = spect.replace("_", "XX0XX")  
                new.write("[&R]"+spect+"\n")
                new.write("[&U]"+genet)
    return ranger_in_base_names, i
    #base name will be something like: PREFIX+RangerIn.txt


#this should turn nodes into a dictionary for a single tree given in string form.
def Newick_Tree_Nodes_To_Dict(genetree_string):
    #genetree_string = (((A,B)n1,(C,D)n2)n3),E)n4
    node_to_tips_dict = {}
    finished = False
    overall_tips_list = []
    while finished is False:   
        node_name = re.sub ("(.*)(\()([^\(|\)]*)(\))([a-z][0-9]*)([^\)]*)(.*)", "\\5", genetree_string)
        node_tips_string = re.sub ("(.*)(\()([^\(|\)]*)(\))([a-z][0-9]*)([^\)]*)(.*)", "\\3", genetree_string)
        #the tips will be seperated by commas
        genetree_string = re.sub ("(.*)(\()([^\(|\)]*)(\))([a-z][0-9]*)([^\)]*)(.*)", "\\1\\3\\6\\7", genetree_string)
        node_tips_string = node_tips_string.strip()
        node_name = node_name.strip()
        tips_list = node_tips_string.split(",")
        #sort so that m1 = a,b is the same as m1 = b,a
        #m1 = [a,b]
        tips_list = sorted(tips_list)
        for item in tips_list:
            if item not in overall_tips_list:
                overall_tips_list.append(item)
        node_to_tips_dict[node_name] = tips_list
        it = 0
        for item in genetree_string:
            if item == "(":
                it += 1
        if it == 0:
            finished = True
    #print(node_to_tips_dict)
    #I think I also want to add add the single tips to the dictionary with node_name = tip and tips_list = tip.
    for tip in overall_tips_list:
        node_to_tips_dict[tip] = [tip]

    return node_to_tips_dict
    #this dict will be used to determine which tips received the gene via a transfer even if the transfer is deeper in the tree.

#obsolete
def Parse_Single_Ranger_Output(rangerfile):
    transfer_tups = []
    transfer_list = []
    loss_tups = []
    dup_tups = []
    gene_tree_dict = {}
    

    #each tuple should be of form "FROM HERE" "TO HERE"
    # eg transfer ("LCA(cat,dog)"),("Methanogen")
    #transfer list should contain TIPS that received thing in transfer. tips only.
    #implement this stuff
    prep = False
    with open (rangerfile) as ranger:
        for line in ranger:
            if "Species Tree:" in line:
                prep = True
            elif "minimum reconciliation cost" in line:
                continue
            elif prep is True:
                prep = False
                gene_tree_dict = Newick_Tree_Nodes_To_Dict(line)
                #print(gene_tree_dict)
            elif "Speciation" in line:
                continue
    #             do a thing
            elif "Transfer" in line:
                #add the transfer to the tupo list     
                donor = re.sub("(.*)(Mapping --> )([^,]*)(.*)" , "\\3",line)
                donor = donor.strip()
                recipient = re.sub("(.*)(Recipient --> )([^,]*)(.*)" , "\\3",line)
                recipient = recipient.strip()
                recipient = recipient
                transfer_tups.append((donor, recipient))
                #print(type(gene_tree_dict))
                #print(gene_tree_dict)
                try:
                    #this works if the recipient is in the dictionary of gene trees (m1, m2, m3, m4)
                    list_of_tips_that_recieved_in_transfer = gene_tree_dict[recipient]
                except KeyError:
                    list_of_tips_that_recieved_in_transfer = [recipient]
                    print("Sanity check... is "+recipient+" a single tip?")
                #add all tips that received the gene as a result of a transfer anywhere upsteam to a list.
                for item in list_of_tips_that_recieved_in_transfer:
                    if item in transfer_list:
                        pass
                    else:
                        transfer_list.append(item)
    #             do other thing
            elif "Loss" in line:
                donor = re.sub("(.*)(Mapping --> )([^,]*)(.*)" , "\\3",line)
                loss_tups.append(donor)
    #             do the third thing
            elif "Duplication" in line:
                donor = re.sub("(.*)(Mapping --> )([^,]*)(.*)" , "\\3",line)
                dup_tups.append(donor)
    #             do the third thing
            elif "Leaf Node" in line:
                continue
                #re.sub(
                #Leaf Node
            else:
                continue
    #             do a fourth thing
    #             do a fourth thing
    #         else:
    #             print(line)
    #             print("^ was not processed")
   
        
    return transfer_tups, loss_tups, dup_tups, transfer_list, gene_tree_dict

#obsolete    
def Parse_Single_Ranger_Output_M_not_N(rangerfile):
    transfer_clades = []
    transfer_list = []
    loss_tups = []
    dup_tups = []
    gene_tree_dict = {}
    speciation_list = []
    new_transfer_list = []
    species_tree_dict = {}
    transfers_mapped_to_species_recips = []

    #each tuple should be of form "FROM HERE" "TO HERE"
    # eg transfer ("LCA(cat,dog)"),("Methanogen")
    #transfer list should contain TIPS that received thing in transfer. tips only.
    #implement this stuff
    prep = False
    sp_prep = False
    with open (rangerfile) as ranger:
        for line in ranger:
            if "Gene Tree:" in line:
                prep = True
            elif "Species Tree:" in line:
                sp_prep = True
            elif "minimum reconciliation cost" in line:
                continue
            elif prep is True:
                prep = False
                gene_tree_dict = Newick_Tree_Nodes_To_Dict(line)
                #rint(gene_tree_dict)
                #rint("gene_d")
            elif sp_prep is True:
                sp_prep = False
                species_tree_dict = Newick_Tree_Nodes_To_Dict(line)
                #print(gene_tree_dict)
            elif "Speciation" in line:
                node = re.sub("([^ ]*)(.*)" , "\\1",line)
                node = node.strip()
                list_of_tips = gene_tree_dict[node]
                str_of_tips = ' '.join(list_of_tips)
                speciation_list.append(str_of_tips)
                #rint(str_of_tips)
                #rint("spec")
                #aise SystemExit
    #             do a thing
    # m8 = LCA[PyrobaculumXX0XXneutrophilum_501318754, PyrobaculumXX0XXoguniense_504114203]: Transfer, Mapping --> n11, Recipient --> n12
            elif "Transfer" in line:
                node = re.sub("([^ ]*)(.*)", "\\1", line)
                node = node.strip()
                r_list = line.split(" ")
                #[m8, LCABLAH, PYROBLAH, Transfer, Mapping, -->, n11, Recipient, -->, n12]
                r_chunk = r_list[-1]
                #n12
                recipient = r_chunk.strip()

                #process the node.
                list_of_tips = gene_tree_dict[node]
                str_of_tips = ' '.join(list_of_tips)
                transfer_clades.append(str_of_tips)
                for item in list_of_tips:
                    if item in transfer_list:
                        #print("uh there are multiple of this transfer?: "+item)
                        pass
                    else:
                        transfer_list.append(item)
                #process the recipint mess
                #make a specient_tree_dict
                list_of_tips = species_tree_dict[recipient]
                #this will give a list of tips a b c that we need to find in one of the tip lists in gene_tree_dict, and extract the key from there.
                poten_list = []
                #this is working for list_of_tips == a single tip

                #right here we have the species tree recipient kept in "list of tips" -- can we just make an output from this?
                str_of_tips = ' '.join(list_of_tips)
                transfers_mapped_to_species_recips.append(str_of_tips)
   
                for potential_clade in gene_tree_dict:
                    #print(list_of_tips)
                    a = Is_X_within_Y(list_of_tips, gene_tree_dict[potential_clade])
                    if a is True:
                        poten_list.append([len(gene_tree_dict[potential_clade]), potential_clade])
                        #print("appending: "+potential_clade)
                        #print(list_of_tips)
                if poten_list == []:
                    print("ERROR")
                #   print("we were looking for a clade containing all of: ")
               #    print(list_of_tips)
              #     print("from the dict:")
             #      print(gene_tree_dict)
            #       print("error with gene transfer in RunRangerline 584")
                    #aise SystemExit
                    continue
                    
                b = sorted(poten_list)
                LCA = b[0]
                LCA = LCA[1]
           #    print(LCA)
          #     print(gene_tree_dict)
                 #process the node.
                list_of_tips = gene_tree_dict[LCA]
                str_of_tips = ' '.join(list_of_tips)
                new_transfer_list.append(str_of_tips)
            






            elif "Loss" in line:
                donor = re.sub("([^ ]*)(.*)" , "\\1",line)
                donor = donor.strip()
                list_of_tips = gene_tree_dict[donor]
                str_of_tips = ' '.join(list_of_tips)
                loss_tups.append(str_of_tips)
    #             do the third thing
            elif "Duplication" in line:
                donor = re.sub("([^ ]*)(.*)" , "\\1",line)
                donor = donor.strip()
                list_of_tips = gene_tree_dict[donor]
                str_of_tips = ' '.join(list_of_tips)
                dup_tups.append(str_of_tips)
    #             do the third thing
            elif "Leaf Node" in line:
                continue
            else:
                #this should ignore the uh initial SPecies Tree": and the species tree itseld and any \n lines
                continue
    print("done with ranger_parseing")
    return speciation_list, transfer_clades, loss_tups, dup_tups, transfer_list, gene_tree_dict, new_transfer_list, transfers_mapped_to_species_recips

def Is_X_within_Y(x_tips, y_tips):
    #for now fix y tips
    y_new = []
    for item in y_tips:
        item_l = item.split("_")
        newitem = item_l[0]
        y_new.append(newitem)
    within = True
    for tip in x_tips:
       #print(tip)
        if tip in y_new:
            pass
        else:
        #   print("not in")
         #  print(y_new)
            within = False
            break
    return within

#obsolete
def Parse_Ranger_Clade(rangerfiles, bipartitions_table):
    clade_transfer_dict = {}
    clade_transfer_tips_dict = {}
    clade_loss_dict = {}
    clade_dup_dict = {}
    transfer_recipients = []
    transfer_recipient_dict = {}
    
    a = bipartitions_table.split("biparts")
    summary_table = a[0]+"summary_table.txt"
    
    bs_num = 1
    for item in rangerfiles:
        transfer_tups, loss_tups, dup_tups, transfer_list, species_tree_dict = Parse_Single_Ranger_Output(item)
        print("checking transfers - RECIPIENTS ONLY")
        for extup in transfer_tups:
            transfer_recipients.append(extup[1])
            #this should add recipient/recip clade to the list. trans_recips should be [n24, n31, bos_taurus, n2] etc
        transfer_recipient_dict = dictionary_check_bipart(transfer_recipient_dict, transfer_recipients)
        
        print("checking transfers - donor + recipients")
        #print(transfer_tups)
        clade_transfer_dict = dictionary_check_bipart(clade_transfer_dict, transfer_tups)
        #print(clade_transfer_dict)
        #this passes dict where the keys are (sender, recipient)
        print(":checking losses")
        #print(loss_tups)
        #removed [1] from all print tups. testing.

        clade_loss_dict = dictionary_check_bipart(clade_loss_dict, loss_tups)
        #print(clade_loss_dict)
        #dict is (loss)
        print("checking duplications")
        #dict is #duplication location
        #print(dup_tups)
        clade_dup_dict = dictionary_check_bipart(clade_dup_dict, dup_tups)
        #print(clade_dup_dict)
        print("checking transfer tip recipients")
        #print(transfer_list)

        #dict is ???
        clade_transfer_tips_dict = dictionary_check_bipart(clade_transfer_tips_dict, transfer_list)
        #print(clade_transfer_tips_dict)

        #each thing in the dict is a list. item 1 = total number of times that transfer has been seen. item 2 is a list of "Y" and "N" where y means the transfer was observewd in that BS.
    #each thing has, uh, [4, ["Y","Y","Y","N","N","Y"] ]
    
    with open(bipartitions_table, "w") as bipart:
        bipart.write("Transfers_By_Recipient_Clade\n")
        for i in transfer_recipient_dict:
            i_string = str(i)
            #i_string = '|'.join(map(str, i))
            str2 = ''.join(transfer_recipient_dict[i][1])
            bipart.write(i_string+"\t"+str(transfer_recipient_dict[i][0])+"\t"+str2+"\n")
        bipart.write("\n\n")
        bipart.write("Transfers\n")
        for i in clade_transfer_dict:

            str1 = ''.join(clade_transfer_dict[i][1])
            #i = str(i)
            i_string = '|'.join(map(str, i))
            #in case i is a tuple
            bipart.write(i_string+"\t"+str(clade_transfer_dict[i][0])+"\t"+str1+"\n")
        bipart.write("Losses\n")
        for i in clade_loss_dict:
            i_string = str(i)
            #i_string = '|'.join(map(str, i))

            str2 = ''.join(clade_loss_dict[i][1])
            bipart.write(i_string+"\t"+str(clade_loss_dict[i][0])+"\t"+str2+"\n")
        bipart.write("Duplications\n")
        for i in clade_dup_dict:
            i_string = str(i)
            #i_string = '|'.join(map(str, i))
            str3 = ''.join(clade_dup_dict[i][1])
            bipart.write(i_string+"\t"+str(clade_dup_dict[i][0])+"\t"+str3+"\n")
        bipart.write("TRANSFERTIPS\n")
        for i in clade_transfer_tips_dict:
            #i = i.strip()
            i_string = str(i)
            #i_string = '|'.join(map(str, i))
            str4 = ''.join(clade_transfer_tips_dict[i][1])
            bipart.write(i_string+"\t"+str(clade_transfer_tips_dict[i][0])+"\t"+str4+"\n")

        ####this should work. pass in a bipartitions_table string.
        transfer_tips = []
        for item in clade_transfer_tips_dict:
            if clade_transfer_tips_dict[item][0] > 80:
                transfer_tips.append(str(item))
        bipart.write("species_tree_dict")
        for node in species_tree_dict:
            bipart.write("\n"+node+"\t"+' '.join(species_tree_dict[node]))
        bipart.write("TIPS THAT RECIEVED TRANSFER AT LEAST 80% of the time: \n")
        for receive_tip in transfer_tips:
            bipart.write(receive_tip+"\n")
        #write the stuff
        #append to the bipart table all the things
        
    print("tips that received the gene via transfer > 80 bs support in the subtree")
    print(transfer_tips)
    summary_dict = {}
    for key in transfer_recipient_dict:
        summary_dict[key] = [transfer_recipient_dict[key][0]]
    for key in clade_dup_dict:
        if key in summary_dict:
            summary_dict[key].append(clade_dup_dict[key][0])
        else:
            summary_dict[key] = [0, clade_dup_dict[key][0]]
    for key in summary_dict:
        try:
            #try returning the second value (dup number) - if not there, add a 0
            a = summary_dict[key][1]
        except:
            summary_dict[key].append(0)
    for key in clade_loss_dict:
        if key in summary_dict:
            summary_dict[key].append(clade_loss_dict[key][0])
        else:
            summary_dict[key] = [0, 0, clade_loss_dict[key][0]]
    for key in summary_dict:
        try:
            #try returning the second value (dup number) - if not there, add a 0
            a = summary_dict[key][2]
        except:
            summary_dict[key].append(0)
    for key in summary_dict:
        sumall = summary_dict[key][0]
        sumall += summary_dict[key][1]
        sumall += summary_dict[key][2]
        speciation = 100 - sumall
        summary_dict[key].append(speciation)
    with open(summary_table, "w") as summ:
        summ.write("node/tip \t transferrecip \t duplication \t loss \t speciation")
        for clade_or_tip in summary_dict:
            summ.write(clade_or_tip+"\t"+str(summary_dict[clade_or_tip][0])+"\t"+str(summary_dict[clade_or_tip][1])+"\t"+str(summary_dict[clade_or_tip][2])+"\t"+str(summary_dict[clade_or_tip][3]))
    print(bipartitions_table)
    print(summary_table)
    return bipartitions_table, transfer_tips

#rangerfiles should be a list of all the bootstraps in a specific clade, bip[artitions_table should be a string
#**********M NOT N

#need to change transfer tracking

#obsolete
def Parse_Ranger_Clade_M_not_N(rangerfiles, bipartitions_table):
    clade_transfer_dict = {}
    clade_transfer_tips_dict = {}
    clade_loss_dict = {}
    clade_dup_dict = {}
    transfer_recipients = []
    transfer_recipient_dict = {}
    clade_speciation_dict = {}
    clade_new_transfer_dict = {}
    transfers_mapped_to_species_recips_dict = {}


    a = bipartitions_table.split("biparts")
    summary_table = a[0]+"summary_table.txt"
    
    bs_num = 1
    for item in rangerfiles:
        speciation_list, transfer_recipients,loss_tups, dup_tups, transfer_list, gene_tree_dict, new_transfer_list, transfers_mapped_to_species_recips = Parse_Single_Ranger_Output_M_not_N(item)
        #print("checking transfers - RECIPIENTS ONLY")
  
            #this should add recipient/recip clade to the list. trans_recips should be [n24, n31, bos_taurus, n2] etc
        transfer_recipient_dict = dictionary_check_bipart(transfer_recipient_dict, transfer_recipients)

        #this one is new for understnading the most fundamental level of ranger parsing
        transfers_mapped_to_species_recips_dict = dictionary_check_bipart(transfers_mapped_to_species_recips_dict, transfers_mapped_to_species_recips)
        #print(":checking losses")
        #print(loss_tups)
        #removed [1] from all print tups. testing.
        clade_loss_dict = dictionary_check_bipart(clade_loss_dict, loss_tups)
        #print(clade_loss_dict)
        #dict is (loss)
        #print("checking duplications")
        #dict is #duplication location
        #print(dup_tups)
        clade_dup_dict = dictionary_check_bipart(clade_dup_dict, dup_tups)
        #print(clade_dup_dict)
        #print("checking transfer tip recipients")
        #print(transfer_list)
        #dict is ???
        clade_transfer_tips_dict = dictionary_check_bipart(clade_transfer_tips_dict, transfer_list)
        #print(clade_transfer_tips_dict)
        #print("checking speciation")
        #print(transfer_list)
        #dict is ???
        clade_speciation_dict = dictionary_check_bipart(clade_speciation_dict, speciation_list)
        #print(clade_transfer_tips_dict)
        clade_new_transfer_dict = dictionary_check_bipart(clade_new_transfer_dict, new_transfer_list)
        #print(clade_transfer_tips_dict)

    with open(bipartitions_table, "w") as bipart:
        #this one is new and untestedz
        bipart.write("Transfers_By_Recipient_Clade\n")
        for i, value in sorted(clade_new_transfer_dict.items(), key=lambda (k, v): v[0], reverse=True):
            i_string = str(i)
            str1 = ''.join(value[1])
            bipart.write(i_string+"\t"+str(clade_new_transfer_dict[i][0])+"\t"+str1+"\n")
        bipart.write("\nTransfers_By_Species_Tree_Recipient\n")
        for i, value in sorted(transfers_mapped_to_species_recips_dict.items(), key=lambda (k, v): v[0], reverse=True):
            i_string = str(i)
            str1 = ''.join(value[1])
            bipart.write(i_string+"\t"+str(transfers_mapped_to_species_recips_dict[i][0])+"\t"+str1+"\n")
        bipart.write("\n")
        bipart.write("Transfers_By_Representative_Node\n")
        for i, value in sorted(transfer_recipient_dict.items(), key=lambda (k, v): v[0], reverse=True):
            i_string = str(i)
            str2 = ''.join(value[1])
            bipart.write(i_string+"\t"+str(transfer_recipient_dict[i][0])+"\t"+str2+"\n")
        bipart.write("\n\n")
        bipart.write("\nLosses\n")
        for i, value in sorted(clade_loss_dict.items(), key=lambda (k, v): v[0], reverse=True):
            i_string = str(i)
            #i_string = '|'.join(map(str, i))
            str2 = ''.join(value[1])
            bipart.write(i_string+"\t"+str(clade_loss_dict[i][0])+"\t"+str2+"\n")
        bipart.write("\nDuplications\n")
        for i, value in sorted(clade_dup_dict.items(), key=lambda (k, v): v[0], reverse=True):
            i_string = str(i)
            #i_string = '|'.join(map(str, i))
            str3 = ''.join(value[1])
            bipart.write(i_string+"\t"+str(clade_dup_dict[i][0])+"\t"+str3+"\n")
        bipart.write("\nSpeciation\n")
        for i, value in sorted(clade_speciation_dict.items(), key=lambda (k, v): v[0], reverse=True):
            i_string = str(i)
            #i_string = '|'.join(map(str, i))
            str3 = ''.join(value[1])
            bipart.write(i_string+"\t"+str(clade_speciation_dict[i][0])+"\t"+str3+"\n")
        bipart.write("\nTRANSFERTIPS\n")
        for i, value in sorted(clade_transfer_tips_dict.items(), key=lambda (k, v): v[0], reverse=True):
            #i = i.strip()
            i_string = str(i)
            #i_string = '|'.join(map(str, i))
            str4 = ''.join(value[1])
            bipart.write(i_string+"\t"+str(clade_transfer_tips_dict[i][0])+"\t"+str4+"\n")
        ####this should work. pass in a bipartitions_table string.
        transfer_tips = []
        for item in clade_transfer_tips_dict:
            if clade_transfer_tips_dict[item][0] > 80:
                transfer_tips.append(str(item))
       # bipart.write("species_tree_dict\n")
        #for node in gene_tree_dict:
         #   bipart.write("\n"+node+"\t"+' '.join(gene_tree_dict[node]))
        bipart.write("\n\nTIPS THAT RECIEVED TRANSFER AT LEAST 80% of the time: \n")
        for receive_tip in transfer_tips:
            bipart.write(receive_tip+"\n")
        #write the stuff
        #append to the bipart table all the things
        
    print("tips that received the gene via transfer > 80 bs support in the subtree")
    print(transfer_tips)
    #this summarized how many times each thing is a transfer, loss, duplication, or speciation. done on NODES eg m1 m2 m3 m4 m5 m6
    #PROBLEM: the nodes are going to be different in each case. so this is going to have to replace nodename with the tips it contains and get really long. UGH.
    #it will work if the replacement is standardized, but it will be UGLY af.
    summary_dict = {}
    for key in clade_new_transfer_dict:
        summary_dict[key] = [clade_new_transfer_dict[key][0]]
    for key in transfer_recipient_dict:
        if key in summary_dict:
            summary_dict[key].append(clade_new_transfer_dict[key][0])
        else:
            summary_dict[key] = [transfer_recipient_dict[key][0]]
    for key in clade_dup_dict:
        if key in summary_dict:
            summary_dict[key].append(clade_dup_dict[key][0])
        else:
            summary_dict[key] = [0, clade_dup_dict[key][0]]
    for key in summary_dict:
        try:
            #try returning the second value (dup number) - if not there, add a 0
            a = summary_dict[key][1]
        except:
            summary_dict[key].append(0)
    for key in clade_loss_dict:
        if key in summary_dict:
            summary_dict[key].append(clade_loss_dict[key][0])
        else:
            summary_dict[key] = [0, 0, clade_loss_dict[key][0]]
    for key in summary_dict:
        try:
            #try returning the second value (dup number) - if not there, add a 0
            a = summary_dict[key][2]
        except:
            summary_dict[key].append(0)
    for key in clade_speciation_dict:
        if key in summary_dict:
            summary_dict[key].append(clade_speciation_dict[key][0])
        else:
            summary_dict[key] = [0, 0, 0, clade_speciation_dict[key][0]]
    for key in summary_dict:
        try:
            #try returning the second value (dup number) - if not there, add a 0
            a = summary_dict[key][3]
        except:
            summary_dict[key].append(0)
            #print("weird this one has 0 speciation?: "+key)
    print("making summary table: "+summary_table)
    with open(summary_table, "w") as summ:
        summ.write("transferrecip \t transfer_rep_node \t duplication \t loss \t speciation \t node_strings")
        for clade_or_tip in summary_dict:
            a = (summary_dict[clade_or_tip])
            if a == 4:
                summ.write("\n"+"0\n"+str(summary_dict[clade_or_tip][0])+"\t"+str(summary_dict[clade_or_tip][1])+"\t"+str(summary_dict[clade_or_tip][2])+"\t"+str(summary_dict[clade_or_tip][3])+"\t"+clade_or_tip)

            elif a == 5:
                summ.write("\n"+str(summary_dict[clade_or_tip][0])+"\t"+str(summary_dict[clade_or_tip][1])+"\t"+str(summary_dict[clade_or_tip][2])+"\t"+str(summary_dict[clade_or_tip][3])+"\t"+str(summary_dict[clade_or_tip][4])+"\t"+clade_or_tip)
    print(bipartitions_table)
    print(summary_table)
    return bipartitions_table, transfer_tips

#this created lists w Y and N
def dictionary_check_bipart(clade_dict, tuples_list):
    #calde_dict is an empty dictionary
    #tuples list is a list of [loss,loss,loss] or [(transd,transr),(trd, trr),(trd,trr)]
    bs_num = 0
    for ctd in clade_dict:
        #track how many N already should exist for use later.
        bs_num = len(clade_dict[ctd][1])
        err = "yes"
        #look to see if we should add a Y or N to each member of the clade dictionary.
        for tt in tuples_list:
            #if transfer exists in the dictionary and in current bootstrap, add a Y and a number.
            if tt == ctd:
                clade_dict[tt][0] += 1
                clade_dict[tt][1].append("Y")
                err = "no"
        #if transfer in dictionary but not in current bootstrap, add a no.
        if err == "yes":
            clade_dict[ctd][1].append("N")


            ####PROBLEM sometimes TUPLES LIST and CLADE DICT KEYS are not the same.
            ###
        
    #if the transfer is NOT in the dictionary, add it with the appropriate number of NOs.
    for tt in tuples_list:
        if tt in clade_dict:
            pass
        else:
            inner_list = []
            for item in range(bs_num):
                inner_list.append("N")
            inner_list.append("Y")
            clade_dict[tt] = [1,inner_list]
    return clade_dict

#____________RESET


def RunRangerOnInputs_boots(ranger_names, projectname, i):
    list_ranger_outs = []
    #ranger names is just the base-name. 1-5
    #for item in ranger_names:
        #output_name = item[:-6]+"Out.txt"
        #list_ranger_outs.append(output_name)
    list_ranger_outs_to_do = []
    ranger_names_done = []
    ranger_names_to_do = []
    allhome = True
    #check if they are all already done
    for item in ranger_names:
        a = os.path.isfile(str(i)+item[:-6]+"Out.txt")
        if a is True:
            ranger_names_done.append(item)
            list_ranger_outs.append(item[:-6]+"Out.txt")
        else:
            print("Not Pre-Done: missing: "+item)
            ranger_names_to_do.append(item)
            list_ranger_outs_to_do.append(item[:-6]+"Out.txt")
            list_ranger_outs.append(item[:-6]+"Out.txt")
            allhome = False
    if allhome is True:
        print("the ranger analysis was already completed!")
        return list_ranger_outs
    else:
        print("running ranger on inputs for the base files:")
        print(list_ranger_outs_to_do)
    
    move_ran_cluster_boots(ranger_names_to_do, i)
    #modify to do 500 instead of 5.
    run_ranger_on_cluster_boots(ranger_names_to_do, list_ranger_outs_to_do, projectname, i)
    #modify
    move_home_ranger_cluster_boots(list_ranger_outs_to_do, i, projectname)
    #modify to use i - 500 instead of 5
    #os.system("rangeru -i "+item+" -o "+output_name)
    
    print(list_ranger_outs)
    return list_ranger_outs

def move_ran_cluster_boots(ranger_names, i):
    clus_path = "Ranger/"
    i = int(i)
    for item in ranger_names:
        os.system("scp *"+item+" "+clus_head+clus_path)
       # for num in range(i):
        #    num+=1
         #   os.system("scp "+str(num)+item+" "+clus_head+clus_path)
    print("finished moving items to cluster")



def gen_ranger_corr_boots(ran_in, ran_out, indexname):
    print(len(ran_in))
    print(len(ran_out))
    print(ran_out)
    with open(indexname, "w") as corr:
        for i in range(len(ran_in)):
            corr.write(str(i+1)+" "+ran_in[i]+" "+ran_out[i]+"\n")
    return indexname

###this needs to be modified to run 1-100 all at once, so there will only by ~6 jobs.

def gen_ranger_script_boots(n, indexname, scriptfile, bootnum):
    if int(bootnum) != 100:
        print("error in number of bootstraps")
    a =  """#!/bin/bash                                                                             
#SBATCH -p sched_mit_g4nier                                                                           
#SBATCH -t 7-00:00:00    
#SBATCH -J Ranger   
#SBATCH -o Ranger.out                                                                                 
#SBATCH --array=1-"""+n+"""

. /etc/profile.d/modules.sh
module add engaging/openmpi/1.8.8

MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
THE_OUTPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $3}' )
THEROOT=${THE_INPUT_FILE/In.txt/}
OPT=OptRoot.txt
OPTOUT=$THEROOT$OPT

RAN=$PostOpt_RangerIn.txt
RANIN=$THEROOT$RAN

for i in {1..100}; do {
./OptRoot.linux -i $i$THE_INPUT_FILE -o $i$OPTOUT
sleep 1m
mpirun python OptRootParser.py -i $i$OPTOUT -o $i$RANIN -s $i$THE_INPUT_FILE
sleep 1m
./Ranger-DTL.linux -i $i$RANIN -o $i$THE_OUTPUT_FILE
} ; done

exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile

def move_home_ranger_cluster_boots(list_ranger_outs, bootnum, projectname):
    #print("currently set to not move anything home for a test.")
    #return "done"

    movehome = []
    direct = os.getcwd()
    clus_path = "/Ranger"
    list_all_ranger_outs = []
    for thing in list_ranger_outs:
        for iteration in range(bootnum):
            iteration +=1
            list_all_ranger_outs.append(str(iteration)+thing)
    
    for i in list_all_ranger_outs:
        movehome.append(i)
    finished = False
    iteration = 0
    while finished is not True:
        #try and move each home.
        if len(movehome) < 100:
            for filename in movehome:
                print("scp "+clus_head[:-1]+clus_path+"/"+filename+" "+direct)
                os.system("scp "+clus_head[:-1]+clus_path+"/"+filename+" "+direct)
        else:
            print("OR WAS IT")
            os.system("scp "+clus_head[:-1]+clus_path+"/*"+projectname+"*RangerOut* "+direct)
        finished = "yes"
        for item in list_all_ranger_outs:  
            #see if everything got moved home.
            exists = os.path.isfile(item)
            #print(exists)
            #print(item)
            if exists is True:
                if item in movehome:
                    movehome.remove(item)
            else:
                finished = False
                print(item+" not found")
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            #wait five minutes and then try again.
            iteration += 1
            if iteration < 10:
                print("checking.... some ranger outputs do not exist yet. sleeping for 5 minutes.")
                time.sleep(300)
            elif iteration > 100:
                a = raw_input("type resume to resume, skip to skip to next step")
                if a == "skip":
                    return "done"
            else:
                print("checking.... some ranger outputs do not exist yet. sleeping for 15 minutes.")
                time.sleep(900)
     
           

def Separate_Bootstraps(boots_file):
    #input: a single bootstrap file
    list_of_boot_trees = []
    with open(boots_file) as bootfile:
        for line in bootfile:
            list_of_boot_trees.append(line)
            
    #output: a list where each element is one bootstrap
    return list_of_boot_trees

def Check_Boots_Exist(list_of_big_clades, base_dir):
    good = []
    bad = []
    for item in list_of_big_clades:
        try:
            a = item.boot_gene_file
        except:
            print("No bootstrap gene file associated with clade:")
            print(item)
            print(item.prefix)
            print(item.boot_gene_file)
            bad.append(item)
            continue
        #see if the file is here. (this will work if a is a full filename)
        b = os.path.isfile(a)
        
        if b is True:
            #findable as is
            good.append(item)
        else:
            #check base_dir+/gene/boots/thing
            projectname = item.projectname
            b = os.path.isfile(base_dir+"/"+projectname+"/Gene_Trees/boots/"+a)
            if b is True:
                item.boot_gene_file = base_dir+"/"+projectname+"/Gene_Trees/boots/"+a
                good.append(item)
            else:
                print("where am i:")
                print(os.getcwd())
                print("looking for:")
                print(base_dir+"/"+projectname+"/Gene_Trees/boots/"+a)
                raise SystemExit
                d = a.split("/")
                e = d[-1]
                c = os.path.isfile(projectname+"/Gene_Trees/boots/"+e)
                if c is True:
                    item.boot_gene_file = projectname+"/Gene_Trees/boots/"+e
                    good.append(item)
                else:
                    print("Couldn't find a bootstraps file for:")
                    print(item.prefix)
                    print(item.boot_gene_file)
                    bad.append(item)
    return good, bad

def Check_Species_Exist(list_of_big_clades, bad, base_dir):
    #this is checking that the best SPECIES tree exists. besttree raxml.
    #are we checking for the filename or for the string containing the content of the file? besttree_str = ((a,b)(c,d),e);
    #so really, A) check that besttree_str isn't empty.
    #then, for double - check you can insure that the file exists... but honestly thats not really important.
    good = []
    for item in list_of_big_clades:
        if item.besttree_str == "":
            print("no species best_tree was associated with the clade: "+item.prefix)
            bad.append(item)
        else:
            good.append(item)
        # a = item.besttree_str
        # b = os.path.isfile(a)
        # if b is True:
        #     good.append(item)
        # else:
        #     projectname = item.projectname
        #     b = os.path.isfile(base_dir+"/"+projectname+"/Species_Trees/trees/"+a)
        #     if b is True:
        #         item.besttree_str = base_dir+"/"+projectname+"/Species_Trees/trees/"+a
        #         print("b was true")
        #         good.append(item)
        #     else:
        #         d = a.split("/")
        #         e = d[-1]
        #         c = os.path.isfile(base_dir+"/"+projectname+"/Species_Trees/trees/"+e)
        #         if c is True:
        #             print("c was true")
        #             item.besttree_str = base_dir+"/"+projectname+"/Species_Trees/trees/"+e
        #             good.append(item)
        #         else:
        #             print("Couldn't find a besttree file for :")
        #             print(item.prefix)
        #             print(item.besttree_str)
        #             print("issue finding species besttree for: "+e)
        #             raise SystemExit
        #             bad.append(item)
    return good, bad

def run_ranger_on_cluster_boots(ranger_in_names, ranger_out_names, projectname, i):
    indexname = "Ranger_"+projectname+"_Corr.txt"
    scriptfile = "Ranger_"+projectname+"_Sc.sh"
    a = gen_ranger_script_boots(str(len(ranger_in_names)), indexname, scriptfile, i)
    b = gen_ranger_corr_boots(ranger_in_names, ranger_out_names, indexname)
    move_ran_cluster([a, b])
#    os.system()
    os.system(ssh_inst+" 'cd ~/Ranger/;echo $PWD;sbatch "+a+"'")



def Make_Species_Tips(boot, subtree_object):
    #boot is a text string of a bootstrap .newick treee
    a = boot
    if subtree_object.fasta_object == "":
        if subtree_object.fasta == "":
            print("this isn't going to work")
            raise SystemExit
        subtree_object.fasta_object = Fasta(subtree_object.fasta)
    b = subtree_object.fasta_object
    c = b.gen_species_lists()
    d = b.gen_gis_list()


    old = b.ids
    newt = a

    newsp = b.species_names
    newgi = b.gis_list
    for item in old:
        index = old.index(item)
        newspeciesname = newsp[index]
        newspeciesname = newspeciesname.replace("_", "XX0XX")
        replacement_name = newspeciesname+"_"+newgi[index]
        newt = newt.replace(item, replacement_name)
    subtree_object.boots_gene_with_species_names.append(newt)
    return newt



# open species tree and save it as dictionary.
# open gene tree and save it as dictionary.

# read the reconciliations: each is defining everything found in the GENE TREE

# GENE TREE:
# leaves. (no info)
# nodes.
# -> speciation + mapping to species tree
# -> transfer + mapping recip and donor to species tree
# -> duplication + mapping to species tree

# so pull:
# 1. list of DTLS mapped to gene tree nodes
# 2. list of DLS T(r) T(d) to species tree nodes

# visualize:
# 1. list of DTLS that are RELEVENT TO THE BESTTREE
# 2. list of DTLS that are closest-match-mapped onto the BESTTREE
# 3. list of DLS T(R/D) mapped to species tree.

# for context, all need to be normalized with bootstrap %???

#this is a copy so feel free to fuck with it for now.

#rangerfiles is a list of all the files to read. eg 1Cyanos 2Cyanos 3Cyanos if we have 3 bootstraps of cyanos.
#biparitions table is a string, used to write the unique result.
#update june29
def Parse_Ranger_Clade_BOTHWAYS(rangerfiles, bipartitions_table):

    #this counts number of times node represents a transfer (in the case of species tree, we count donor and recip seperately)
    node_transfers_sp_recip = {}
    node_transfers_sp_donor = {}
    node_transfers_ge = {}


    #this counts number of times node represents a loss (currently doesnt work)
    node_loss_sp = {}
    node_loss_ge = {}


    #this counts number of times node represents a duplication
    node_duplication_sp = {}
    node_duplication_ge = {}

    #this counts number of times node represents a speciation
    node_speciation_sp = {}
    node_speciation_ge = {}

    #this counts how many times node is ever found
    node_overall_sp = {}
    node_overall_ge = {}


    #setting up file names i think
    a = bipartitions_table.split("biparts")
    summary_table = a[0]+"summary_table.txt"
    
    bs_num = 1
    #will iterate through all available bootstraps
    for singlebs in rangerfiles:
        #first, we parse a single bootstrap file.
#todo edit this
        results_sp, results_ge = Parse_Single_Ranger_Output_BOTHWAYS(singlebs)

        #parse results to species tree
        speciation_res_sp = results_sp[0]
        transfer_recip_res_sp = results_sp[1]
        transfer_donor_res_sp = results_sp[2]
        loss_res_sp = results_sp[3]
        duplicate_res_sp = results_sp[4]
        overall_res_sp = results_sp[5]

        #parse results to gene tree
        speciation_res_ge = results_ge[0]
        transfer_res_ge = results_ge[1]
        loss_res_ge = results_ge[2]
        duplicate_res_ge = results_ge[3]
        overall_res_ge = results_ge[4]

        #now, update the species_tree dictionaries
        node_transfers_sp_recip = dictionary_check_bipart_BOTHWAYS(node_transfers_sp_recip, transfer_recip_res_sp)
        node_transfers_sp_donor  = dictionary_check_bipart_BOTHWAYS(node_transfers_sp_donor, transfer_donor_res_sp)
        node_loss_sp  = dictionary_check_bipart_BOTHWAYS(node_loss_sp, loss_res_sp)
        node_duplication_sp  = dictionary_check_bipart_BOTHWAYS(node_duplication_sp, duplicate_res_sp)
        node_speciation_sp  = dictionary_check_bipart_BOTHWAYS(node_speciation_sp, speciation_res_sp)
        node_overall_sp  = dictionary_check_bipart_BOTHWAYS(node_overall_sp, overall_res_sp)

        #and the gene tree dicitonaries
        node_transfers_ge  = dictionary_check_bipart_BOTHWAYS(node_transfers_ge, transfer_res_ge)
        node_loss_ge  = dictionary_check_bipart_BOTHWAYS(node_loss_ge, loss_res_ge)
        node_duplication_ge  = dictionary_check_bipart_BOTHWAYS(node_duplication_ge, duplicate_res_ge)
        node_speciation_ge  = dictionary_check_bipart_BOTHWAYS(node_speciation_ge, speciation_res_ge)
        node_overall_ge  = dictionary_check_bipart_BOTHWAYS(node_overall_ge, overall_res_ge)

    #now write the bipartitions file! haha.
    with open(bipartitions_table, "w") as bipart:
        #this one is new and untestedz
        bipart.write("BASED ON SPECIES TREE\n")
        bipart.write("Transfer_Recipients_Sp\n")
        sort_write_the_dict(bipart, node_transfers_sp_recip)
        bipart.write("Transfer_Donors_Sp\n")
        sort_write_the_dict(bipart, node_transfers_sp_donor)
        bipart.write("Loss_Sp\n")
        sort_write_the_dict(bipart, node_loss_sp)
        bipart.write("Duplicate_Sp\n")
        sort_write_the_dict(bipart, node_duplication_sp)
        bipart.write("Speciation_Sp\n")
        sort_write_the_dict(bipart, node_speciation_sp)
        bipart.write("Overall_Sp\n")
        sort_write_the_dict(bipart, node_overall_sp)
        bipart.write("BASED ON GENE TREE\n")
        bipart.write("Transfer_Ge\n")
        sort_write_the_dict(bipart, node_transfers_ge)
        bipart.write("Loss_Ge\n")
        sort_write_the_dict(bipart, node_loss_ge)
        bipart.write("Duplicate_Ge\n")
        sort_write_the_dict(bipart, node_duplication_ge)
        bipart.write("Speciation_Ge\n")
        sort_write_the_dict(bipart, node_speciation_ge)
        bipart.write("Overall_Ge\n")
        sort_write_the_dict(bipart, node_overall_ge)

    #the summary table has been removed for right now becauase i never used it for anything??

    print("completed parsing and wrote to: "+bipartitions_table)

    return bipartitions_table, []

#update june 29
#this created lists w Y and N
def dictionary_check_bipart_BOTHWAYS(node_dict, list_res):
    #node_dict is something like node_loss_dict = [node1:1, node2:4]
    #list_res is something like loss_list = [node1, node2, node4]
    #tracks each bootstrap y/n for each node's presence. also have an overall value (#of times node was yes)
    #0 = total
    #1 = YNNNNYYYYYNNN
    bs_num = 0
    for each_node in node_dict:
        #track how many N already should exist for use later.
        bs_num = len(node_dict[each_node][1])
        err = "yes"
        #look to see if we should add a Y or N to each member of the clade dictionary.
        for result_node in list_res:
            #if transfer exists in the dictionary and in current bootstrap, add a Y and a number.
            if result_node == each_node:
                node_dict[result_node][0] += 1
                node_dict[result_node][1].append("Y")
                err = "no"
        #if transfer in dictionary but not in current bootstrap, add a N and no number.
        if err == "yes":
            node_dict[each_node][1].append("N")        
    #if the transfer is NOT in the dictionary, add it with the appropriate number of NOs.
    for result_node in list_res:
        if result_node in node_dict:
            pass
        else:
            inner_list = []
            for item in range(bs_num):
                inner_list.append("N")
            inner_list.append("Y")
            node_dict[result_node] = [1,inner_list]
    return node_dict

#update june29
def sort_write_the_dict(openfile, diction):
    for i, value in sorted(diction.items(), key=lambda (k, v): v[0], reverse=True):
        i_string = str(i)
        str1 = ''.join(value[1])
        openfile.write(i_string+"\t"+str(diction[i][0])+"\t"+str1+"\n")


#update june29 
def Parse_Single_Ranger_Output_BOTHWAYS(rangerfile):

   

    #initial species lists
    speciation_res_sp = []
    transfer_recip_res_sp = []
    transfer_donor_res_sp = []
    loss_res_sp = []
    duplicate_res_sp = []
    overall_res_sp = []

    #initial gene lists
    speciation_res_ge = []
    transfer_res_ge = []
    loss_res_ge = []
    duplicate_res_ge = []
    overall_res_ge = []

     #initialize lists
    results_sp = [speciation_res_sp, transfer_recip_res_sp, transfer_donor_res_sp, loss_res_sp, duplicate_res_sp, overall_res_sp]
    results_ge = [speciation_res_ge, transfer_res_ge, loss_res_ge, duplicate_res_ge, overall_res_ge]

    #set to recognize gene and species trees
    prep = False
    sp_prep = False
    with open (rangerfile) as ranger:
        for line in ranger:
            if "Gene Tree:" in line:
                prep = True
            elif "Species Tree:" in line:
                sp_prep = True
            elif "minimum reconciliation cost" in line:
                continue
            elif prep is True:
                #make the gene tree dict
                prep = False
                gene_tree_dict = Newick_Tree_Nodes_To_Dict(line)
            elif sp_prep is True:
                #make the species tree dict
                sp_prep = False
                species_tree_dict = Newick_Tree_Nodes_To_Dict(line)
            elif "Speciation" in line:
                #store info one speciation
                #m3 = LCA[MeiothermusXX0XXcerbereus_654400501, MeiothermusXX0XXtimidus_517278376]: Speciation, Mapping --> n7
                node_finder = line.split(" ")
                node_ge = node_finder[0]
                node_sp = node_finder[-1]
                node_ge = node_ge.strip()
                node_sp = node_sp.strip()
                #convert and add gene version
                list_of_tips_ge = gene_tree_dict[node_ge]
                str_of_tips_ge = ' '.join(list_of_tips_ge)
                speciation_res_ge.append(str_of_tips_ge)
                #convert and add species version
                list_of_tips_sp = species_tree_dict[node_sp]
                str_of_tips_sp = ' '.join(list_of_tips_sp)
                speciation_res_sp.append(str_of_tips_sp)
            #
             # m8 = LCA[PyrobaculumXX0XXneutrophilum_501318754, PyrobaculumXX0XXoguniense_504114203]: Transfer, Mapping --> n11, Recipient --> n12
            elif "Transfer" in line:
                node_finder = line.split(" ")
                node_ge = node_finder[0]
                node_sp_recip = node_finder[-1]
                node_sp_donor = node_finder[-4]
                node_sp_donor = node_sp_donor.strip(",")
                node_ge = node_ge.strip()
                node_sp_donor = node_sp_donor.strip()
                node_sp_recip = node_sp_recip.strip()

                #convert and add gene version
                list_of_tips_ge = gene_tree_dict[node_ge]
                str_of_tips_ge = ' '.join(list_of_tips_ge)
                transfer_res_ge.append(str_of_tips_ge)
                #convert and add species version
                list_of_tips_sp = species_tree_dict[node_sp_donor]
                str_of_tips_sp = ' '.join(list_of_tips_sp)
                transfer_donor_res_sp.append(str_of_tips_sp)
                list_of_tips_sp = species_tree_dict[node_sp_recip]
                str_of_tips_sp = ' '.join(list_of_tips_sp)
                transfer_recip_res_sp.append(str_of_tips_sp)
                ################OLD BULLSHIT FOR MAPPING ONTO THE GENE BESTTREE
          #       #process the recipint mess
          #       #make a specient_tree_dict
          #       list_of_tips = species_tree_dict[recipient]
          #       #this will give a list of tips a b c that we need to find in one of the tip lists in gene_tree_dict, and extract the key from there.
          #       poten_list = []
          #       #this is working for list_of_tips == a single tip

          #       #right here we have the species tree recipient kept in "list of tips" -- can we just make an output from this?
          #       str_of_tips = ' '.join(list_of_tips)
          #       transfers_mapped_to_species_recips.append(str_of_tips)
   
          #       for potential_clade in gene_tree_dict:
          #           #print(list_of_tips)
          #           a = Is_X_within_Y(list_of_tips, gene_tree_dict[potential_clade])
          #           if a is True:
          #               poten_list.append([len(gene_tree_dict[potential_clade]), potential_clade])
          #               #print("appending: "+potential_clade)
          #               #print(list_of_tips)
          #       if poten_list == []:
          #           print("ERROR")
          #       #   print("we were looking for a clade containing all of: ")
          #      #    print(list_of_tips)
          #     #     print("from the dict:")
          #    #      print(gene_tree_dict)
          #   #       print("error with gene transfer in RunRangerline 584")
          #           #aise SystemExit
          #           continue
                    
          #       b = sorted(poten_list)
          #       LCA = b[0]
          #       LCA = LCA[1]
          #  #    print(LCA)
          # #     print(gene_tree_dict)
          #        #process the node.
          #       list_of_tips = gene_tree_dict[LCA]
          #       str_of_tips = ' '.join(list_of_tips)
          #       new_transfer_list.append(str_of_tips)
            elif "Loss" in line:
                node_finder = line.split(" ")
                node_ge = node_finder[0]
                node_sp = node_finder[-1]
                node_ge = node_ge.strip()
                node_sp = node_sp.strip()
                #convert and add gene version
                list_of_tips_ge = gene_tree_dict[node_ge]
                str_of_tips_ge = ' '.join(list_of_tips_ge)
                loss_res_ge.append(str_of_tips_ge)
                #convert and add species version
                list_of_tips_sp = species_tree_dict[node_sp]
                str_of_tips_sp = ' '.join(list_of_tips_sp)
                loss_res_sp.append(str_of_tips_sp)

    #             do the third thing
            elif "Duplication" in line:
                node_finder = line.split(" ")
                node_ge = node_finder[0]
                node_sp = node_finder[-1]
                node_ge = node_ge.strip()
                node_sp = node_sp.strip()
                #convert and add gene version
                list_of_tips_ge = gene_tree_dict[node_ge]
                str_of_tips_ge = ' '.join(list_of_tips_ge)
                duplicate_res_ge.append(str_of_tips_ge)
                #convert and add species version
                list_of_tips_sp = species_tree_dict[node_sp]
                str_of_tips_sp = ' '.join(list_of_tips_sp)
                duplicate_res_sp.append(str_of_tips_sp)
    #             do the third thing
            elif "Leaf Node" in line:
                continue
            else:
                #this should ignore the uh initial SPecies Tree": and the species tree itseld and any \n lines
                continue
    #print("done with single_ranger_parseing_bothways")
    return results_sp, results_ge
#____________RESET

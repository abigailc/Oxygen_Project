# #!/usr/bin/python

# last edit abigailc@Actaeon on july 6 2017


#we are reorganizing the massive script to use modules
#because it is, quite frankly, getting ridiculous.

#PERSONAL_SETTINGS
ssh_inst = "ssh -l abigailc -i ~/.ssh/id_rsa eofe4.mit.edu"
clus_head = "abigailc@eofe4.mit.edu:/home/abigailc/"
#######gene_queries_for_gene_trees#####

#change this to the modern version. NT.

#Arc_Gene_List = "/Users/abigailc/Dropbox/Phylo/Scripts/Oxy_SS/Archaeal_Ribo.fasta"
#Bac_Gene_List = "/Users/abigailc/Dropbox/Phylo/Scripts/Oxy_SS/Bacterial_Ribo.fasta"
#Euk_Gene_List = "/Users/abigailc/Dropbox/Phylo/Scripts/Oxy_SS/Eukaryal_Ribo.fasta"
#All_Gene_List = "/Users/abigailc/Dropbox/Phylo/Scripts/Oxy_SS/Al_Ribo.fasta"
#New_Gene_List = "/Users/abigailc/Dropbox/Phylo/Scripts/Oxy_SS/New_Ribo.fasta"
#Arc_Gene_List = "/Users/abigailc/Dropbox/Phylo/Scripts/Oxy_SS/New_Ribo.fasta"
#Bac_Gene_List = "/Users/abigailc/Dropbox/Phylo/Scripts/Oxy_SS/New_Ribo.fasta"
#Euk_Gene_List = "/Users/abigailc/Dropbox/Phylo/Scripts/Oxy_SS/New_Ribo.fasta"

All_Gene_List = "/Users/abigailc/Documents/New_MST/Newest_Ribo.fasta"

#experiment?
Arc_Gene_List = All_Gene_List
Bac_Gene_List = All_Gene_List
Euk_Gene_List = All_Gene_List
New_Gene_List = All_Gene_List
#IMPORTS

#CALL SOMETHING TO MAKE SUBTREES, GET CLADE INFORMATION
from oxy_mods.Subtree_Clade_Gen import *

#CALL SOMETHING TO MAKE SPECIES TREES
from oxy_mods.MakeSpeciesTrees import *

#CALL SOMETHING TO MAKE GENE TREES
from oxy_mods.MakeGeneTrees import *

#CALL SOMETHING TO RUN RANGER DTL
from oxy_mods.Run_Ranger import *

#CALL SOMETHING TO DO THE SUBSAMPLING
from oxy_mods.Preform_Subsampling import *

#CLASSES
from oxy_mods.Classes_DTL_Detector import Fasta
from oxy_mods.Classes_DTL_Detector import Subtree
#includes #FASTA and #Subtree

#PARSER

def Write_Info(projectname, tabl):
    information = projectname+"_cladeINFO.txt"
    print("************************************************************")
    for item in tabl:
        print(item)
    with open (information, "w") as info:
        for titem in tabl:
            if type(titem) == str:
                info.write("\n"+titem+"\n")
            else:
                titem_text = titem.get_string()
                info.write(titem_text)
        try:
            info.write(str(args))
        except:
            pass

def Create_Subtrees(subtrees_new_dict, st_trees,inputfasta, projectname, typedict, stdict):
    list_of_subtree_class_objects = []
    list_of_short_clades = []
    list_of_big_clades = []
    #from these messy returns, create objects of class subtree to pass on to the next function.
    #each subtree in sub-new-dict must be represented.
    for item_name in subtrees_new_dict:
        #item should be the subtree NUMBER identifier
        #this adds an object of class Subtree with name = number identifier to the list "list_of_subtree_class_objects"
        this_object = Subtree(item_name)
        list_of_subtree_class_objects.append(this_object)
        
        this_object.set_fasttree(st_trees[item_name])
        #cleaning up the group name to never include lead/trail whitespace or bars.
        string_a = subtrees_new_dict[item_name]
        string_a = string_a.strip("|")
        string_a = string_a.strip()
        type_a = typedict[item_name]
        this_object.set_type(type_a)
        this_object.set_string(string_a)
        this_object.set_tips(stdict[item_name])
        #finished generating the basics
    #load the input fasta file
    print(inputfasta)
    orig_fasta = Fasta(inputfasta)
    orig_fasta.gen_original_lists(inputfasta)
    #create the sub-fasta files and name them something reasonable?
    #create the PREFIXES which will be:
    #string(first six letters) + stnumber (first four numbers) +
    list_of_fasta_to_align = []
    for subtree_object in list_of_subtree_class_objects:
        #generate a list of " keeps" which needs to exactly match up with the sequence IDs as stored in the Fasta object.
        #we have the tips loaded as list from subtree_object.ret_tips
        alltips = []
        for tip in subtree_object.ret_tips():
            newtip = tip.strip()
            newtip = newtip.strip(">")
            newtip = newtip.strip("'")
            alltips.append(newtip)
        #alltips is a list of 'cleaned' tip names from the tips stored in stdict. so it will match those from original .fasta file.
        #this modifies the current tips and seqs lists in orig_fasta to only include those in this specific subtree.
        orig_fasta.extract(alltips)
        #this creates the prefix for this specific subtree
        st_num = subtree_object.ret_name()[:6]
        group_nam = subtree_object.ret_string()[:4]
        prefix = projectname+group_nam+st_num
        #this generates the new .fasta file (literally makes a fasta containing only the sequences seen in this subtree)
        orig_fasta.gen_new_fasta(prefix+"_gene.fasta")
        #this links the newly created fasta with this subtree instance.
        subtree_object.fasta = prefix+"_gene.fasta"
        subtree_object.set_fasta(prefix+"_gene.fasta")
        subtree_object.fastaname = prefix+"_gene.fasta"
        #this links the prefix used with this specific subtree instance.
        subtree_object.set_prefix(prefix)
        #CHECK they have enough tips
        if len(subtree_object.tips) < 5:
            if len(subtree_object.tips) == 0:
                print ("zero?????")
            else:
                print(subtree_object.prefix)
                print(subtree_object.fasta)
            subtree_object.set_fasta_object()
            list_of_short_clades.append(subtree_object)
        else:
            #check the tips have enough DIFFERENT species names.
            subtree_object.set_fasta_object()
            sub_species = subtree_object.fasta_object.gen_species_lists()
            sub_species_reduced = []
            for item in sub_species:
                if item in sub_species_reduced:
                    pass
                else:
                    sub_species_reduced.append(item)
            if len(sub_species_reduced) < 5:
                list_of_short_clades.append(subtree_object)
            else:
                list_of_big_clades.append(subtree_object)
                list_of_fasta_to_align.append(prefix+"_gene.fasta")
    return(list_of_big_clades, list_of_fasta_to_align, list_of_short_clades, orig_fasta)
    
def del_slurm_files(ssh_inst, prefix):
    os.system(ssh_inst+" \' cd Species_Trees;rm -r *"+prefix+"*;rm *"+prefix+"* \'")
    os.system(ssh_inst+" \' cd Gene_Trees;rm -r *"+prefix+"*;rm *"+prefix+"* \'")
    print("All your cluster files are now purged")

def Remove_Extra_Gene_Tree_Species(list_of_big_clades):
    #this needs to remove sequences from the genetree fasta that we could not find species_tree matches for.
    for item in list_of_big_clades:
        fasbase = item.fasta.split(".")[-1]
        #print(item.species_fasta_object)
        #print(item.species_fasta_object.original_ids)
        exclude = []

        #this should be unneccessary and deleted after testing
        try:
            a = item.species_fasta_object.original_ids
        except:
            print (item.prefix)
            print ("has no original IDS? trying gen original lists but maybe this should be a small clade...?")
            item.species_fasta_object.gen_original_lists()
            print(item.species_fasta_object.original_ids)
        ind = 0
        for sp in item.species_list_original:
            if sp in item.species_fasta_object.original_ids:
                pass
            else:
                print("excluding "+sp+" for not being present in species tree")
                exclude.append(ind)
            ind += 1
        if exclude == []:
            pass
            #    os.system("cp "+item.fasta+" "+fasbase+".new")
        else:
            #generate a new .fasta file with the offending entries removed, and add its name.
            os.system("mv "+item.fasta+" "+item.fasta+".old")
            item.old_fasta_name = item.fasta+".old"
            item.old_fasta_object = item.fasta_object
            print(exclude)
            print("len old fasta:"+str(len(item.fasta_object.ids)))
            exclude.sort()
            exclude.reverse()
            for num in exclude:
                item.fasta_object.ids.pop(num)
                item.fasta_object.seqs.pop(num)
            print("For clade : "+item.prefix+" , we had to remove "+str(len(exclude))+" sequences because the species was not found in the generated species tree.")
            item.fasta_object.gen_new_fasta(item.fasta)
            #this printed a new fasta. but did not set it as the new object. lets do that.
            item.fasta = item.fasta
            item.set_fasta_object()
            print("len new fasta:"+str(len(item.fasta_object.ids)))

            #should result in .fasta and .fasta_object being the new version, without anything not in species trees.
                #make this named what the other one should have been/is, and rename that one something else.
        

#this is the main function to call when starting a run of OverallDTL detacter. should be called in the command parser thing
def Run_Detection(stringlist, treefile, artcutoff, fewnum, treesubfile, inputfasta, verbose, projectname = "NA"):
    
    #initialized for file moving later
    start_objects = [treefile, treesubfile, inputfasta]

    #parse the subtrees
    tabl, subtrees_new_dict, st_trees, stdict, typedict = PerformScan_Complex_SubSample(stringlist, treefile, artcutoff, fewnum, treesubfile, verbose)

    #make a bunch of organizational directories that will eventually be used, and save information on this run.
    os.system("mkdir "+projectname+"; cd "+projectname+"; mkdir Species_Trees; mkdir Gene_Trees; cd Species_Trees; mkdir trees; mkdir construction; mkdir boots; mkdir concat; cd -; cd Gene_Trees; mkdir fasta; mkdir muscle; mkdir boots;  mkdir trees")
    Write_Info(projectname, tabl)

    #initializing lists of Subtree objects, gene fastas.
    list_of_big_clades, list_of_fasta_to_align, list_of_short_clades, orig_fasta = Create_Subtrees(subtrees_new_dict, st_trees,inputfasta, projectname, typedict, stdict)
    for item in list_of_big_clades:
        item.projectname = projectname
 
    #making the species trees
    rax_sp_out, cc_files, incomplete_species, boots_species = Overall_Species_Trees(list_of_big_clades, projectname, Euk_Gene_List, Bac_Gene_List, Arc_Gene_List, All_Gene_List)

    #Overall_Species_Trees(list_of_big_clades, projectname, Euk_Gene_List, Bac_Gene_List, Arc_Gene_List, All_Gene_List)

    #list of big clades needs to remove anything that doesn't have a matched species tree.
    #testcase: gast
    to_remove_large = []
    for item in list_of_big_clades:
        try:
            a = item.species_besttree_name
            if a == "":
                print("A IS EMPTY")
        except:
            print("didn't find a species_besttree")
            to_remove_large.append(item)
    print("small clades from MST")
    print(to_remove_large)
    for item in to_remove_large:
        list_of_big_clades.remove(item)
        list_of_short_clades.append(item)
   # print(list_of_big_clades)

    
    #now each thing in big list has a species_list_original and a species_list_after_removal.

    #compare the two; if there are ANY species in the ORIGINAL but not also in the AFTER REMOVAL; those ID and SEQS need to be removed.

    #remove from big_clades gene.fasta file any sequences whose species name returned no hits in the species tree generator.
    Remove_Extra_Gene_Tree_Species(list_of_big_clades)
    
    #GENE TREE HAPPENS
    rax_list, muscle_list, list_of_fasta_to_align, incomplete_gene, boots_gene = Overall_Gene_Trees(list_of_fasta_to_align, projectname, list_of_big_clades)

    #again, check that all things had enough sequences left to be aligned in a reasonable manner.
    #if not, move them to small clades.
    print("removing any clades that could not complete this step")

    print(incomplete_gene)
    for item in incomplete_gene:
        list_of_big_clades.remove(item)
        list_of_short_clades.append(item)
    #print(list_of_big_clades)


    
    #print("raxlist")
    #print(rax_list)
    
    #print(boots_species)
    #print(boots_gene)
    #print("boots lists")
    for item in boots_gene:
        a = os.path.isfile(item)
        if a is False:
            print("couldn't find: "+item+" trying to scp again... NOT IMPLEMENTED")
            pass
        
    #relocating objects to be organized and find-able
    print("Moving objects to better locations...")
    for item in list_of_fasta_to_align:
        os.system("mv "+item+" ./"+projectname+"/Gene_Trees/fasta")
    for item in muscle_list:
        os.system("mv "+item+" ./"+projectname+"/Gene_Trees/muscle")
    for item in cc_files:
        os.system("mv "+item+" ./"+projectname+"/Species_Trees/concat")
        #rax_list doesnt appear to be the correct list of gene_raxml_outputs file
    for item in rax_list:
        #why is _gene.fasta here??
        os.system("mv "+item+" ./"+projectname+"/Gene_Trees/trees")
    for item in rax_sp_out:
        os.system("mv "+item+" ./"+projectname+"/Species_Trees/trees")
    os.system("mkdir "+projectname+"/Species_Trees/boots")
    for item in boots_species:
        os.system("mv "+item+" ./"+projectname+"/Species_Trees/boots")
    os.system("mkdir "+projectname+"/Gene_Trees/boots")
    for item in boots_gene:
        os.system("mv "+item+" ./"+projectname+"/Gene_Trees/boots")
        
    #cleanup. delete extra trees, move folders of junk into misc in case we want them for troubleshooting.
    os.system("cd "+projectname+";rm RAxML_bip*")
    os.system("mv "+projectname+"*_gene.fasta ./"+projectname+"/Gene_Trees/fasta")
    os.system("mv "+projectname+"*INFO.txt ./"+projectname)
    #move start objects to SOD7/
    for s_o in start_objects:
        os.system("mv "+s_o+ " ./"+projectname)
    #move anything else with sod7 to CONSTRUCTION
    os.system("mv "+projectname+"*Rax_Bipart.newick ./"+projectname+"/Species_Trees/trees")
    os.system("mv "+projectname+"* ./"+projectname+"/Species_Trees/construction")
    #copy the start objects back to MAKESPECIESTREE
    #because if we need to re-run this is where they should be.
    for s_o in start_objects:
        os.system("cp ./"+projectname+"/"+s_o+" ./")
    #del_slurm_files(ssh_inst, projectname)
    #commented out for now b/c correlation keeps failing and then i need to rerun everything.
    print("Now it's time to run everything through rangerdtl and it's parser!")
    chosen_seqs = []
    #move_slurm_files(ssh_inst, projectname)
    print("Creating input files with matching species_names to make ranger_in files from")
    for subtree_object1 in list_of_big_clades:
        #make sure that the gene_tree and species_tree files are where we expect them to be.
        #RAxML_bestTree.SOD7Inse###_CC is missing, errored out.
        subtree_object1.set_gene_tree_species_tips()
        subtree_object1.set_alignment_with_species_tips()
       # try:
        gene_tree = subtree_object1.ret_gene_tree_species_tips()
           # print(gene_tree)
        with open("./"+projectname+"/Gene_Trees/trees/"+subtree_object1.ret_prefix()+"_gene_tree_sp_tips.newick", "w") as new:
            new.write(gene_tree)
            #species_tree = subtree_object1.ret_species_tree()
            ####use BESTTREE instead of bipart species tree.
        #subtree_object1.species_besttree_name
        ex = os.path.isfile(subtree_object1.species_besttree_name)
        if ex is False:
            #check is the besttree file was moved during file organizing.
            newsbn = "./"+projectname+"/Species_Trees/trees/"+subtree_object1.species_besttree_name
            ex2 = os.path.isfile(newsbn)

            if ex2 is True:
                subtree_object1.species_besttree_name = newsbn
            else:
                print ("couldnt find species_besttree_name in "+newsbn)

        with open(subtree_object1.species_besttree_name) as st:
            for line in st:
                if line == "":
                    pass
                else:
                    species_tree = line
                    print("setting species tree")
                    #print(species_tree)
                    subtree_object1.besttree_str = species_tree
    #is orig_fasta a str or a fasta_object?
    orig_fasta = Fasta(inputfasta)
    orig_fasta.gen_original_lists(inputfasta)
    print(len(orig_fasta.ids))
    Master_Ranger_SS(list_of_big_clades, projectname, list_of_short_clades, orig_fasta)
    subfasta = Master_Preform_Subsampling(list_of_big_clades, list_of_short_clades, orig_fasta, projectname)


    print("cleanup")
    os.system(ssh_inst+" \' cd Ranger;rm *"+projectname+"* \'")


    print("the newly subsampled file is at: "+subfasta)
    print("Everything done.")
    print("lets try generating a species tree for the SS file and running it through everything")
    print("how do we root this species tree???????")
    #calling fnx right below
    run_on_ss_fasta(subfasta, projectname)
    print("really done")

#this is the function that is called when doing the FINAL PASS of an OverallDTL detector run. it should be fed the subsampled, final dataset.
def run_on_ss_fasta(subfasta,projectname):
    
    a = os.getcwd()
    print(a)
    print("^ THE CURRENT DIRECTORY")
    
    #os.system("mkdir "+projectname+"; cd "+projectname+"; mkdir Species_Trees; mkdir Gene_Trees; cd Species_Trees; mkdir trees; mkdir construction; mkdir boots; mkdir concat; cd -; cd Gene_Trees; mkdir fasta; mkdir muscle; mkdir boots;  mkdir trees")
  
    #generate the subtree and fasta
    SS_Subtree = Subtree("SubSampled")
    SS_Subtree.fasta = subfasta
    SS_Subtree.fasta_object = Fasta(subfasta)
    SS_Subtree.fasta_object.gen_original_lists(subfasta)

    SS_Subtree.fasta_object.prefix = projectname+"_SS"
    SS_Subtree.prefix = projectname+"_SS"
    SS_Subtree.projectname = projectname+"_SS"
    projectname = projectname+"_SS"
    os.system("mkdir "+projectname+"; cd "+projectname+"; mkdir Species_Trees; mkdir Gene_Trees; cd Species_Trees; mkdir trees; mkdir construction; mkdir boots; mkdir concat; cd -; cd Gene_Trees; mkdir fasta; mkdir muscle; mkdir boots;  mkdir trees; cd ..; cd ..")
    a = os.getcwd()
    print(a)
    print("^ THE CURRENT DIRECTORY")

    # for now, use Bacterial queries to make the SubSampled Species Tree
    cat_file = All_Gene_List


    #optional NEW CAT FILE... including queries from multiple candidates.
    #make a species tree from this. how do I root or add loss candidates here???
    #currently: no rooting implemented.
    #if you want to add rooting, add a variable containing the rooting species name, and pass it through:
    #species tree for the subsampled fasta
    #run everything SS
    #run raxml on cluster
    #instead of "FALSE" going into run raxml - pass the tip name in a list [species_name]
    print("making a species tree")
    rax_out, cc_files, incomplete, boots_species = Species_Tree_For_The_Subsampled_Fasta(SS_Subtree, SS_Subtree.projectname, cat_file)
    print("finished Species_Tree_for_ss_fasta")
    print(rax_out)
    rax_b = rax_out[0].split(".")
    rax_b.pop(0)
    rax_be = "".join(rax_b)
    rax_best = "RAxML_bestTree."+rax_be
    print(rax_best)
    SS_Subtree.set_species_tree_name(rax_best)
    SS_Subtree.besttree_str = rax_best
    print("besttree is :"+rax_best)

    Remove_Extra_Gene_Tree_Species([SS_Subtree])

    print("Making a gene tree")
    rax_list, muscle_list, list_of_fasta_to_align, incomplete, boots_gene = Overall_Gene_Trees([SS_Subtree.fasta], SS_Subtree.projectname, [SS_Subtree])
    print(rax_list)
    SS_Subtree.set_gene_tree_name(rax_list[0])

    print("^ is rax_list")
    SS_Subtree.set_gene_tree_species_tips()
    print(boots_gene)

    #the boots file isn't being downloaded or assigned - problem in MGT or here?


    #moving things
    print("Moving objects to better locations...")
    for item in list_of_fasta_to_align:
        os.system("mv "+item+" ./"+projectname+"/Gene_Trees/fasta")
    for item in muscle_list:
        os.system("mv "+item+" ./"+projectname+"/Gene_Trees/muscle")
    os.system("mv "+projectname+"*gene.fasta_Muscle.fasta ./"+projectname+"/Gene_Trees/muscle")
    #in case... whatever.
    for item in cc_files:
        os.system("mv "+item+" ./"+projectname+"/Species_Trees/concat")
        #rax_list doesnt appear to be the correct list of gene_raxml_outputs file
    for item in rax_list:
        #why is _gene.fasta here??
        os.system("mv "+item+" ./"+projectname+"/Gene_Trees/trees")
    for item in rax_out:
        os.system("mv "+item+" ./"+projectname+"/Species_Trees/trees")
    os.system("mkdir "+projectname+"/Species_Trees/boots")
    for item in boots_species:
        os.system("mv "+item+" ./"+projectname+"/Species_Trees/boots")
    os.system("mkdir "+projectname+"/Gene_Trees/boots")
    for item in boots_gene:
        os.system("mv "+item+" ./"+projectname+"/Gene_Trees/boots")
        
    #cleanup. delete extra trees, move folders of junk into misc in case we want them for troubleshooting.
    os.system("cd "+projectname+";rm RAxML_bip*")
    os.system("mv "+projectname+"*_gene.fasta ./"+projectname+"/Gene_Trees/fasta")
    os.system("mv "+projectname+"*INFO.txt ./"+projectname)

    #end moving things


    print("running rangerdtl")
    Master_Ranger_Single(SS_Subtree, SS_Subtree.projectname)
    
    print("uhh so you should probably have a ranger output now. and next maybe attach it to trees visually???")
    #






#parser
    
#todo: simplify this, we don't need all the options I think.

if __name__ == "__main__":

    #IMPORTS again?
    #CALL SOMETHING TO MAKE SUBTREES, GET CLADE INFORMATION
    from oxy_mods.Subtree_Clade_Gen import *

    #CALL SOMETHING TO MAKE SPECIES TREES
    from oxy_mods.MakeSpeciesTrees_oxy_branch import *

    #CALL SOMETHING TO MAKE GENE TREES
    from oxy_mods.MakeGeneTrees_oxy import *

    #CALL SOMETHING TO RUN RANGER DTL
    from oxy_mods.Run_Ranger import *

    #CALL SOMETHING TO DO THE SUBSAMPLING
    #from oxy_mods.???? import *

    #CLASSES
    from oxy_mods.Classes_DTL_Detector import Fasta
    from oxy_mods.Classes_DTL_Detector import Subtree
    
    print("Running in terminal")  
    import sys
    import argparse
    import os
    import re
    import time

    #AN EXAMPLE CALL: $python Overall_DTL_Detector.py ~/Documents/New_MST/ example_aligned_fasttree.nexus -r 4 -m 5 (-t if treefile) (-ss or -c ???) 
    #-f example.fasta -v -n EXAMpleproject1
 
    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in (where .nex resides)")
    parser.add_argument("NEXUS_TREE", type=str, help="type the name of your .nex file")
    parser.add_argument("-st", "--string", action = "store", help="to parse only given string, or strings space-seperated within quotes")
    parser.add_argument("-r", "--rank", action = "store", help="to parse all strings of given depth in tree. Euk|Meta|Chor|Mam = 1|2|3|4. Type one number or several space-sep in quotes.")
    parser.add_argument("-po", "--polyclade", type = float, action = "store", default = 0.1, help = "grouptips(in-clade/in-tree) required in subtree to count as a distinct clade. default is 0.1")
    parser.add_argument("-m", "--minimum", action = "store", default = 10, help="number of tips below which a group will not be considered. default is 10")
    parser.add_argument("-t", "--treesfile", action = "store", default = "NA", help="if you already have a subtrees file and dont want to rerun it")
    parser.add_argument("-ss", "--subsampling", action = "store_true", help="toggles one rank lower subsampling within each clade")
    parser.add_argument("-d", "--detection_program", action = "store_true", help="toggles run the detection program(main point)")
    parser.add_argument("-f", "--fasta", action = "store", help="if using subsampling, provide a fasta file to pull sequences from")
    parser.add_argument("-v", "--verbose", action = "store_true", help="prints more information - for debugging mostly.")
    parser.add_argument("-n", "--nameofproject", action = "store", help="name for youyr project - will be included in output file names. short. eg: SQMO")
    parser.add_argument("-pr", "--progress_read", action = "store", default = False, help="input progress file name if need be")
    parser.add_argument("-fp", "--final_pass", action = "store", default = False, help="input a fasta file. i will make a species tree (currently unrooted) and do ranger DTL on the BEST raxml genetree. output biparts + (future) tree w transfer values at nodes)")
    
    args = parser.parse_args()
    
#making input list
    try:
        os.chdir(args.directory)
    except:
        print ("didn't change dir")

    if args.final_pass != False:
        #do the final pass
        #and then close
        
        a = run_on_ss_fasta(subfasta)

        print("done")
        raise SystemExit
    if ".newick" in args.NEXUS_TREE:
        try:
            os.system("ConvertPhylo . "+args.NEXUS_TREE+" newick nexus")
            name, newick = args.NEXUS_TREE.split(".")
            args.NEXUS_TREE = name+".nexus"
        except:
            pass

    #how are we parsing the clades - by rank or by string presencce/absence?
    if str(args.string) == "None":
        if str(args.rank) == "None":
            print("You need to either specify a string with -s, or a rank-depth with -r")
            raise SystemExit
        arank = args.rank.split()
        stringlist = ["PerformGetRankFromSubtreesFile", args.rank]
    elif str(args.rank) == "None":
        stringlist = args.string.split()
    else:
        print("You may not use both -s and -r at the same time (yet), sorry!")

    #how are we creating the subtrees-by-(rank or string) file? what paramenters?
    artcutoff = 1.0 - float(args.polyclade)
    print("Distincty clade = "+str(args.polyclade)+" of group-tips")
    print("Minimum group tips for clade consideration: "+str(args.minimum))
    
    #complex subsampling needs to exist. creates a .raxml of the subtree and also makes a raxml species tree for those taxa.
    if args.subsampling == True:
        print("Beginning Smart SubSampling!")
        SubSampling_Master(stringlist, args.NEXUS_TREE, artcutoff, int(args.minimum), args.treesfile, args.fasta, args.verbose)
        #call master subsampling
    if args.detection_program == True:
        if args.progress_read != False:
            list_big, list_small = Read_Progress_File(args.progress_read, args.nameofproject)
        else:
            list_big, list_small = Run_Detection(stringlist, args.NEXUS_TREE, artcutoff, int(args.minimum), args.treesfile, args.fasta, args.verbose, args.nameofproject)
        #seqs_to_keep = Master_Ranger_SS(list_big, args.nameofproject, list_small, args.fasta)
        
    else:
        tabl = PerformScan(stringlist, args.NEXUS_TREE, artcutoff, int(args.minimum), args.treesfile, args.verbose)      
        information = args.NEXUS_TREE+"_INFO.txt"
        for item in tabl:
            print(item)
        with open (information, "w") as info:
            for titem in tabl:
                if type(titem) == str:
                    info.write("\n"+titem+"\n")
                else:
                    titem_text = titem.get_string()
                    info.write(titem_text)
            info.write(str(args))
        print("finished writing table info")
    print("Analysis Finished!")
    print("************************************************************")
          

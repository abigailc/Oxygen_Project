#!/usr/bin/python

#abigailc@Actaeon Jan 4 2017
#edit july 6 2017
#this is a module to run gene trees on the cluster.
##SET
names_file = "/Users/abigailc/Documents/Taxonomy_Stuff/taxdump/names.dmp"
nodes_file = "/Users/abigailc/Documents/Taxonomy_Stuff/taxdump/nodes.dmp"

######### PERSONAL_SETTINGS #########
ssh_inst = "ssh -l abigailc -i ~/.ssh/id_rsa eofe4.mit.edu"
clus_head = "abigailc@eofe4.mit.edu:/home/abigailc/"
Path_Blast = "/Users/abigailc/blast/"

import time
import os
from oxy_mods.Classes_DTL_Detector import Fasta
def gen_correlate_file(list_of_input_files, corr_file):
    #this should be in form
    #1 name
    #2 name
    #3 name
    #requires 1: list of files 2. name for corr_file.
    i = 0
    with open(corr_file, "w") as corr:
        for item in list_of_input_files:
            i += 1
            #make sure the \n spacing works correctly.
            corr.write(str(i)+" "+item+"\n")
    return corr_file

def move_to_cluster(list_of_files, clus_path):
    #requires os.system
    #requires scp(?)
    #do the thing
    for item in list_of_files:
        os.system("scp "+item+" "+clus_head+clus_path)
    print("Finished moving files to cluster in place:"+clus_path)
#todo make cluster stuff into it\s own module?

#takes input of all of the fastas you want to align, and the group named
#here end-file-list is actually just a list of all the .fasta files you want to align.
#here "prefix" refers to the project name, not the individual prefixes. Confusing, right? I know.
def muscle_align_on_cluster(end_file_list, prefix):
    #this creates dir you will use on the cluster.
    aligned_list = []
    remove_list = []
    for item in end_file_list:
        testf = Fasta(item)
        testf.gen_original_lists(item)
        a = len(testf.ids)
        if a < 5:
            remove_list.append(item)
        else:
            aligned_list.append(item+"_Muscle.fasta")
    for athing in remove_list:
        end_file_list.remove(athing)
    if end_file_list == []:
        print("Too few files for to make tree")
        return "NA"
    print("removed : "+str(len(remove_list))+"files for having too few sequences. will align "+str(len(end_file_list))+" files.")
    check_directory_existance(prefix, ssh_inst)
    clus_path = "/Gene_Trees/"+prefix
    a = gen_muscle_script(prefix+"_Gene_Sc.sh", "~"+clus_path+"/"+prefix+"_Gene_Corr.txt", str(len(end_file_list)), prefix+"job")
    b = gen_correlate_file(end_file_list, prefix+"_Gene_Corr.txt")
    end_file_list.append(a)
    end_file_list.append(b)
    direct = os.getcwd()
    move_to_cluster(end_file_list, clus_path)
    print("everything should be generated and on the cluster")
    n = str(len(end_file_list))
    x = str(int(n)-2)
    print("there are "+x+" files that should be aligning right now")
    os.system(ssh_inst+" 'cd ~/Gene_Trees/"+prefix+";echo $PWD;sbatch "+a+"'")
    finished = "start"
    #to see if the run is complete, see if each new file has been generated. check every 5 minutes for muscle.
    #initialize list of things to move home. initially equal to aligned_list.
    movehome = []
    for i in aligned_list:
        movehome.append(i)
    while finished is not True:
        #try and move each home.
        for filename in movehome:
            os.system("scp "+clus_head[:-1]+clus_path+"/"+filename+" "+direct)
        finished = "yes"
        for item in aligned_list:
           
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
            print("checking.... some alignment outputs do not exist yet. sleeping for 5 minutes.")
            time.sleep(300)
           
    #now, concatenate the output file. when doing large batches, might want to move each to new folder.... but maybe not for now.
    print("Your files have been aligned! They are located at "+direct)
    return aligned_list


#cc_file_list is your list of aligned files.
#prefix is actually projectname.
def raxml_run_on_cluster(cc_file_list, prefix):
    #this creates dir you will use on the cluster.
    tree_list = []
    print(cc_file_list)
    for item in cc_file_list:
        if item == "NA":
            pass
        #this is going to be something different... like raxml_bipartitions.blah
        alist = item.split(".")
        athing = alist[0]
        tree_list.append("RAxML_bipartitions."+athing)
        tree_list.append("RAxML_bestTree."+athing)
    print(tree_list)
    check_directory_existance(prefix, ssh_inst)
    clus_path = "/Gene_Trees/"+prefix
    a = gen_raxml_script(prefix+"_Rax_Gene_Sc.sh", "~"+clus_path+"/"+prefix+"_Rax_Gene_Corr.txt", str(len(cc_file_list)), prefix+"job")
    b = gen_correlate_file(cc_file_list, prefix+"_Rax_Gene_Corr.txt")
    cc_file_list.append(a)
    cc_file_list.append(b)
    direct = os.getcwd()
    move_to_cluster(cc_file_list, clus_path)
    n = str(len(cc_file_list))
    print("everything should be generated and on the cluster. filenumber: "+n+" starting raxml.")
    os.system(ssh_inst+" 'cd ~/Gene_Trees/"+prefix+";echo $PWD;sbatch "+a+"'")
    finished = "start"
    #to see if the run is complete, see if each new file has been generated. check every 5 minutes for muscle.
       #initialize list of things to move home. initially equal to aligned_list.
    movehome = []
    #we are looking for stuff we shouldn't!
    for i in tree_list:
        movehome.append(i)
    while finished is not True:
        #try and move each home.
        for filename in movehome:
            os.system("scp "+clus_head[:-1]+clus_path+"/"+filename+" "+direct)
        finished = "yes"
        for item in tree_list:
            #see if it got moved home.
           
            exists = os.path.isfile(item)
            if exists is True:
                if item in movehome:
                    movehome.remove(item)
                
            else:
                finished = False
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            #wait an hour and then try again.
            #this value should vary based on the size of the tree we are running. 125 tips - > five minutes.
            #thousands? will take much longer, use maybe an hour between checks.
            print("Waiting 15 mins.")
            print("Still need to bring home:")
            print(movehome)
            time.sleep(900)
           
         
    c = tree_list
    print(tree_list)
    print("RAXML finished, tree(s) should exist locally")
    return c


def gen_raxml_script(scriptfile, indexname, n, Jobname):
    #currently assuming you are running in the dir that files are in and should be returned to.
    #thats it. just print the script. return its filename, which will need to be added to list of things to be moved to the cluster.
    

##example script
    a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 7-00:00:00    
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20                                                                                 
#SBATCH -J RAX"""+Jobname+"""   
#SBATCH -o RAX"""+Jobname+""".out                                                                                         
#SBATCH --array=1-"""+n+"""

. /etc/profile.d/modules.sh
module add engaging/RAxML/8.2.9
##gets my array id, which is needed for use below. this will be, i think, a number like 1,2,3 etc
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
echo $MY_ARRAY_ID

## given an index file formatted                                                                        
## <index> <filename>                                                                                   
## produce the filename for given index                                                                 
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
echo $THE_INPUT_FILE
NEW=${THE_INPUT_FILE%%.*}
echo $NEW
  
raxmlHPC-PTHREADS-AVX -T 20 -f a -m PROTGAMMALGF -p 12345 -x 12345 -#100 -n $NEW -s $THE_INPUT_FILE         

exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile



   ####needs editing
def check_directory_existance(prefix, ssh_inst):
    os.system(ssh_inst+" \' mkdir Gene_Trees;cd Gene_Trees;mkdir "+prefix+"\'")

def move_slurm_files(ssh_inst, prefix):
    os.system(ssh_inst+" \' mkdir Species_Trees;cd Species_Trees;mkdir "+prefix+"\'")
    os.system(ssh_inst+" \' cd Species_Trees;cd "+prefix+"; mv "+prefix+"* "+prefix+"\'")
    print("All your cluster files are now within /home/username/Species_Trees/"+prefix)
    
def gen_muscle_script(scriptfile, indexname, n, Jobname):
    #currently assuming you are running in the dir that files are in and should be returned to.
    # direct = os.getcwd()
    # host = socket.gethostname()
    # user = getpass.getuser()
    # addr = user+"@"+host+":"+direct
    #figure out how many need to be run. n = len(listoffilesinindex). n
    #figure out a name - scriptfile
    #figure out path to index file. indexname
    #thats it. just print the script. return its filename, which will need to be added to list of things to be moved to the cluster.
    addr = "PLACEHOLDER"
    

##example script
    a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 0-20:00:00                                                                                   
#SBATCH -J Mus"""+Jobname+"""                                                                                         
##SBATCH -o Mus"""+Jobname+""".out
#SBATCH --array=1-"""+n+"""

. /etc/profile.d/modules.sh
module add engaging/muscle/3.8.31
##gets my array id, which is needed for use below. this will be, i think, a number like 1,2,3 etc
MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
echo $MY_ARRAY_ID

## given an index file formatted                                                                        
## <index> <filename>                                                                                   
## produce the filename for given index                                                                 
THE_INDEX="""+indexname+"""
THE_INPUT_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
echo $THE_INPUT_FILE
ENDING=_Muscle.fasta
echo $THE_INPUT_FILE$ENDING

muscle -in $THE_INPUT_FILE -out $THE_INPUT_FILE$ENDING

exit"""
    with open(scriptfile, "w") as script:
        script.write(a)
    return scriptfile

def Overall_Gene_Trees(list_of_fasta_to_align, projectname, list_of_big_clades):
    print("Beginning to make all gene_trees")
    origlen = len(list_of_fasta_to_align)
    already_done_list = []
    for item in list_of_fasta_to_align:
        outputfile_exists = os.path.isfile("RAxML_bestTree."+item[:-6])
        if outputfile_exists is True:
            already_done_list.append(item)
        else:
            print(os.getcwd())
            print("looking for : "+projectname+"/Gene_Trees/trees/RAxML_bestTree."+item[:-6])
            outputfile_exists = os.path.isfile(projectname+"/Gene_Trees/trees/RAxML_bestTree."+item[:-6])
            #this didnt fix it
            if outputfile_exists is True:
                print("found")
                already_done_list.append(item)
    for remov in already_done_list:
        list_of_fasta_to_align.remove(remov)
    #now here is a list of fastas that have not been completed. many seem like they will error too-small. fix this?
    newlen = len(list_of_fasta_to_align)
    if already_done_list != []:
        print("done:")
        print (already_done_list)
    print("of an original: "+str(origlen)+" files to align, finished trees were found for "+str(len(already_done_list))+". We will now align "+str(newlen)+" files.")
    #should have been fixed by changing 4 to 5 as cutoff below value
    if list_of_fasta_to_align == []:
        rax_list = []
        muscle_list = []
    else:
        muscle_list = muscle_align_on_cluster(list_of_fasta_to_align, projectname)
        if muscle_list == "NA":
            muscle_list = []
    #muscle_list = []
    #rax_list = []
    ##TEMPORART
        mus_parsed_list = []
        for muscle_n in muscle_list:
            if muscle_n == "NA":
                pass
            else:
                mus_parsed_list.append(muscle_n)
        rax_list = raxml_run_on_cluster(mus_parsed_list, projectname)

        #rax list - contains only bipartitions or also bestTree ??
        #here rax list is mus_parsed_list - only new stuff
    #add back in the ones you already found
    #print(rax_list)
    for item in already_done_list:
        #item is the fasta file. Cyano0000_gene.fasta for example
        #print("added to rax_list: RAxML_bestTree."+item[:-6])
        outputfile_exists = os.path.isfile("RAxML_bestTree." + item[:-6])
        if outputfile_exists is True:
            #print("output exists")
            rax_list.append("RAxML_bestTree."+item[:-6])
        else:
            outputfile_exists = os.path.isfile(projectname + "/Gene_Trees/trees/RAxML_bestTree." + item[:-6])
            # this didnt fix it
            if outputfile_exists is True:
                
                rax_list.append(projectname + "/Gene_Trees/trees/RAxML_bestTree." + item[:-6])


    #add appropriate raxml filename to each subtrees instance. they should be called.... originailname+_Muscle.fasta+??

        #then, send everything over the MakeSpeciesTree.py to generate the species trees for everything.
        #-0still need determination of 
        
    #correlate gene trees
    raxset = 0
    total = len(rax_list)
    print("raxlist: "+str(total))
    for item in rax_list: 
        if "bipartitions" in item:
            rax_list.remove(item)
    total = len(rax_list)
    print("looking to match: "+str(total))
    #rax list should be quite long and include new and old stuff.
    for tree in rax_list:
        #print(tree)
      
      # it should be RAxML_bestTree.prefix_gene or bipartitions.
        #rax_list has all of the names, regardless of if the file actually exists or not.
        #print(tree)
        #RAxML_bestTree.pre_fix_gene
        rblah, rest = tree.split(".")
        if "_gene" in rest:
        #pre_fix_gene
            listunder = rest.split("_")
            #[pre,fix,gene]
            listunder.pop(-1)
            #[pre,fix]
            pref = "_".join(listunder)
            #pre_fix
        else:
            pref = rest
        match = "no"
        for st_o in list_of_big_clades:
            #print (st_o.ret_prefix())
            if st_o.ret_prefix() == pref:
                st_o.set_gene_tree_name(tree)
                #print("set to gene_tree_name: "+tree)
                raxset += 1
                match = "yes"
        if match == "no":
            print("no gene tree found for the tree: "+pref)
    #wtf why is raxset and total both ZERO?!?!??!
    if raxset == total:
        print("All RAXML trees are complete! Total new: "+str(raxset))
        incomplete = []

    else:
        print("tried to generate "+str(newlen)+" NEW gene raxml trees. "+str(raxset)+" gene trees were successfully matched with the correct subtree object!")
        print("we had "+str(len(already_done_list))+" already existing gene trees.")
        incomplete = []
        for set_tree in list_of_big_clades:
            if set_tree.gene_tree_name == "":
                incomplete.append(set_tree)
                print("no gene tree found for: "+set_tree.prefix)
    boots_list = Get_Attach_Boots(list_of_big_clades, projectname, incomplete)
    print(rax_list)
    print("genetrees raxlist ^")
    print(boots_list)
    print("genetrees bootslist")
    return rax_list, muscle_list, list_of_fasta_to_align, incomplete, boots_list


def Get_Attach_Boots(list_of_big_clades, projectname, incomplete):
    print("attaching boots")
    boots_list = []
    for item in list_of_big_clades:
        if item in incomplete:
            continue
        boot_name = "RAxML_bootstrap."+item.prefix+"_gene"
        #fix for the final round -- which doesn't have _gene associated! for whatever reason!
        if "_SS" in boot_name:
            boot_name = boot_name[:-5]
        #print("looking for: "+boot_name)
        #exists?
        a = os.path.isfile(boot_name)
        #if a is true, the boot file is already in the current folder, and should be properly moved (now or later?). 
        if a is True:
            os.system("mv "+boot_name+" "+projectname+"/Gene_Trees/boots/"+boot_name)
        #otherwise:
        if a is False:
            #see if the boot file is already in the /gene/boots final location:
            inplace = "./"+projectname+"/Gene_Trees/boots/"+boot_name
            b = os.path.isfile(inplace)
            #if so, leave it there!
            #otherwise, we need to fetch it from the cluster!
            if b is False:
                os.system("scp "+clus_head+"Gene_Trees/"+projectname+"/"+boot_name+" "+projectname+"/Gene_Trees/boots/"+boot_name)
                #scp the boot home
        item.boot_gene_file = boot_name
        
        #scp the boot home
        boots_list.append(boot_name)
        #add the boot name to a boots list
    return boots_list

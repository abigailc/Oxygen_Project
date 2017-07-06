# #!/usr/bin/python

# last edit abigailc@Actaeon on jan 27 2017

#pulling the taxonomy functions out of makespeciestree because I need to make them faster...
#insects is running for literally >20 hours.


names_file = "/Users/abigailc/Documents/Taxonomy_Stuff/taxdump/names.dmp"
nodes_file = "/Users/abigailc/Documents/Taxonomy_Stuff/taxdump/nodes.dmp"



######### PERSONAL_SETTINGS #########
ssh_inst = "ssh -l abigailc -i ~/.ssh/id_rsa eofe4.mit.edu"
clus_head = "abigailc@eofe4.mit.edu:/home/abigailc/"
Path_Blast = "/Users/abigailc/blast/"

import os
import re
import time
import sys
from Oxy_SS_modules.Classes_DTL_Detector import Fasta

#BASIC OPERATIONS
def Str_To_Taxid(string, names_file):
    #init done
    #turns a string to its taxon id NCBI
    #this is easier than expected. just open names.dmp and find the first hit. format:
    found = False
    #print("strtotaxid")
    #print(string+" str to taxid")
    string = string.replace("_", " ")
    #print(string)
    with open (names_file) as names:
        for line in names:
            
            if "\t"+string+"\t" in line:
                #print("got:"+line)
                taxid_int = re.sub ("(\d*)(\t\|\t)("+string+")(\t)(.*)", "\\1", line)
                found = True
                break
    if found is False:
        print("Error finding string: "+string+" in file: "+names_file)
        taxid_int = "NA"
    return taxid_int

def Taxid_To_Children(taxid, nodes_file):

    #goes one level deeper. finds all taxids that list the given taxid as "parent", returns as a list
    childlist = []
    child_rank_list = []
    with open (nodes_file) as nodes:
        for line in nodes:
            if "\t"+taxid+"\t" in line:
                #print("gotcha")
                #print(line)
                #the thing matches, do the re.sub.
                #includes the tab bc otherwise taxid 12 would match 12, 123, 12345355, etc
                baby_taxid_rank = re.sub("(\d*)(\t\|\t)("+taxid+")(\t\|\t)([a-z]*)(.*)", "\\1~\\5", line)
                if "\t" in baby_taxid_rank:
                    #this happens if the re.sub does not occur - eg if \ttaxid\t occured somewhere in the line other than where it should've. 
                    pass
                else:
                    baby_taxid, baby_rank = baby_taxid_rank.split("~")
                    #add to list of bbys
                    baby_taxid = baby_taxid.strip()
                    baby_rank = baby_rank.strip()
                    childlist.append(baby_taxid)
                    child_rank_list.append((baby_taxid, baby_rank))
    return child_rank_list

def Get_Taxid_Rank(taxid, nodes_file):
    taxid = taxid.strip()
    ranklist = []
    len_tax = len(taxid)
    len_tax_t = len_tax+1
    #given taxid = 100, len_tax = 3, len_tax_t = 5
    with open (nodes_file) as nodes:
        for line in nodes:
            #print(line[:len_tax_t])
            #print(taxid+"\t")
            if line[:len_tax_t] == taxid+"\t":
                #the thing matches, do the re.sub.
                #includes the tab bc otherwise taxid 12 would match 12, 123, 12345355, etc
                apparent_rank = re.sub("("+taxid+")(\t\|\t)(\d*)(\t\|\t)([a-z]*)(.*)", "\\5", line)
                apparent_rank = apparent_rank.strip()
                if "\t" in apparent_rank:
                    pass
                else:
                    return apparent_rank
        return "NA"
    #returns the rank (eg, "order" of a taxid")

def One_Rank_Lower(rank):
    print("looking one level lower than"+rank)
    if rank == "species":
        print("is species!")
        return "NA"
    ordered_str = "superkingdom kingdom phylum class order family genus species"
    ordered_list = ordered_str.split()
    if rank in ordered_list:
        pass
    elif rank == "NA":
        return "NA"
    else:
        print(rank+" is weird")
        return "NA"
    current = ordered_list.index(rank)
    lowindex = current + 1
    one_lower = ordered_list[lowindex]
    return one_lower
    
    #given phylum, returns class. given class, return order. etc.

# rank = "class"
# string = "cyanobacteria"
# taxid = "12345"

def Return_Parent(taxid, nodes_file):
    #eg for a given rank taxid, find it's up-one-level (not rank) taxid, and return it.
    len_tax = len(taxid.strip())
    len_tax_t = len_tax+1
    #given taxid = 100, len_tax = 3, len_tax_t = 5
    #print("searching for one level above taxid:"+str(taxid))
    #print("tiud: "+taxid)
    with open (nodes_file) as nodes:
        for line in nodes:
            #print(taxid.strip()+"\t")
            #print(line[:len_tax_t])
            if line[:len_tax_t] == taxid.strip()+"\t":
                
               # print("got: "+line)
                #the thing matches, do the re.sub.
                #includes the tab bc otherwise taxid 12 would match 12, 123, 12345355, etc
                parent_taxid = re.sub("("+taxid.strip()+")(\t\|\t)(\d*)(\t\|\t)([a-z]*)(.*)", "\\3", line)
                #print(parent_taxid)
                if "\t" in parent_taxid:
                    pass
                else:
                    return parent_taxid
    print("error finding parent taxa")
    return("NA")


#COMPLEX OPERATIONS

def Ret_A_Valid_Species_Below_LESS_EFFICIENTLY(taxid, nodes_file, names_file, acc_list):
    children = []
    list_ch_remove = []
    child_list_a = []
    #this is a list of children : [ [child_taxid, child_rank], [child2_taxid, child2_rank] ] 
    child_list_atup = Taxid_To_Children(taxid, nodes_file)
    #this is a list of children TAXIDS ONLY

    #print("initial pass")
    #print(child_list_atup)
    #print(child_list_a)
    done = False
    saved_top_level = []
    #we're going to do one at a time, so save all, and load them one-by-one.
    for itema in child_list_atup:
        saved_top_level.append(itema)
    maxi = len(saved_top_level)
   # print("maxi: "+str(maxi))
    atup = saved_top_level[0]
    saved_top_level.remove(atup)
    child_list_atup = [atup]
    for item in child_list_atup:
        child_list_a.append(item[0])
    i = 1
    #also lets implement a saved second level... for  further spe.
    while done is False:
        for item in child_list_atup:
            if item[1] == "species":
                #add the taxid to the list of species_level_children
                children.append(item[0])
                sis_spec_name = Taxid_To_Name(item[0], names_file)
                if sis_spec_name[0].islower() is False:
                    in_blast = Check_Spec_Name_Acceptable_List(sis_spec_name, acc_list)
                    if in_blast is True:
                        return sis_spec_name
                list_ch_remove.append(item)
        #remove taxids that were saved at the species level
        #print(list_ch_remove)
        for rem in list_ch_remove:
            child_list_atup.remove(rem)
            child_list_a.remove(rem[0])
        #if all tips have terminated at the species level: you are done.
        if child_list_a == []:
            if i == maxi:
                #print("found none")
                return "NA"
                done = True
            else:
                i += 1
                #print(i)
                list_ch_remove = []
                
                atup = saved_top_level[0]
                #print(atup)
                saved_top_level.remove(atup)
                child_list_atup = [atup]
                #print(child_list_atup)
                for item in child_list_atup:
                    child_list_a.append(item[0])
                continue
        list_ch_remove = []
        child_list_b = []
        child_list_c = []
        for parent in child_list_a:
            child_list_btup = Taxid_To_Children(parent, nodes_file)
            for item in child_list_btup:
                child_list_b.append(item[0])
            if child_list_btup == []:
                pass
            else:
                for bitem in child_list_btup:
                    child_list_c.append(bitem)
        child_list_atup = child_list_c
        #print("New parent list:")
        #print(child_list_atup)
        child_list_a = []
        for itup in child_list_atup:
            child_list_a.append(itup[0])
        #print(child_list_a)
    #children is a list of all species-level TAXIDS that belong to the given group. 
    return "NA"

#WHY ARE THERE TWO OF THESE???????
def Ret_A_Valid_Species_Below(taxid, nodes_file, names_file, acc_list):
    masterlist = []
    #this is a list of children : [ [child_taxid, child_rank], [child2_taxid, child2_rank] ] 
    complete = False
    masterlist.append([(taxid, "starter")])
    while complete is False:
        #print(masterlist)
        if masterlist == []:
            return("NA")
        #now lookat is the last member of the last list in masterlist.
        now_list = masterlist[-1]
        if now_list == []:
            while [] in masterlist:   
                masterlist.remove([])
            if masterlist == []:
                return("NA")
            now_list = masterlist[-1]
        #lookat first member of that list.
        now_tup = now_list[0]
        now_taxid, now_rank = now_tup[0], now_tup[1]
        #see if its a species
        if now_rank == "species":
            now_list.remove(now_tup)
            now_name = Taxid_To_Name(now_taxid, names_file)
            if now_name[0].islower() is False:
                in_blast = Check_Spec_Name_Acceptable_List(now_name,acc_list)
                if in_blast is True:
                    #now_name is a species_name
                    return now_name
            #check if now_tup is valid. if so, return.
        else:
            now_list.remove(now_tup)
            #generate a new list - of the descendents of this one.
            newlist = Taxid_To_Children(now_taxid, nodes_file)
            #print(newlist)
            if newlist == "NA":
                pass
            else:
            #add it to masterlist.
                masterlist.append(newlist)
    return("Uh, what?")

def Ret_All_Species_Below_Less_Efficiently(taxid, nodes_file):
    children = []
    list_ch_remove = []
    child_list_a = []
    #this is a list of children : [ [child_taxid, child_rank], [child2_taxid, child2_rank] ] 
    child_list_atup = Taxid_To_Children(taxid, nodes_file)
    #this is a list of children TAXIDS ONLY
    for item in child_list_atup:
        child_list_a.append(item[0])
    #print("initial pass")
    #print(child_list_atup)
    #print(child_list_a)
    done = False
    while done is False:
        for item in child_list_atup:
            if item[1] == "species":
                #add the taxid to the list of species_level_children
                children.append(item[0])
                list_ch_remove.append(item)
        #remove taxids that were saved at the species level
        for rem in list_ch_remove:
            child_list_atup.remove(rem)
            child_list_a.remove(rem[0])
        #if all tips have terminated at the species level: you are done.
        if child_list_a == []:
            done = True
        list_ch_remove = []
        child_list_b = []
        child_list_c = []
        #for remaining non-species level taxids in lista:
        #   -get their children (listb)
        #   -add their children to a persistant list(listc)
        #   -then set lista(the list to check and remove species-level-entries) to be  ==  listc.
        for parent in child_list_a:
            child_list_btup = Taxid_To_Children(parent, nodes_file)
            for item in child_list_btup:
                child_list_b.append(item[0])
            if child_list_btup == []:
                pass
            else:
                for bitem in child_list_btup:
                    child_list_c.append(bitem)
        child_list_atup = child_list_c
        #print("New parent list:")
        #print(child_list_atup)
        child_list_a = []
        for itup in child_list_atup:
            child_list_a.append(itup[0])
        #print(child_list_a)
    #children is a list of all species-level TAXIDS that belong to the given group. 
    return children


def Ret_All_Groups_One_Rank_Below(taxid, nodes_file):
    taxid = taxid.strip()
    print("looking for taxid:"+str(taxid))
    rank = Get_Taxid_Rank(taxid, nodes_file)
    print(rank)
    #raise SystemExit
    target_rank = One_Rank_Lower(rank)
    if target_rank == "NA":
        return("NA")
    removal_ranks = "superkingdom kingdom phylum class order family genus species"
    garbage, remove_string = removal_ranks.split(target_rank)
    remove_rank_list = remove_string.split()
    children = []
    list_ch_remove = []
    #print(remove_rank_list)
    #this is a list of children : [ [child_taxid, child_rank], [child2_taxid, child2_rank] ] 
    child_list_a = Taxid_To_Children(taxid, nodes_file)
    done = False
    while done is False:
        for item in child_list_a:
            if item[1] == target_rank:
                #add the taxid to the list of species_level_children
                children.append(item[0])
                list_ch_remove.append(item)
            if item[1] in remove_rank_list:
                list_ch_remove.append(item)
        #remove taxids that were saved at the species level
        for rem in list_ch_remove:
            child_list_a.remove(rem)
        #if all tips have terminated at the target species level: you are done.
        if child_list_a == []:
            done = True
        list_ch_remove = []
        child_list_b = []
        child_list_c = []
        #for remaining non-species level taxids in lista:
        #   -get their children (listb)
        #   -add their children to a persistant list(listc)
        #   -then set lista(the list to check and remove species-level-entries) to be  ==  listc.
        for parent in child_list_a:
            child_list_b = Taxid_To_Children(parent[0], nodes_file)
            if child_list_b == []:
                pass
            else:
                for bitem in child_list_b:
                    child_list_c.append(bitem)
        child_list_a = child_list_c
        #print(child_list_a)
    #children is a list of all ONE-RANK-BELOW level TAXIDS that belong to the given group. 
    return children
    #runs until all children are found of one rank below. eg (CLASS -> [order1, order 2, order3, order 4)
    #for checking loss candidates, i will want to 1) run this 2) run a species_level_children generation for each member of the output list. 3) chose one member of each of those output lists to go in the species tree. hopefully checking that we have data for the chosen species.

    

def Ret_Sister_Same_Rank(string, nodes_file, names_file):
    #from str rank - get current taxid, go up one level, then get all descendents in a list, remove the current taxid, and return the resulting sister list
    print(string)
    interest_taxid = Str_To_Taxid(string, names_file)
    print(interest_taxid)
    up_taxid = Return_Parent(interest_taxid, nodes_file)

    up_taxid = up_taxid.strip()
    interest_taxid = interest_taxid.strip()
    sis_self_tuples = Taxid_To_Children(up_taxid, nodes_file)
    sister_and_self = []
    for tup in sis_self_tuples:
        sister_and_self.append(tup[0])
        #sis_and_self is a list of TAXIDS ONLY
    print(sister_and_self)
    print(interest_taxid)
    sister_and_self.remove(interest_taxid)
    sisterlist = sister_and_self
    print(sisterlist)
    return sisterlist
#sisterlist will be a list of taxids for the sister clades to the current thing. by level, not by rank.
#todo = implement something to redo if sisterlist is empty.

def Taxid_To_Name(taxid, names_file):
    #this needs to be the  backwards version of Str to Taxid.
    found = False
    taxid = taxid.strip()
    len_tax = len(taxid)
    len_tax_t = len_tax+1
    with open (names_file) as names:
        for line in names:
            if line[:len_tax_t] == taxid+"\t":
               # print("got here")
                name_wanted = re.sub ("(\d*)(\t\|\t)([^\t]*)(\t\|\t)(.*)(\t\|\t)(scientific name)(.*)", "\\3", line)
                if "\t" in name_wanted:
                    pass
                else:
                    found = True
                    break
    if found is False:
        print("Error finding name for: "+taxid+" in file: "+names_file)
        name_wanted = "NA"
    if found is True:
        #print(name_wanted)
        name_wanted = name_wanted.strip()
    return name_wanted
                    
def Choose_One_OG_Seq(string, species_list, names_file, acc_list, nodes_file):
    print("one og sequence choser initiating")
    if "_" in string:
        string = string.replace("_", " ")
    sislist = Ret_Sister_Same_Rank(string, nodes_file, names_file)
    print("Sisterlist")
    print(sislist)
    if sislist == []:
        go = True
    else:
        go = False    
    my_taxid = Str_To_Taxid(string, names_file)
    while go is True:
        parent_of_me_taxid = Return_Parent(my_taxid, nodes_file)
        parent_of_me = Taxid_To_Name(parent_of_me_taxid, names_file)
        sislist = Ret_Sister_Same_Rank(parent_of_me, nodes_file, names_file)
        my_taxid = parent_of_me_taxid
        if sislist == []:
            pass
        else:
            go = False
    for item in sislist:
        #spec_sis_list = Ret_All_Species_Below(item, nodes_file)
        test = Ret_A_Valid_Species_Below(item, nodes_file, names_file, acc_list)
        if test == "NA":
            pass
        else:
            print(test)
            return test
    #if test == "None":
    #    return "None"
    #if nothing in the first level sister list is a valid hit, keep moving up the tree until you get one.

    while test == "NA":
        sislist = []
        go = True
        if my_taxid == 1:
            break
        while go is True:
            parent_of_me_taxid = Return_Parent(my_taxid, nodes_file)
            parent_of_me = Taxid_To_Name(parent_of_me_taxid, names_file)
            sislist = Ret_Sister_Same_Rank(parent_of_me, nodes_file, names_file)
            my_taxid = parent_of_me_taxid
            if sislist == []:
                pass
            else:
                go = False
        for item in sislist:
            test = Ret_A_Valid_Species_Below(item, nodes_file, names_file, acc_list)
            if test != "NA":
                pass
            else:
                return test
    return test
        
        #print (spec_sis_list)
        #for sis_spec_taxid in spec_sis_list:
         #   sis_spec_name = Taxid_To_Name(sis_spec_taxid, names_file)
          #  in_blast = Check_Spec_Name_Blast_File(sis_spec_name, blast_file)
           # if in_blast is True:
            #    print("Outgroup sequence chosen:"+sis_spec_name)
             #   return sis_spec_name

            

    #double break so we only keep ONE sequence.
        #go all the way down the first one until you get a species-level entry.
        #check if the species-level entry is found in your .blast file (if that is where we are implementing this??? )
        #if not, continue... check each species-level thing you find.
        #this would then need to be included in make_species_trees... and only called if the request is sent directly from Parser_blah_master.
def Check_If_We_Have_A_Rep_Already(species_list, tid_list, rank):
    print("Checking for reps... target rank is: "+rank)
    list_of_correct_rank = []
    found = []
    removal_ranks = "superkingdom kingdom phylum class order family genus species"
    remove_string, garbage = removal_ranks.split(rank)
    remove_rank_list = remove_string.split()
    for species in species_list:
        nid = Str_To_Taxid(species, names_file)
        #go up the ladder
        go = True
        while go is True:
            #get parent taxid
            rp = Return_Parent(nid, nodes_file)
            #if its 1, we're done.
            if rp == "NA":
                list_of_correct_rank.append(rp)
                go = False
            if rp.strip() == 1:
                rp = "NA"
                list_of_correct_rank.append(rp)
                go = False
            #get rank for that new taxid
            par_rank = Get_Taxid_Rank(rp, nodes_file)
            #if it's what we want it to be, add to list.
            if par_rank == rank:
                rp = rp.strip()
                list_of_correct_rank.append(rp)
                go = False
            #if its a step too high, terminate - we went too far somehow
            elif par_rank in remove_rank_list:
                rp = "NA"
                list_of_correct_rank.append(rp)
                go = False
            #else, go up another level and test that one!
            else:
                nid = rp
    print(tid_list)
    print(list_of_correct_rank)
    for item in tid_list:
        if item in list_of_correct_rank:
            a = tid_list.index(item)
            found.append(tid_list[a])
    return found

#@blast_file should actually be a list of raw_blast_FASTA objects
def Choose_Loss_Candidates(string, species_list, names_file, acc_list, nodes_file):
    print("loss search initiating")
    if "_" in string:
        print(string)
        string = string.replace("_", " ")
        print(string)
    taxid = Str_To_Taxid(string, names_file)
     #for checking loss candidates, i will want to 1) run this 2) run a species_level_children generation for each member of the output list. 3) chose one member of each of those output lists to go in the species tree. hopefully checking that we have data for the chosen species.
    sub_taxids = Ret_All_Groups_One_Rank_Below(taxid, nodes_file)
    if sub_taxids == "NA":
        print("Error getting loss candidates for string:"+string)
        return([])
    subgroup_names = []
    for item in sub_taxids:
        subgroup_names.append(Taxid_To_Name(item, names_file))
    b = Get_Taxid_Rank(taxid, nodes_file)
    a = One_Rank_Lower(b) 
    found = Check_If_We_Have_A_Rep_Already(species_list, sub_taxids, a)
    print("Representatives already exist for:")
    found_names = []
    for foundtid in found:
        foundtid = foundtid.strip()
        index1 = sub_taxids.index(foundtid)
        found_names.append(subgroup_names.pop(index1))
        del sub_taxids[index1]
    print(found_names)
    print("Looking for one representative from each of the following:")
    print(subgroup_names)
    loss_list = []
    ite = 0
    # #first check if it is in the output loss list.
    # for item in sub_taxids:
    #     with open(saved_loss_candidates) as saved:
    #         for line in saved:
    #             if item in line:
    #                 #newthing will be a species name.
    #                 newthing = re.sub("("item")(\t)(.*)", "\\3", line))
    #                 loss_list.append(newthing)
    #                 found2.append(item)
    #                 break
    #remove those found from file from the search list.
    # for item in found2:
    #     sub_taxids.pop(item)
    for item in sub_taxids:
        test = Ret_A_Valid_Species_Below(item, nodes_file, names_file, acc_list)
        #print(test)
        print(subgroup_names[ite]+" : "+test)
        ite+=1
        loss_list.append(test)
        continue
    print("Loss candidates will be added:")
    na = 0
    for item in loss_list:
        if item == "NA":
            na +=1
    while "NA" in loss_list: loss_list.remove("NA")
    
    print(loss_list)
    print("there were "+str(na)+" "+a+"s that no suitable loss candidate was found for.")
    return loss_list
    #either one per next-level-down
    #or one per next-rank-down
    
def Check_Spec_Name_Acceptable_List(ssp_name, acc_list):
    if ssp_name in acc_list:
        return True
    else:
        
        result = next((True for item in acc_list if ssp_name in item), False)
        if result is True:
            print("Err in match spec name - gen list: "+ ssp_name +" "+ item)
        return result

    
    
def Check_Spec_Name_Blast_File(ssp_name, blast_fasta_list):
    lf = (len(blast_fasta_list))
    half = lf/2
    yes = 0
    att = 0
    #print("Checking :"+ssp_name)
    ssp_name = ssp_name.replace(" ", "_")
    ssp_name = ssp_name.strip()
    for current_blast in blast_fasta_list:
        att += 1
        if att > 6:
            if yes < att/3:
                return False
        if ssp_name in current_blast.species_names:
            
            yes += 1
            continue
        else:
        
            #print(ssp_name)
            #print(current_blast.species_names[0])
            for spec in current_blast.species_names:
                if ssp_name in spec:
                    yes +=1
                    break
            continue
    #print(yes)
    #print(half)
    if yes > half:
        #print("validated: "+ssp_name)
        return True
    else:
        return False

def gen_acceptable_species_list(list_raw_gene_fastas, acc_name):
    #this is printing an empty file. why?
    names_list_acc = []
    numbers_list_acc = []
    
    for raw in list_raw_gene_fastas:
        #do they have species lists?
        raw.gen_species_lists()
        raw_sl = raw.species_names
        print(raw_sl[0])
        for rawsp in raw_sl:
            if rawsp in names_list_acc:
                ind = names_list_acc.index(rawsp)
                numbers_list_acc[ind] = numbers_list_acc[ind]+1
            else:
                names_list_acc.append(rawsp)
                numbers_list_acc.append(1)
    #the numbers list can specify a cut off that is necesary for the thing being acceptable
    #for now let's be consistant and use 1/2 of lsit of raw fastas?
    cutoff_num = (len(list_raw_gene_fastas)/2)
    print(cutoff_num)
    #this will be 15 currently. might be .5 sometimes.
    list_of_rem = []
    index = 0
    for n in numbers_list_acc:
        if n > cutoff_num:
            #means that we dont care if its a decimal or not. 1 will pass .5
            pass
        else:
            list_of_rem.append(names_list_acc[index])
            #add the index to be removed to a list. index into names and numbers should be identicle
        index +=1
    print(len(list_of_rem))
    list_of_rem.sort(reverse=True)
    for remove_me in list_of_rem:
        #uhhhhh i think we need to sort the numbers so removal of the largest number happens first so as to not fuck up list order.
        #sorting now. should be good.
        names_list_acc.remove(remove_me)
    a = write_acc_list(names_list_acc, acc_name)   
    return a
                
def write_acc_list(acc_list, acc_name):
    with open(acc_name, "w") as acc_list_file:
        for item in acc_list:
            acc_list_file.write(item+"\n")
    return acc_name

def write_spc_list(spc_list, spcname):
    with open(spcname, "w") as spc_list_file:
        for item in spc_list:
            #stripiing strain data from this version of the species_list such that it will 
            if "_" in item:
                dash_sep = item.split("_")
                item = dash_sep[0]+"_"+dash_sep[1]
            spc_list_file.write(item+"\n")
    return spcname

#parser stuff


def Run_OG_LOSS_ON_CLUSTER(script_name,all_files, all_result_files):
    #here acc list is the name of the acc_list_current_file
    #auto gen an sbatch script
    os.system(ssh_inst+" \'mkdir Taxonomy\'")
    sb_script = script_name
    #scp it over
    
    print(all_files)
    for item in all_files:
        os.system("scp "+item+" "+clus_head+"Taxonomy")
    #run it

    #edit the script on the cluster to deal with my mistakes

    os.system(ssh_inst+" 'cd ~/Taxonomy; sbatch "+sb_script+"'")
    #scp it back and verify
    direct = os.getcwd()
    exists = False
    #now it should exist locally
    movehome = []
    finished = "start"
    #bring home the d
    for i in all_result_files:
        movehome.append(i)
    while finished is not True:
        for filename in movehome:
            os.system("scp "+clus_head+"Taxonomy/"+filename+" "+direct)
        for item in all_result_files:
            #see if it got moved home.
            exists = os.path.isfile(item)
            if exists is True:
                if item in movehome:
                    movehome.remove(item)
                finished = "yes"
            else:
                finished = False
                print("Tax not done yet. could not locate : "+item+"checking again in 5 minutes")
                break
        if finished == "yes":
            print("Should be done!")
            finished = True
        else:
            #wait ten minutes and then try again.
            time.sleep(600)
            finished = "yes"
#TEMPORARILY REMOVED result file deletion from the cluster to make testing progress faster.
    #for item in all_result_files:
    #    os.system(ssh_inst+" 'cd ~/Taxonomy; rm "+item+"'")
    #for item in all_files:
    #    os.system(ssh_inst+" 'cd ~/Taxonomy; rm "+item+"'")
    print("Taxonomy parsing complete")
    #remove the script and the og loss file from cluster




            
def Get_OG_LOSS_DATA(list_of_clades, projectname):
    #the acceptable list should be a list of taxa that are present in at least 50% (?) of the blast hit files for the genes given.

    #get all gene-query-files to look at
    list_catfiles = []
    list_of_lists_of_raw_blast_files = []
    for item in list_of_clades:
        catfile = item.cat_file
        list_of_raw_blast_files = item.blast_raw
        if catfile in list_catfiles:
            pass
        else:
            list_catfiles.append(catfile)
            list_of_lists_of_raw_blast_files.append(list_of_raw_blast_files)
    cat_acc_dict = {}

    #for each, create an acceptable list output name
    for i in range(len(list_catfiles)):
        item = list_catfiles[i]
        list_raws = list_of_lists_of_raw_blast_files[i]
        gsflist = item.split(".")
        gsf_a = gsflist[0]
        gsf_b = gsf_a.split("/")[-1]
        acc_file = gsf_b+"_Acc_List.txt"
        #print("Looking for loss-candidates and a rooting sequence to add....")
        acc_exists = os.path.isfile(acc_file)
        if acc_exists is True:
            pass
        #if not already done, actually make the output acceptable list.
        else:
            print("....initializing all_acceptables from gene_seq_query file: "+gsf_b+". this should only happen once...")
            #generate it
            #should be passing in A LIST OF ALL THE BLAST_FILES ASSOCIATED WITH THE GENE. eg the things in Raw_Blasts that were consulted.
            #are these stored in each subtree? should pass a list of fasta objects.
            #ist_raw_objects = []
            #rint(list_raws)
            #or raw in list_raws:
            #   print(raw.name)

            acc_file = gen_acceptable_species_list(list_raws, acc_file)
            #this is returning "NONE" which is super not okay.
        cat_acc_dict[item] = acc_file
    
    list_of_species_files = Gen_Species_File(list_of_clades, projectname)

    #check if we already ran the taxonomy and have data downloaded. (this is mostly for while fixing errors; i keep getting stuck at this point & ity is a waste of time to re-run the taxonomy parser.
    list_to_tax_clades = []
    for item in list_of_clades:
        exists_result = os.path.isfile(item.result)
        if exists_result is False:
            list_to_tax_clades.append(item)
    #sets species_file and result to each subtree.
    corr_file_name, results_list = Generate_Cat_File_OGLOSS(list_to_tax_clades, cat_acc_dict, projectname)
    #makes the correlation file.
    #for each clade, generate a species_list, result name, acc_file_name, string_name and print them all to a corr.file
    n = len(list_to_tax_clades)
    #gen the script
    script_name = projectname+"_OGLScript.sh"
    scriptfile = Generate_Script_File_OGLOSS(n, corr_file_name, script_name)
    all_files = []
    for item in cat_acc_dict.values():
        all_files.append(item)
    for item in list_of_species_files:
        all_files.append(item)
    all_files.append(scriptfile)
    all_files.append(corr_file_name)
    if len(results_list) is 0:
        pass
    else:
        Run_OG_LOSS_ON_CLUSTER(scriptfile,all_files, results_list)
    
    #run the script

    #add loss_species, root_species to each subtree as a value and also add them to the species_list going forward.
    for item in list_of_clades:
        results_file = item.result
        loss_species = []
        print(item.string_name)
            #open the file and get loss and species results.
        with open(results_file) as res:
        #    print("opened")
            a=0
            for line in res:
                #get loss results
                if a == 0:
                    loss_species = line.strip()
                    loss_species = loss_species.split("~")
                    print("loss candidates")
                    if "" in loss_species:
                        loss_species.remove ("")
                    if "\n" in loss_species:
                        loss_species.remove("\n")
                    item.loss_species_list = loss_species 
                    print(loss_species)
                #get root results
                if a == 1:
                    root_species = line.strip()
                    item.root_species = root_species
                    print("root: "+root_species)
                #get how long it took
                if a == 2:
                    print("time:")
                    print(line)
                a += 1
        #if no loss, do nothing

        item.species_list_plus_og_loss = []
        for thing in item.species_list_original:
            item.species_list_plus_og_loss.append(thing)
        if loss_species == []:
            pass
        #else, add them to the species list, and also track them(?)
        else:
            for ls in loss_species:
                item.species_list_plus_og_loss.append(ls)
                
                
        if root_species == "":
            pass
        
        else:
            item.species_list_plus_og_loss.append(root_species)
    return results_list
#        os.system("rm "+results_file)

        #done

def Generate_Cat_File_OGLOSS(list_of_clades, cat_acc_dict, projectname):
    
    corr_file_name = "Corr_"+projectname+".txt"
    results_list = []
    with open(corr_file_name, "w") as corr:
        for n in range(len(list_of_clades)):
            corr.write(str(n+1)+" "+list_of_clades[n].species_file+" "+list_of_clades[n].string_name+" "+cat_acc_dict[list_of_clades[n].cat_file]+" "+list_of_clades[n].result+"\n")
            results_list.append(list_of_clades[n].result)
    return corr_file_name, results_list

def Generate_Script_File_OGLOSS(n, indexname, scriptname):
    n = str(n)
    a =  """#!/bin/bash                                                                                             
#SBATCH -p sched_mit_g4nier                                                                             
#SBATCH -t 2-00:00:00                                                                                   
#SBATCH -J Tax
  
#SBATCH --array=1-"""+n+"""                                                                                        

. /etc/profile.d/modules.sh
module add engaging/openmpi/1.8.8

MY_ARRAY_ID=$SLURM_ARRAY_TASK_ID
THE_INDEX="""+indexname+"""
SPECIES_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $2}' )
STRING_NAME=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $3}' )
ACC_FILE=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $4}' )
RESULT=$( cat $THE_INDEX | grep "^$MY_ARRAY_ID " | awk '{print $5}' )

echo $SPECIES_FILE
echo $STRING_NAME
echo $ACC_FILE

mpirun python Online_Taxon_Parse.py -s $SPECIES_FILE -g $STRING_NAME -b $ACC_FILE -n $RESULT

exit"""
    with open(scriptname, "w") as script:
        script.write(a)
    return scriptname

def Gen_Species_File(list_of_clades, projectname):
    list_sp_files = []
    for item in list_of_clades:
        species_list = item.species_list_original
        species_file_name = item.prefix+"_Species_List.txt"
        species_list2 = []
        for sl2 in species_list:
            sl2 = sl2.strip("\"")
            species_list2.append(sl2)
        spc_file = write_spc_list(species_list2, species_file_name)
        item.species_file = species_file_name
        list_sp_files.append(species_file_name)
        item.result = item.prefix+"_OGL_Result.txt"
    return list_sp_files


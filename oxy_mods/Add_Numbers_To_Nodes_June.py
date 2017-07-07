#!/usr/bin/python

#abigailc@Actaeon Mar 20 2017
#update abigailc@Actaeon April19 -- working

#currently not run by the overall program. 


#idea: open a newick file that has no node_labels (eg: a  best_tree)

#you have a list of numbers you want to associate with each node.
#node is defined by it's downstream tips

#eg tree (A:0.05,(B:0.09,((C:0.01,D:0.02):0.38,((D:0.02,(E:0.0052,F:0.001):0.010):0.10,G:0.19):0.24)100:0.28):0.14,H:0.05);
#has a node called [C,D]
#if we want to add the number 8 to [C,D]
#output should be:
# (A:0.05,(B:0.09,((C:0.01,D:0.02)8:0.38,((D:0.02,(E:0.0052,F:0.001):0.010):0.10,G:0.19):0.24)100:0.28):0.14,H:0.05);
#adding 100s everywhere gives:
#(A:0.05,(B:0.09,((C:0.01,D:0.02)100:0.38,((D:0.02,(E:0.0052,F:0.001)100:0.010)100:0.10,G:0.19)100:0.24)100:0.28)100:0.14,H:0.05);

#this will be used to add transfer recipient frequency to clades
#is the given gene tree in format >Bac|Cyano|Nosto|Nostoc_punctiforme_PCC_73102|##### or
#                                 >Nostoc_punctiforme_PCC_73102
#                                 >NostocXX0XXpunctiformeXX0XXPCCXX0XX73102
#generate a new test file

import os
import re
import glob
import operator

#todo:
#move besttrees to better location
#figure out how to write fileinformation to a nexus file as appropriate. make this a definition.

def Get_All_Files_To_Run_On(directory, projectname, reroot = True):
    biparts_files = []
    bestree_files = []
    list_of_prefixes = []
    rootstatus = ""
    os.chdir(directory)
    for name in glob.glob(projectname+"/Ranger_DTL/*biparts_table.txt"):
        biparts_files.append(name)
    print(biparts_files)
    #old version: 
    if reroot is False:
        for name in glob.glob(projectname+"/Gene_Trees/trees/RAxML_bestTree.*"):
            bestree_files.append(name)
    #with optimal rooting:
    else:
        rootstatus = "_rerooted_"
        for name in glob.glob(projectname+"/Gene_Trees/trees/Best_Rooted_Renamed_Gene_*"):
        #this is currently broken due to the new_named BIPARTS files (i think)
            bestree_files.append(name)

    if bestree_files == []:
        rootstatus = "_rerooted_"
        bestree_files = Populate_Bestree(directory, projectname)


###################################

    print(bestree_files)
    print(biparts_files)
    print("above")
    biparts_files.sort()
    bestree_files.sort()
    a = len(projectname)
    for item in biparts_files:
        iteml = item.split("/")
        item = iteml[-1]
        if "\\" in item:
            item1 = item.split("\\")
            item = item1[-1]
        if "newBS" in item:
            item = item.split("RangerOut")
            list_of_prefixes.append(item[0])
        #print(item)
        else:
            item = item[a:-17]
        #print(item)

            list_of_prefixes.append(item)
    #@raise SystemExit
    print(list_of_prefixes)
    new_bestree_files = []
    new_biparts_files = []
    rem_pref = []
    if len(bestree_files) == len(biparts_files):
        print("equal lists")
    else:
    
        print("uh oh unequal. why?")
        for item in list_of_prefixes:
            #print("matching : "+item)
            got = False
            get = False
            for thing in bestree_files:
                #print(thing)
                if item in thing:
                    got = True
                    new_bestree_files.append(thing)
                    break
            if got is False:
                print("couldnt match gene file for: "+item)
                rem_pref.append(item)
                continue
            for thing in biparts_files:
                if item in thing:
                    get = True
                    new_biparts_files.append(thing)
                    break
            if get is False:
                print("No species file for: "+item)
                new_bestree_files.remove(thing)
                rem_pref.append(item)
        for pre in rem_pref:
            list_of_prefixes.remove(pre)
        #print(new_besttree_files)
        bestree_files = new_bestree_files
        biparts_files = new_biparts_files
    print(len(bestree_files))
    print(len(biparts_files))
    print(len(list_of_prefixes))
    if len(bestree_files) != len(biparts_files):
        print(bestree_files)
        print(biparts_files)
        print("still an error")
        raise SystemExit
    print("made file lists containing "+str(len(biparts_files))+" files")
    for i in range(len(bestree_files)):
        Overall_Add_Numbers(bestree_files[i], biparts_files[i], projectname, projectname+list_of_prefixes[i]+rootstatus, directory)

#gene_tree should be a FILENAME, bipart_file should be a FILENAME, projectname should be a STRING
def Overall_Add_Numbers(gene_tree, biparts_file, projectname, prefix, directory):
    with open(gene_tree) as gene:
        for line in gene:
            if line != "":
                gene_tree_string = line
                break
    gene_tree_string = gene_tree_string.strip()
    #print(gene_tree_string)
    print("beginning to make dictionary of string to included tips for file: "+gene_tree)
    str_tips_dict, overall_tiplist = String_to_Tips(gene_tree_string)
    #print(str_tips_dict)
    #tips are of type Archaea|Eury|Halo|Blah_blah|gi|940578439857
   # print(str_tips_dict)
    #this dictionary holds: (nost1, (nost2, nost3)): [nost1,nost2,nost3]
    interest = "Transfers_By_Recipient_Clade"
    print("Opening biparts file: "+biparts_file)
    print("Making dictionaries of number of losses, speciations, duplications, transfers by clade")
    #each of these holds a list of tips : number of times that happened in the BSs
    tips_trans_recip_dict = Parse_biparts_file(biparts_file, interest, overall_tiplist)
    tips_trans_sp_dict = Parse_biparts_file(biparts_file, "Transfers_By_Species_Tree_Recipient", overall_tiplist)
    tips_trans_dict = Parse_biparts_file(biparts_file, "Transfers_By_Representative_Node", overall_tiplist)
    tips_loss_dict = Parse_biparts_file(biparts_file, "Losses", overall_tiplist)
    tips_spec_dict = Parse_biparts_file(biparts_file, "Speciation", overall_tiplist)
    tips_dups_dict = Parse_biparts_file(biparts_file, "Duplications", overall_tiplist)

    #each of these holds a string (subtree) : number of times that happened in the BSs
    #print("AHHHHHH need to convert tips from bi[parts file to tips from besttree format!!!!")
    print("making dictionary of clades in the given gene tree to events in the species tree")

    str_val_dict_trans_recip = MakeStrValDict(str_tips_dict, tips_trans_recip_dict)
    str_val_dict_trans = MakeStrValDict(str_tips_dict, tips_trans_dict)
    #print(str_val_dict_trans)
    str_val_dict_loss = MakeStrValDict(str_tips_dict, tips_loss_dict)
    str_val_dict_dups = MakeStrValDict(str_tips_dict, tips_dups_dict)
    str_val_dict_spec = MakeStrValDict(str_tips_dict, tips_spec_dict)
    str_val_dict_trans_sp = MakeStrValDict(str_tips_dict, tips_trans_sp_dict)
    #this holds each subtree string: number of times it happened in the BSs
    #eg (C,(D,E)): 20

    print("making a dictionary of the number of times a given clade-topology was found across bootstraps")
    topology_dict = topology_support_dict([str_val_dict_spec, str_val_dict_dups, str_val_dict_loss,str_val_dict_trans], overall_tiplist)

    #NOW, MAKE THE OUTPUT TREES


    #write a tree with just the transfer ##s
    print("writing the transfer tree")
    trans_gene_tree = AddValToTree(gene_tree_string, str_val_dict_trans_recip)
    os.system("mkdir "+projectname+"/AnnotatedTrees")
    with open(projectname+"/AnnotatedTrees/"+prefix+"_transfers.newick", "w") as new:
        new.write(trans_gene_tree)
    print("located at: "+projectname+"/AnnotatedTrees/"+prefix+"_transfers.newick")
    #write a tree with transfer number / overall topology found number
    # eg 21/30
    print("writing a topology support tree... this should be == bipartitions tree????")
    top_gene_tree = AddValToTree(gene_tree_string, topology_dict)
    with open(projectname+"/AnnotatedTrees/"+prefix+"_topology.newick","w") as new:
        new.write(top_gene_tree)

    print("writing the transfers / overall tree")
    trans_total_dict = {}
    heat_dict = {}
    for item in topology_dict:
        #ADD TYO TRANS TOTAL DICT:  NAME OF TOPOLOGY DICT : 
        #print(item)
        #will keyerror if str not associated with a transfer event
        try:
         #   print(str_val_dict_trans[item])
            total_trans = str_val_dict_trans[item]
        except KeyError:
            total_trans = 0
        #print(topology_dict[item])
        trans_total_dict[item] = str(total_trans)+"/"+topology_dict[item]
        division = float(total_trans)/float(topology_dict[item])
        heat_dict[item] =str(division)
    #print("heatdict")
    #print(heat_dict)
    trans_total_gene_tree = AddValToTree(gene_tree_string, trans_total_dict)
    with open(projectname+"/AnnotatedTrees/"+prefix+"_transfers_topology.newick","w") as new:
        new.write(trans_total_gene_tree)

    print("writing a heatmapable file")
    heat_total_gene_tree = AddValToTree(gene_tree_string, heat_dict)
    with open(projectname+"/AnnotatedTrees/"+prefix+"_transfers_heat.newick","w") as new:
        new.write(heat_total_gene_tree)
    nexus_file = projectname+"/AnnotatedTrees/"+prefix+"_NEXUS.nexus"
    list_of_dicts_for_write_nexus = [str_val_dict_trans_recip, str_val_dict_trans, topology_dict]
    Write_Annotated_Nexus(gene_tree_string, list_of_dicts_for_write_nexus, overall_tiplist, nexus_file)
    print("wrote a nexus")

    #copy these with dups or loss  --- waiting to implement until I know if trans or trans/total is more effecting. could also show all depending.

    #
    #TODO check if XX0XX is still in the tips replacing the underscores
    #
    
def Write_Annotated_Nexus(gene_tree_string, list_of_dicts_for_write_nexus, overall_tiplist, nexus_file):
    #I think I want to add after each tree [&label="THETHINGTOADD",!name="ANNOTATION"]
    #so this is basically the same as the already existing function.
    #in fact, I can just pass a dictionary of string : [&label="THETHINGTOADD",!name="ANNOTATION"] probably.
    trans_recip = list_of_dicts_for_write_nexus[0]
    trans_node = list_of_dicts_for_write_nexus[1]
    topo = list_of_dicts_for_write_nexus[2]
    complex_dict = {}
    for item in topo:
        complex_dict[item] = "[&topoBS=\""+topo[item]+"\""
        if item in trans_recip:
            complex_dict[item] = complex_dict[item]+",!name=\""+trans_recip[item]+"\""
        if item in trans_node:
            complex_dict[item] = complex_dict[item]+",&transnode=\""+trans_node[item]+"\""
        complex_dict[item] = complex_dict[item]+"]"
    
    #print(complex_dict)
    lentax = len(overall_tiplist)
    print("are there single tips above?")

    THETREESTRING = AddValToTree(gene_tree_string, complex_dict)
    with open(nexus_file, "w") as new:
        new.write("#NEXUS\nbegin taxa;\n\tdimensions ntax="+str(lentax)+";\n\ttaxlabels\n")
        for item in overall_tiplist:
            new.write("\t"+item+"\n")
        new.write(";\nend;\n\nbegin trees;\n\ttree tree_1 = [&R] ")
        new.write(THETREESTRING)
        new.write("\nend;")



#ok so this is fucked rn.
def topology_support_dict(list_inputs, overall_tiplist):
    #input is a list of.. dictionaries... that map STRING: #transfers


    top_dict = {}
    for item in list_inputs:
        for key in item:
            if key not in top_dict:
                top_dict[key] = item[key]
                # (A,b) = 5
            else:
                # (A,b) = 5
                # (A,B) = 5+2
                # (A,B) = 7
                a = (int(top_dict[key])+int(item[key]))
                #print (a)
                #raise SystemExit
                top_dict[key] = str(a)
    #print(top_dict)
    for item in overall_tiplist:
        top_dict[item] = "100"
    return top_dict


def AddValToTree(gene_tree, str_val_dict):
    #gene_tree is a string of a newick
    #str val dict is substring of newick + value to insert
    #
    #ERROR with the orderingh of the replacement strings --- needs to go insiude out.
    #dict is inherentlly unordered. lets make it ordered then i guess.
    list_of_items = []
    for item in str_val_dict:
        list_of_items.append(item)

    comma_list = Order_List_By_Commas(list_of_items)

    for item in comma_list:
        gene_tree = gene_tree.replace(item, item+str_val_dict[item])
    return gene_tree


def Order_List_By_Commas(list_input):
    comma_dict = {}
    for item in list_input:
        commanum = 0
        for character in item:
            if character == ",":
                commanum +=1
        comma_dict[item] = commanum
 
    sorted_x = sorted(comma_dict.items(), key=operator.itemgetter(1), reverse=True)
    #will be a list of sorted tuples, those with MOST commas in front.

    #now turn it into a list
    comma_list = []
    for item in sorted_x:
        comma_list.append(item[0])
    return comma_list


def MakeStrValDict(str_tips_dict, tips_val_dict):
    #str_tips = {"(A,b),C)":[A,B,C]}
    #tips_val_dict = {"A B C":8}
    #some wont match because the given topology either had no transfers, or never gave that topology.
    str_val_dict = {}
    for item in str_tips_dict:
        match = False
        tips = str_tips_dict[item]
        d = sorted(tips)
        string_strtips  = " ".join(d)

        #a is list of tips
        for tips_vald in tips_val_dict:
            if tips_vald  == string_strtips:
                match = True
                #print ("ofound a match")
                #if the string of the sorted lists of tips are identicle;
                #match the string (C,(B,A)) with the value (8)
                #newdict[(c(ba))] = 8 
                str_val_dict[item] = tips_val_dict[tips_vald]
        #if match is False:
            #print("didnt find a match for: "+string_strtips+" in")
            #print(tips_val_dict)
    return str_val_dict

#
                
    #interest should be : "Transfers_By_Recipient_Clade" "Speciation" "Losses" "Duplications" "TRANSFERTIPS" "TIPS THAT RECIEVED TRANSFER AT LEAST 80% of the time:"
def Parse_biparts_file(biparts_file, interest, list_new_options):
    #biparts file is XX0XX_gi
    start = False
    end = False
    tips_val_dict = {}
    print(interest)
    with open(biparts_file) as bi:
        for line in bi:
            if end is True:
                break
            elif line.strip() == "Transfers_By_Representative_Node":
                if start is True:
                    end = True
                    break
                elif line.strip() == interest.strip():
                    start = True
                    print("found: "+interest)
                else:
                    continue
            elif line == "\n":
                if start is True:
                    end = True
                    break

            elif line.strip() == interest.strip():
                start = True
                print("found: "+interest)
            elif start is True:

                line = line.strip()
                try:
                    tips, value, order = line.split("\t")
                except:
                    print(line)
                    raise SystemExit
                #process the tips.
                list_of_old_tips = tips.split(" ")
                list_of_new_tips = []
                for oldtip in list_of_old_tips:
                    newtip = convert_XO_to_fulltipname(oldtip, list_new_options)
                    list_of_new_tips.append(newtip)
                list_of_new_tips.sort()
                #if len(list_of_new_tips) == 1:
                    #print("adding a singletip in biparts parse")
                new_tips_sorted_string = " ".join(list_of_new_tips)
                tips_val_dict[new_tips_sorted_string] = value
             
                #tips, here, are space seperated string like "A B C" = 17
#                MethanosarcinaXX0XXacetivorans MethanosarcinaXX0XXsiciliaeXX0XXHI350 MethanosarcinaXX0XXsiciliaeXX0XXT4XX0XXM MethanosarcinaXX0XXspXX0XX2XX0XXHXX0XXAXX0XX1BXX0XX4 74  NYYYYYYYYYNYYYYYNYYYYYYNNYNNYNYYYNYYYYYYNNNYYYYYYYYYNYNYYNYYYNYNNYYYNYYNYYYYNYYYYYYYYYYYYNYNNYNNYYYY
    if end is False:
        print("something went wrong looking for interst:"+interest+" in file: "+biparts_file)
        
    return tips_val_dict            

#gradschoolthoughts
#friedman+smith worseplace

#bhullar+wagner+near nofish
#coates+slater+westneat+devopeople 


#this is currently set to work on trees that are the output of rangerDTL.
#we... dont, I think, want to do this one them.
#... was this copied over from Run Ranger? YES. ok. fuck. change it. there is no such thing as a node name here. the name is the ENTIRE BLOODY STRING.

def String_to_Tips(string):
    depth = 0
    depth_to_index = {}
    index = 0
    str_to_tips = {}
    list_of_strings = []
    overall_tiplist = []
    for item in string:
        if item == "(":
            depth += 1
            depth_to_index[depth] = index
        if item == ")":
            list_of_strings.append(string[depth_to_index[depth]:index+1])
            depth -= 1
        index += 1
    #print(list_of_strings)
    for item in list_of_strings:
        tiplist = []
        item2 = item.replace("(", "")
        item2 = item2.replace(")", "")
        items2 = item2.split(",")
        #print(items)
        for tip in items2:
            tipsplits = tip.split(":")
            tip = tipsplits[0]
            tiplist.append(tip)
            if tip not in overall_tiplist:
                overall_tiplist.append(tip)
        str_to_tips[item] = tiplist
    #print(str_to_tips)
    #add single tips str mapping to singletip in list.
    for item in overall_tiplist:
        str_to_tips[item] = [item]
    return str_to_tips, overall_tiplist


def convert_XO_to_fulltipname(oldtip, list_new_options):
    
    if "_" in oldtip:
        spec, gi = oldtip.split("_")
        spec = spec.replace("XX0XX", "_")
        for item in list_new_options:
            if spec in item:
                if gi in item:
                    return item
    #old version
    else:
        spec = oldtip.replace("XX0XX", "_")
        for item in list_new_options:
            if spec in item:
                return item
    print("Error, no tip matching:" + oldtip + " was found in " + list_new_options)
    raise SystemExit



def Populate_Bestree(directory, projectname):
    species_files = []
    gene_files = []
    list_of_prefixes = []
    for name in glob.glob(projectname+"/Species_Trees/trees/RAxML_bestTree.*"):
        species_files.append(name)
    for name in glob.glob(projectname+"/Gene_Trees/trees/RAxML_bestTree.*"):
        gene_files.append(name)
    gene_files.sort()
    species_files.sort()
    a = len(projectname)
    print(a)
    for item in species_files:
        print(item)
        iteml = item.split("/")
        item = iteml[-1]
        if "\\" in item:
            x = item.split("\\") 
            item = x[-1]
        #print(item)
        if "_gene" in item:
            item = item[:-5]
        if "_CC" in item:
            item = item[:-3]
        if "RAxML_bestTree" in item:
            item = item[15+a:]
        #print(item)
        list_of_prefixes.append(item)
        print(item)
    #@raise SystemExit
    print(list_of_prefixes)
    new_gene_files = []
    new_species_files = []
    rem_pref = []
    if len(gene_files) == len(species_files):
        print("equal lists")
    else:
        print("uh oh unequal. why?")
        for item in list_of_prefixes:
            #print("matching : "+item)
            got = False
            get = False
            for thing in gene_files:
                #print(thing)
                if item in thing:
                    got = True
                    new_gene_files.append(thing)
                    break
            if got is False:
                print("couldnt match gene file for: "+item)
                rem_pref.append(item)
                continue
            for thing in species_files:
                if item in thing:
                    get = True
                    new_species_files.append(thing)
                    break
            if get is False:
                print("No species file for: "+item)
                new_gene_files.remove(thing)
                rem_pref.append(item)
        for pre in rem_pref:
            list_of_prefixes.remove(pre)
        #print(new_besttree_files)
        gene_files = new_gene_files
        species_files = new_species_files
    print(len(gene_files))
    print(len(species_files))
    print(len(list_of_prefixes))

    if len(gene_files) != len(species_files):
        print(gene_files)
        print(species_files)
        print("still an error in Add_numbers_To_Nodes // Get files // reroot gene besttree to match species besttree")
        raise SystemExit
    #now we've got lists. all of same length and sorted
    #write a ranger in file.
    pop_best_list = []
    for i in range(len(gene_files)):
        prefix = list_of_prefixes[i]
        newname = projectname+prefix+"OptRootInput.txt"
        midname = projectname+prefix+"OptRootOutput.txt"
        outsimpname = projectname+prefix+"simple_rerooted_besttree.newick"
        outcompname = "Best_Rooted_Renamed_Gene_"+projectname+prefix+".newick"

        finalname = One_File_At_A_Time(gene_files[i], species_files[i], newname, midname, outsimpname, outcompname, projectname)
        if finalname != "error":
            pop_best_list.append(finalname)
    return pop_best_list

def One_File_At_A_Time(gene_file, species_file, newname, midname, outsimpname, outcompname, projectname):

    with open(gene_file) as gene:
        for line in gene:
            if line != "":
                gene_tree_string = line
                break
    with open(species_file) as sp:
        for line in sp:
            if line != "":
                species_tree_string = line
                break
    species_tree_string = species_tree_string.strip()
    gene_tree_string = gene_tree_string.strip()
    #print(gene_tree_string)
    print("beginning to make dictionary of string to included tips for file: "+gene_file)
    str_tips_dict, overall_tiplist = String_to_Tips(gene_tree_string)
    new_tipslist, simple_gene_tree_string = Convert_tree_tips_to_simple(overall_tiplist, gene_tree_string)
    print("writing a file from")
    print(os.getcwd())
    print("called")
    print(newname)
    #print(newname)
    with open(newname, "w") as new:
        species_tree_string_XO = species_tree_string.replace("_", "XX0XX")
        species_tree_string_XO = species_tree_string_XO.replace("#", "")
        simple_gene_tree_string = simple_gene_tree_string.replace("#", "")
        new.write("[&R]"+species_tree_string_XO+"\n")
        new.write("[&U]"+simple_gene_tree_string)
    
    #if len(overall_tiplist) > 400:
        #print("too many tips for testrun. skipping for now:"+gene_file)
        #return "error"

    #optroot it:
    if os.path.isfile(projectname+"/"+midname) is False:
        print("NO FILE EXISTS:")
        print(os.getcwd())

        print(projectname+"/"+midname)
        #this depends on system! .win or .mac
        os.system("OptRoot.mac -i "+newname+" -o "+midname)
    else:
        midname = projectname+"/"+midname

    #see if file errored:
    ex = os.path.isfile(midname)
  
    if ex is False:
        print("OPTROOT FAILED. dropping the file: "+midname)
        print("this might just be due to running on windows...? works fine on mac & linux.")
        raise SystemExit
        return "error"
    else:
        empty = os.path.getsize(midname)
        if empty == 0:
            ex = False
    #parse it
    with open (midname) as optout:
        print("Opening file: "+ midname+" to read trees from...")
        tree = ""
        for line in optout:
            if line[0] == "(":
                if tree == "":
                    tree = line
                else:
                    break
    if tree == "":
        print("optrooted tree not found")
        raise SystemExit
    with open(outsimpname, "w") as final:
        final.write(tree)
    complex_tree = Convert_tree_tips_to_complex(overall_tiplist, new_tipslist, tree)
    with open(outcompname, "w") as comp:
        comp.write(complex_tree)
    return outcompname
    

 
def Convert_tree_tips_to_simple(tipslist, tree_string):
    new_tipslist = []
    gilist = gen_gis_list(tipslist)
    speclist = gen_species_lists(tipslist)
    sp_list_XO = []
    for spect in speclist:
        spect2 = spect.replace("_", "XX0XX")
        sp_list_XO.append(spect2)

    if len(gilist) != len(sp_list_XO):
        print("gilist and speclist do not match error")
        raise SystemExit
    for i in range(len(gilist)):
        new_tipslist.append(sp_list_XO[i]+"_"+gilist[i])
    if len(tipslist) != len(new_tipslist):
        print("old and new tipslist do not match error")
        raise SystemExit
    for i in range(len(tipslist)):
        tree_string = tree_string.replace(tipslist[i], new_tipslist[i])
    return new_tipslist, tree_string

def Convert_tree_tips_to_complex(tipslist, new_tipslist, tree_string):
    for i in range(len(tipslist)):
        tree_string = tree_string.replace(new_tipslist[i], tipslist[i])
    return tree_string

def gen_gis_list(idlist):
    gilist = []
    for item in idlist:
                #print(item)
        taxon = re.sub("(.*)(gi#\|?)([0-9]*)(.*)", "\\3", item)
                #print(taxon)
        if "|" in taxon:
            print("TAXON error in gen_gis_lists():" + taxon)
            bleh, taxon = taxon.split("#")
                #print(taxon)
        gilist.append(taxon)
    return gilist

def gen_species_lists(idlist):
    speclist = []
    for item in idlist:
        taxon = re.sub("([^_]*)([A-Z][a-z]*_?[A-Z]?[a-z]*[^\|]*)(.*)", "\\2", item)
        if "|" in taxon:
            tlist = item.split("|")
            taxon = tlist[-2]
            if "|" in taxon:
                print ("TAXON error in gen_species_lists():" + taxon)
        speclist.append(taxon)
    return speclist

#String_to_Tips("(((Python_bivittatus:0.01745612626516041,(Corvus_brachyrhynchos:0.03015243683938841,Alligator_sinensis:0.0189645691768291):0.003921658577364456):0.005184393235545953,(Rattus_norvegicus:0.011688812842066007,(Felis_catus:0.0016682570148515916,Canis_lupus_familiaris:0.0024726969392472042):0.008197203346173576):0.025744562272010385):0.039211054917268215,Poecilia_reticulata:0.039211054917268215);")


#Get_All_Files_To_Run_On("/Users/acaro/Documents/Plane_Boots/", "SOD8")

#Get_All_Files_To_Run_On("/Users/abigailc/Documents/MakeSpeciesTrees/", "SMCEurCyOrder")


print("DONE! PARSER IS BEGINNING!!!!!!!!!")

if __name__ == "__main__":

    print("Running in terminal")
    #imports (not sure we use them all)
    import sys
    import argparse
    import os
    import re
    parser = argparse.ArgumentParser(description="All")

    #optional directory
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in eg MakeSpeciesTree")
    parser.add_argument("-p", "--projectname", action = "store", help="type projectname eg SOD8")
    args = parser.parse_args()


    print("single tree runthrough currently disabled.")
    if args.projectname is False:
        print("hopefully you weren't trying to run in terminal")
        raise SystemExit

    Get_All_Files_To_Run_On(args.directory, args.projectname)


    print("Run is concluded! Hope your files look alright!")
    







#!/usr/bin/python

#abigailc@Actaeon Mar 8 2017
#edit July 7 abigailc@Actaeon 2017

#issue may2: for SOD we are missing verification of anything that was "fix your mistakes-d" like cyanobacteria! FIXED

#this will take input of (list of big clades) (list_of_small_clades)(original_fasta)
#for all big clades, variables include:
#    transfer_recipient_tips
#    other stuff


#end result is a subsampled .fasta file generated from the initial input .fasta that includes:
#-for any majopr clade - 2 divergent tips if no other large or small clades within
# if divergent clade: tips number depends on how deep the transfer is. shallowest: keeps 2. one deeper: keeps 4. more than that: keeps shallow 4 and 2 from sister to transfer.
#-for any minor clade - 2 representatives if possible. CHECK THIS.

import re

def Master_Preform_Subsampling(list_of_big_clades, list_of_small_clades, original_fasta_object, projectname):
    #original fasta tips are all fine.


    #fix your mistakes again
    if projectname[:3] == "SOD":
        sod = True
        print("SOD IS TRUE")
    else:
        sod = False

    #Make Lists of Tips to subsample, to preserve all transfers between groups of the level specified (eg - order.... looking at each run's created gene tree.)
    print("Defining BigClades")
    Define_Bigclade(list_of_big_clades, list_of_small_clades)
    print("Defining SmallClades")
    Define_Smallclade(list_of_big_clades, list_of_small_clades)
    SUBSAMPLED_SEQUENCE_IDS = []
    print("Looking through each large idenfified clade")
    print("taking between 2 and n sequences from each, depending on the number of deep potential transfers identified")
    for item in list_of_big_clades:
        illegals = []
        listtree = item.gene_tree
        needed = item.string_name
        print(needed+" selections:")
        listset = [item.take0, item.take1, item.take2]
        If_Empty_Add_First_Split_To_1(listset, listtree)
        #there will be duplicate lists if it found multiple subclades within a subtree.
        #also maybe we want to not sample deepest split if there are more specified samples requested. 
        #so remove lists that are a superset of any other list
        #and also identical copies
        listset[0] = Remove_Excess_Lists(listset[0])
        listset[1] = Remove_Excess_Lists(listset[1])
        listset[2] = Remove_Excess_Lists(listset[2])

        #this is not tracked anymore, but maybe we want to re-implement it? 
        #Picking tips within a subtree that are likely recipients of transfer
        #removed because ALL TIPS were called "recipients of transfer" and we were seriously messing up due to this.
        
    #    trans_recips = item.transfer_recipient_tips


        #AS A TEST IMMA DISABLE THIS
   #     #for tip in trans_recips:
  #      #    newtip = convert_XO_to_fulltipname(tip, original_fasta_object)
 #       #    illegals.append(newtip)

        #don't take the same tip twice. this shouldn't ever happen with string-matching, but, uh, better safe than sorry
        for tip in SUBSAMPLED_SEQUENCE_IDS:
            illegals.append(tip)

        #make illegal anything set in the "take 0" category. (why?)
        #if listset[0] != []:
        #    for a_list in listset[0]:
        #        for tip in a_list:
        #            illegals.append(tip)
        #take one from each list in listset1
        if listset[1] != []:
            for a_list in listset[1]:
                verified = False
                for thistip in a_list:
                    verified = Verify(thistip, needed, illegals, sod)
                    if verified is True:
                        break
                print("appending "+thistip)
                SUBSAMPLED_SEQUENCE_IDS.append(thistip)
                illegals.append(thistip)
        #take one from each list in listset2
        if listset[2] != []:
            for a_list in listset[2]:
                onedone = False
                for thistip in a_list:
                    verified = Verify(thistip, needed, illegals, sod)
                    if verified is True:
                        if onedone is False:
                            SUBSAMPLED_SEQUENCE_IDS.append(thistip)
                            print("appending "+thistip)
                            illegals.append(thistip)
                            onedone = True
                        else:
                            SUBSAMPLED_SEQUENCE_IDS.append(thistip)
                            print("appending "+thistip)
                            illegals.append(thistip)
                            break
                if onedone is False:
                    print("no tips could be verified? list is:")
                    print(a_list)
                    #raise SystemExit

    #this is a very basic TAKE 2 for the small clades --- because we don't ahve any trees for them!
    #maybe output two files - one with and one without small clades??
    print("Now looking through the list of small clades... identified groups too sparse to make a subtree of")
    print("Keeping two per identified group")
    for item in list_of_small_clades:
        onedone = False
        needed = item.string_name
        print(needed+" small selections")
        for tip in item.tips:
            tip = tip.replace("\'", "")
            tip = tip.replace("\"", "")
            verified = Verify(tip, needed, SUBSAMPLED_SEQUENCE_IDS, sod)
            if verified is True:  
                if onedone is False:
                    SUBSAMPLED_SEQUENCE_IDS.append(tip)
                    print("appending "+tip)
            
                    onedone = True
                else:
                    SUBSAMPLED_SEQUENCE_IDS.append(tip)
                    print("appending "+tip)
       
                    break
           
    print("subsampled a total of "+str(len(SUBSAMPLED_SEQUENCE_IDS))+" sequences")
    #print("an example subsampled seq before fixing")
    #print (SUBSAMPLED_SEQUENCE_IDS[0])
    new_ss = []
    for item in SUBSAMPLED_SEQUENCE_IDS:
        if "\n" in item:
            il = item.split("\n")
            item = il[0]
        item=item.strip("\"")
        item=item.strip("\'")
        print(item)
        new_ss.append(item)
    #print("an example subsampled seq")
    #print(new_ss[0])
    print(len(original_fasta_object.ids))
    print(len(new_ss))
    print("^ number seqs in original fasta")
    #do the extraction and write the new .fasta file.
    original_fasta_object.extract(new_ss)
    original_fasta_object.gen_new_fasta(projectname+"_SS.fasta")
    print("your final subsampled tree (trying to preserve as many deep transfers as possible) is located at: "+projectname+"Overall_Subsampled.fasta")
    #$raise SystemExit
    return projectname+"_SS.fasta"
    


          
######################SUBSAMPLING FUNCTIONS########

def convert_XO_to_fulltipname(oldtip, fasta_object):
    list_new_options = fasta_object.original_ids
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
    print("Error, no tip matching:" + oldtip + " was found in " + fasta_object.str_name)
    raise SystemExit

def Is_X_within_Y(x_tips, y_tips):
    within = True
    for tip in x_tips:
        if tip in y_tips:
            pass
        else:
            within = False
            break
    return within

def Is_X_Exactly_Y(x_tips, y_tips):
    if len(x_tips) != len(y_tips):
           return False
    for tip in x_tips:
        if tip in y_tips:
            pass
        else:
            return False
    for tip in y_tips:
        if tip in x_tips:
            pass
        else:
            return False
    return True

def Define_Bigclade(list_of_big_clades, list_of_small_clades):
    for bigc in list_of_big_clades:
        print("________BEGINNING ON subtree object: "+bigc.prefix)
        bigc.accounted_for = []
        bigc.take0 = []
        bigc.take1 = []
        bigc.take2 = []
        #see if another bigclade resides within
        for potential_sub in list_of_big_clades:
            if potential_sub == bigc:
                continue
            within = Is_X_within_Y(potential_sub.fasta_object.ids, bigc.fasta_object.ids)
            if within is True:
                #see if potential_sub is half of deepest split exactly.
                #see if potential_sub is half of second deepest split exactly.
                #find where deepest split is and record the sistergroup to it.
                #deepest split is generating these tips purely from a gene_tree.
                #maybe the gene_tree is not what i expect it to be.
                print("there is another big clade within this one! "+potential_sub.prefix+" is being tested for optimal subsampling")
                clade1, clade2, clade1_tre, clade2_tre= DeepestSplit(bigc.gene_tree)
                #print("define bigclades: clade1 after deepestsplit[0]")
                #print(clade1[0])
                if clade1 in bigc.accounted_for:
                    print("1 already accounted for - must also be sub to previous bigclade.")
                    pass
                else:
                    within_1 = Is_X_within_Y(potential_sub.fasta_object.ids,clade1)
                    if within_1 is True:
                    #check and see if it is exactly half
                        a1 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade1)
                        if a1 is True:
                            bigc.take2.append(clade2)
                            bigc.accounted_for.append(clade1)
                            bigc.take0.append(clade1)
                        if a1 is False:
                            #go one level deeper
                            ######IS CLADE 1 THE TREE OR THE LIST OF TIPS. I NEED BOTH.
                            clade3, clade4, clade3_tre, clade4_tre = DeepestSplit(clade1_tre)
                            if clade3 in bigc.accounted_for:
                                print("already accounted for")
                                pass
                            else:
                                within_3 = Is_X_within_Y(potential_sub.fasta_object.ids,clade3)
                                if within_3 is True:
                                    #check and see if it is exactly half
                                    a3 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade3)
                                    if a3 is True:
                                        if bigc.take2 == []:
                                            bigc.take2.append(clade2)
                                        bigc.take2.append(clade4)
                                        bigc.take0.append(clade3)
                                    else:
                                        #iterate down until X is Exactly Y.
                                        #answer needs to be either "error" or ([list_of_strings_in_sister_clade],[list_of_string_in_current_clade])
                                        answer = Iterate_Until_X_Is_Y_Or_Fail(potential_sub.fasta_object.ids, clade3_tre)
                                        if answer == "error":
                                            print("it's never exactly half. keeping two from deepest split, two from each of second deepest/")
                                            if bigc.take2 == []:
                                                bigc.take2.append(clade2)
                                            bigc.take2.append(clade3)
                                            bigc.take2.append(clade4)
                                            #done
                                        else:
                                            if bigc.take2 == []:
                                                bigc.take2.append(clade2)
                                            bigc.take2.append(clade4)
                                            bigc.take2.append(answer[0])
                                            bigc.take0.append(answer[1])
                                            
                                else:
                                    if clade4 in bigc.accounted_for:
                                        print("already accounted for")
                                        pass
                                    within_4 = Is_X_within_Y(potential_sub.fasta_object.ids,clade4)
                                    if within_4 is True:
                                        #check and see if it is exactly half
                                        a4 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade4)
                                        if a4 is True:
                                            if bigc.take2 == []:
                                                bigc.take2.append(clade2)
                                            bigc.take2.append(clade3)
                                            bigc.take0.append(clade4)
                                        else:
                                            #iterate down until X is Exactly Y.
                                            #answer needs to be either "error" or ([list_of_strings_in_sister_clade],[list_of_string_in_current_clade])
                                            answer = Iterate_Until_X_Is_Y_Or_Fail(potential_sub.fasta_object.ids, clade3_tre)
                                            if answer == "error":
                                                print("it's never exactly half. keeping two from deepest split, + 2 from each second deepest")
                                                if bigc.take2 == []:
                                                    bigc.take2.append(clade2)
                                                bigc.take2.append(clade3)
                                                bigc.take2.append(clade4)
                                                
                                                #done
                                            else:
                                                if bigc.take2 == []:
                                                    bigc.take2.append(clade2)
                                                bigc.take2.append(clade4)
                                                bigc.take2.append(answer[0])
                                                bigc.take0.append(answer[1])
                    else:
                        if clade2 in bigc.accounted_for:
                            print("2 already accounted for - must be sub to previous bigclade.")
                            pass
                        else:
                            within_2 = Is_X_within_Y(potential_sub.fasta_object.ids,clade2)
                            if within_2 is True:
                            #check and see if it is exactly half
                                a2 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade2)
                                if a2 is True:
                                    if bigc.take2 == []:
                                        bigc.take2.append(clade1)
                                    bigc.accounted_for.append(clade2)
                                    bigc.take0.append(clade2)
                                if a2 is False:
                                #go one level deeper
                                ######IS CLADE 1 THE TREE OR THE LIST OF TIPS. I NEED BOTH.
                                    clade5, clade6, clade5_tre, clade6_tre = DeepestSplit(clade2_tre)
                                    if clade5 in bigc.accounted_for:
                                        print("already accounted for")
                                        pass
                                    else:
                                        within_5 = Is_X_within_Y(potential_sub.fasta_object.ids,clade5)
                                        if within_5 is True:
                                        #check and see if it is exactly half
                                            a5 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade5)
                                            if a5 is True:
                                                if bigc.take2 == []:
                                                    bigc.take2.append(clade1)
                                                bigc.take2.append(clade6)
                                                bigc.take0.append(clade5)
                                            else:
                                            #iterate down until X is Exactly Y.
                                            #answer needs to be either "error" or ([list_of_strings_in_sister_clade],[list_of_string_in_current_clade])
                                                answer = Iterate_Until_X_Is_Y_Or_Fail(potential_sub.fasta_object.ids, clade5_tre)
                                                if answer == "error":
                                                    print("it's never exactly half. keeping two from deepest split, and two from each of second deepests")
                                                    if bigc.take2 == []:
                                                        bigc.take2.append(clade1)
                                                    bigc.take2.append(clade5)
                                                    bigc.take2.append(clade6)
                                                #done
                                                else:
                                                    if bigc.take2 == []:
                                                        bigc.take2.append(clade1)
                                                    bigc.take2.append(clade6)
                                                    bigc.take2.append(answer[0])
                                                    bigc.take0.append(answer[1])

                                        else:
                                            if clade6 in bigc.accounted_for:
                                                print("already accounted for")
                                                pass
                                            within_6 = Is_X_within_Y(potential_sub.fasta_object.ids,clade6)
                                            if within_6 is True:
                                            #check and see if it is exactly half
                                                a6 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade6)
                                                if a6 is True:
                                                    if bigc.take2 == []:
                                                        bigc.take2.append(clade1)
                                                    bigc.take2.append(clade5)
                                                    bigc.take0.append(clade6)
                                                else:
                                                #iterate down until X is Exactly Y.
                                                #answer needs to be either "error" or ([list_of_strings_in_sister_clade],[list_of_string_in_current_clade])
                                                    answer = Iterate_Until_X_Is_Y_Or_Fail(potential_sub.fasta_object.ids, clade6_tre)
                                                    if answer == "error":
                        
                                                        print("it's never exactly half. keeping two from deepest split, and two from each of second deepest splits")
                                                        if bigc.take2 == []:
                                                            bigc.take2.append(clade1)
                                                        bigc.take2.append(clade5)
                                                        bigc.take2.append(clade6)
                                                    #done
                                                    else:
                                                        if bigc.take2 == []:
                                                            bigc.take2.append(clade1)
                                                        bigc.take2.append(clade5)
                                                        bigc.take2.append(answer[0])
                                                        bigc.take0.append(answer[1])
                            else:
                                print("Uh, the inner bigclade was split up even at the deepest split! Taking two from each side of deepest!")
                                bigc.take2.append(clade1)
                                bigc.take2.append(clade2)
        if bigc.take2 == []:
            clade1, clade2, clade1_tre, clade2_tre= DeepestSplit(bigc.gene_tree)
            bigc.take2.append(clade1)
            bigc.take2.append(clade2)

    #results:
    #if nothing within bigclade: two from each side.
    #otherwise: somewhere between 2 and 8 from each side (or more if multiple shallow clades I guess)

    #other bigclades will do their own sampling.

def Iterate_Until_X_Is_Y_Or_Fail(xtips, ytree):
    #x is within y but not exactly it
    #so do deepest split, check if in, check if half
    #return list_to_take_2, list_to_take_0
    done = False
    current_clade = ytree
    while done is False:
        cladeA, cladeB, cladeA_tre, cladeB_tre = DeepestSplit(current_clade)
        withinA = Is_X_within_Y(xtips, cladeA)
        withinB = Is_X_within_Y(xtips, cladeB)
        if withinA is True:
            exactlyA = Is_X_Exactly_Y(xtips, cladeA)
            if exactlyA is True:
                return cladeB,cladeA
                #do the save thing
                #and return
            else:
                current_clade = cladeA_tre
                continue
        if withinB is True:
            exactlyB = Is_X_Exactly_Y(xtips, cladeB)
            if exactlyB is True:
                return cladeA, cladeB
                #do the save thing
                #and return
            else:
                current_clade = cladeB_tre
                continue
        else:
            return("error")

#takes in a list of lists.
#removes lists that contain all elements of a smaller list.
def Remove_Excess_Lists(listset):
    #ERROR if two lists are the same, don't nuke both of them!!!!
    #thus let us sort our lists so they will be exact matches.


    to_remove = []
    to_remove_all_but_one = []
    for current_list in listset:
        #test if the current list is within any other list.
        for potential_large_list in listset:
            if len(potential_large_list) < len(current_list):
                continue
            #if they are the same remove only one!
            current_in_potential = all(x in potential_large_list for x in current_list)
            if current_in_potential is True:
                if sorted(potential_large_list) == sorted(current_list):
                    if sorted(potential_large_list) in to_remove_all_but_one:
                        pass
                    else:
                        to_remove_all_but_one.append(sorted(potential_large_list))
                elif potential_large_list in to_remove:
                    pass
                else:
                    to_remove.append(potential_large_list)
    #remove items that are a superset of a more exact list
    for item in to_remove:
        listset.remove(item)
    #remove items that are exactly identical to a previous list
    #this is going to cause trouble with skipping, since removal from a list that is iterating. make a copy list to run removal on.
    listset2 = []
    #^copy. return this.
    for item in listset:
        listset2.append(item)

    for item in to_remove_all_but_one:
        foundone = False
        for ls in listset:
            if item == sorted(ls):   
                if foundone is True:
                    listset2.remove(ls)
                   
                else:
                    #this is matching with itself
                    foundone = True
              

    return listset2

#if no subtree-identified transfers within this clade, just subsample one from each side of the deepest split
def If_Empty_Add_First_Split_To_1(list_of_lists, subtree_tre):
    #list of lists needs to be passed in [list0, list1, list2]
    empty = True
    #see if all of the lists are totally empty.
    for item in list_of_lists:
        if item == []:
            pass
        else:
            empty = False
            break
    #if they aren't empty, thats fine go ahead
    if empty is False:
        return list_of_lists
    #if they are ALL empty, add deepest split to the get-1 list.
    if empty is True:
        cladeA, cladeB, cladeA_tre, cladeB_tre = DeepestSplit(subtree_tre)
        list_of_lists[1].append(cladeA, cladeB)
        return list_of_lists


#gets the stuff AROUND the small clade, but does not directly sample from it.    
#this is the same as bigclade but draws potential subs from list_small isntead of list_large
def Define_Smallclade(list_of_big_clades, list_of_small_clades):
    #all the lists were already generated by define_bigclade
    #all we are looking to add is any necessary polarizing sequences.
    for bigc in list_of_big_clades:
        print("looking for small clades residing within: "+bigc.prefix)
        #see if another bigclade resides within
        for potential_sub in list_of_small_clades:
            #print(potential_sub.fastaname)
            #print(potential_sub.fasta_object)
            #is there is an error here, it is because the specified small clades doesn't have fasta object ids pre
            within = Is_X_within_Y(potential_sub.fasta_object.ids, bigc.fasta_object.ids)
            if within is True:
                #see if potential_sub is half of deepest split exactly.
                #see if potential_sub is half of second deepest split exactly.
                #find where deepest split is and record the sistergroup to it.
                clade1, clade2, clade1_tre, clade2_tre= DeepestSplit(bigc.gene_tree)
                if clade1 in bigc.accounted_for:
                    print("already accounted for - must be sub to previous bigclade.")
                    pass
                else:
                    within_1 = Is_X_within_Y(potential_sub.fasta_object.ids,clade1)
                    if within_1 is True:
                    #check and see if it is exactly half
                        a1 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade1)
                        if a1 is True:
                            bigc.take2.append(clade2)
                            bigc.accounted_for.append(clade1)
                            bigc.take0.append(clade1)
                        if a1 is False:
                            #go one level deeper
                            ######IS CLADE 1 THE TREE OR THE LIST OF TIPS. I NEED BOTH.
                            clade3, clade4, clade3_tre, clade4_tre = DeepestSplit(clade1_tre)
                            if clade3 in bigc.accounted_for:
                                print("already accounted for")
                                pass
                            else:
                                within_3 = Is_X_within_Y(potential_sub.fasta_object.ids,clade3)
                                if within_3 is True:
                                    #check and see if it is exactly half
                                    a3 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade3)
                                    if a3 is True:
                                        if bigc.take2 == []:
                                            bigc.take2.append(clade2)
                                        bigc.take2.append(clade4)
                                        bigc.take0.append(clade3)
                                    else:
                                        #iterate down until X is Exactly Y.
                                        #answer needs to be either "error" or ([list_of_strings_in_sister_clade],[list_of_string_in_current_clade])
                                        answer = Iterate_Until_X_Is_Y_Or_Fail(potential_sub.fasta_object.ids, clade3_tre)
                                        if answer == "error":
                                            print("it's never exactly half. keeping two from deepest split, and two from each of second deepest splits")
                                            if bigc.take2 == []:
                                                bigc.take2.append(clade2)
                                            bigc.take2.append(clade3)
                                            bigc.take2.append(clade4)
                                        else:
                                            if bigc.take2 == []:
                                                bigc.take2.append(clade2)
                                            bigc.take2.append(clade4)
                                            bigc.take2.append(answer[0])
                                            bigc.take0.append(answer[1])                                         
                                else:
                                    if clade4 in bigc.accounted_for:
                                        print("already accounted for")
                                        pass
                                    within_4 = Is_X_within_Y(potential_sub.fasta_object.ids,clade4)
                                    if within_4 is True:
                                        #check and see if it is exactly half
                                        a4 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade4)
                                        if a4 is True:
                                            if bigc.take2 == []:
                                                bigc.take2.append(clade2)
                                            bigc.take2.append(clade3)
                                            bigc.take0.append(clade4)
                                        else:
                                            #iterate down until X is Exactly Y.
                                            #answer needs to be either "error" or ([list_of_strings_in_sister_clade],[list_of_string_in_current_clade])
                                            answer = Iterate_Until_X_Is_Y_Or_Fail(potential_sub.fasta_object.ids, clade3_tre)
                                            if answer == "error":
                                                print("it's never exactly half. keeping two from deepest split, and two from each of second deepest splits")
                                                if bigc.take2 == []:
                                                    bigc.take2.append(clade2)
                                                bigc.take2.append(clade3)
                                                bigc.take2.append(clade4)
                                                #done
                                            else:
                                                if bigc.take2 == []:
                                                    bigc.take2.append(clade2)
                                                bigc.take2.append(clade4)
                                                bigc.take2.append(answer[0])
                                                bigc.take0.append(answer[1])
                    else:
                        if clade1 in bigc.accounted_for:
                            print("already accounted for - must be sub to previous bigclade.")
                            pass
                        else:
                            within_2 = Is_X_within_Y(potential_sub.fasta_object.ids,clade2)
                            if within_2 is True:
                            #check and see if it is exactly half
                                a2 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade2)
                                if a2 is True:
                                    if bigc.take2 == []:
                                        bigc.take2.append(clade1)
                                    bigc.accounted_for.append(clade2)
                                    bigc.take0.append(clade2)
                                if a2 is False:
                                    #go one level deeper
                                    ######IS CLADE 1 THE TREE OR THE LIST OF TIPS. I NEED BOTH.
                                    clade5, clade6, clade5_tre, clade6_tre = DeepestSplit(clade2_tre)
                                    if clade5 in bigc.accounted_for:
                                        print("already accounted for")
                                        pass
                                    else:
                                        within_5 = Is_X_within_Y(potential_sub.fasta_object.ids,clade5)
                                        if within_5 is True:
                                            #check and see if it is exactly half
                                            a5 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade5)
                                            if a5 is True:
                                                if bigc.take2 == []:
                                                    bigc.take2.append(clade1)
                                                bigc.take2.append(clade6)
                                                bigc.take0.append(clade5)
                                            else:
                                                #iterate down until X is Exactly Y.
                                                #answer needs to be either "error" or ([list_of_strings_in_sister_clade],[list_of_string_in_current_clade])
                                                answer = Iterate_Until_X_Is_Y_Or_Fail(potential_sub.fasta_object.ids, clade5_tre)
                                                if answer == "error":
                                                    print("it's never exactly half. keeping two from deepest split, and two from each of second deepest splits")
                                                    if bigc.take2 == []:
                                                        bigc.take2.append(clade1)
                                                    bigc.take2.append(clade5)
                                                    bigc.take2.append(clade6)
                                                else:
                                                    if bigc.take2 == []:
                                                        bigc.take2.append(clade1)
                                                    bigc.take2.append(clade6)
                                                    bigc.take2.append(answer[0])
                                                    bigc.take0.append(answer[1])

                                        else:
                                            if clade6 in bigc.accounted_for:
                                                print("already accounted for")
                                                pass
                                            within_6 = Is_X_within_Y(potential_sub.fasta_object.ids,clade6)
                                            if within_6 is True:
                                                #check and see if it is exactly half
                                                a6 = Is_X_Exactly_Y(potential_sub.fasta_object.ids, clade6)
                                                if a6 is True:
                                                    if bigc.take2 == []:
                                                        bigc.take2.append(clade1)
                                                    bigc.take2.append(clade5)
                                                    bigc.take0.append(clade6)
                                                else:
                                                #iterate down until X is Exactly Y.
                                                #answer needs to be either "error" or ([list_of_strings_in_sister_clade],[list_of_string_in_current_clade])
                                                    answer = Iterate_Until_X_Is_Y_Or_Fail(potential_sub.fasta_object.ids, clade6_tre)
                                                    if answer == "error":
                                                        print("it's never exactly half. keeping two from deepest split, and two from each of second deepest splits")
                                                        if bigc.take2 == []:
                                                            bigc.take2.append(clade1)
                                                        bigc.take2.append(clade5)
                                                        bigc.take2.append(clade6)
                                                    #done
                                                    else:
                                                        if bigc.take2 == []:
                                                            bigc.take2.append(clade1)
                                                        bigc.take2.append(clade5)
                                                        bigc.take2.append(answer[0])
                                                        bigc.take0.append(answer[1])
                            else:
                                print("small clade split up even in first split; taking two from each side of deepest split")
                                bigc.take2.append(clade1)
                                bigc.take2.append(clade2)


def DeepestSplit(tree):
    #provided subtree will be text. like
    #(((((((((((Elephantulus_edwardii:0.02656660969,Orycteropus_afer:0.006886559777):0.002186555754,Trichechus_manatus:0.0207867919)89:0.0042892984,(Erinaceus_europaeus:0.0150776695,Sus_scrofa:0.01871033399)44:0.002137367876)98:0.001261277757,Equus_asinus:0.006995816459)11:0.0003479426504,((Pteropus_vampyrus:0.01260656735,(Pan_paniscus:0.007271418717,Manis_javanica:0.01036292681)13:0.001322734885)9:0.001552612205,(Dasypus_novemcinctus:0.009902206268,Tupaia_chinensis:0.02113855241)30:0.001364556134)12:0.001840714523)4:0.001979025907,Galeopterus_variegatus:0.00671926141)25:0.004692110764,Ochotona_princeps:0.02065394194)71:0.005193088649,Mus_musculus:0.01442576505)60:0.02036347729,(Monodelphis_domestica:0.007867166076,Sarcophilus_harrisii:0.01059171724)100:0.02083374525)100:0.007086087413,Ornithorhynchus_anatinus:0.03185159824)65:0,Gekko_japonicus:0.04818786247)Root;
    #returns: ["tipa1","tipa2", "tipa3"],["tipb1","tipb2","tipb3"]
    # where the first list represents strings in the first clade made by deepest split, and second list represents strings in the other clade.
    
    indent = 0
    first = ""
    switch = False
    last = ""
    if "Root;" in tree:
        tree.replace("Root;", ";")
    for character in tree:
        if character == "(":
            indent +=1
        if character == ",":
            if indent == 1:
                switch = True
            else:
                indent = indent - 1
        if switch is False:
            first = first+character
        else:
            last = last+character
    first = first[1:] #cuts the initial paren
    last = last[1:-2] #cuts seperating comma and final paren +;
    #now we have the two deepest subtrees.
    #lets extract the tip names.
    first_edit = re.sub("(:[^,]*)", "", first)
    first_edit = re.sub("[\(\)]", "", first_edit)
    first_edit_list = first_edit.split(",")
    last_edit = re.sub("(:[^,]*)", "", last)
    last_edit = re.sub("[\(\)]", "", last_edit)
    last_edit_list = last_edit.split(",")
    #returns list of tip_names in side1 of deepest split, then side2.
    # rets [a,b,c], [d,e,f], "(a,(b,c)", "(d,e),f)"
    first = first+";"
    last = last+";"


    if "\n" in first:
        print("seems like deepest split is still wonky.")
        print("first is: "+first)
    return first_edit_list, last_edit_list, first, last

# verified = Verify(thistip, needed, illegals)
#eg tip = "Cyano|Nostocales|Nostoc_punctiforme_PCC_73102|gi#|12345"
# needed = "Nostocales"
# illegals = list of tips that are transfer recipients
def Verify(tip, needs, illegals, sod = False):
    #print(illegals)
    #make needs FIX YOUR MISTAKES.
   # if "anobac" in needs:
    #    verb = True
     #   print(needs)
    #else:
    verb = False

    #fix your mistakes
    if sod is True:
        if len(needs) > 6:
            if verb is True:
                print (needs[-7:])
            if "mycetes" in needs[-7:]:
                needs = re.sub("mycetes", "myc", needs)
        if len(needs) > 7:    
            if verb is True:
                print (needs[-8:])        
            if "bacteria" in needs[-8:]:
                needs = re.sub("bacteria", "bac", needs)
        if len(needs) > 8:
            if verb is True:
                print (needs[-9:])
            if "mycetales" in needs[-9:]:
                needs = re.sub("mycetales", "mycl", needs)
        if len(needs) > 10:
            if verb is True:
                print (needs[-11:])
            if "bacteriales" in needs[-11:]:
                needs = re.sub("bacteriales", "bacl", needs)
        
    if verb is True:
        print("testing: "+tip+" for : "+needs)
        a = raw_input("waiting...")
    if needs in tip:
        if tip in illegals:
            print(tip+" in illegals")
            return False
        else:
            return True
    return False
    

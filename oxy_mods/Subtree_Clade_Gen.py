#!/usr/bin/python

#abigailc@Actaeon Jan4 2017

#module to handle the subtree parsing and clade determination.
import sys
import argparse
import os
import re
import time

#keep
def MakeSubtreesFile(nexus, subtrees, verbose = False):
    print("Making subtreesfile")
    isfile = os.path.isfile("./"+subtrees)
    if isfile == True:
        print("There is already a subtrees file made from this subtree! \n it is called: "+subtrees+".")
        print("You can indicate it with the -t flag to save time")
        answer = raw_input("Type \"continue\" if you would like to proceed anyways.  :")
        print(answer)
        if answer == "continue":
              pass
        else:
              print("Quitting")
              raise SystemExit
    os.system("subtreegen "+nexus+" "+subtrees)
    return subtrees

#keep
def GetRankFromSubtreesFile(subtrees, ranknum, verbose = False):
    print("Getting ranks from subtrees file seqIDS")
    ranklist = []
    lines = []
    new = "no"
    newtree = "no"
    with open(subtrees) as old:
        for line in old:
            if line == "SUBTREES_BEGIN\n":
                print("beginning...")
            elif new == "yes":
                new = "no"
                newtree = "yes"
            elif newtree == "yes":
                newtree = "no"
            elif line == "subtree#\n":
                new = "yes"
            else:
                if "\n" in line:
                    line = line[:-1]
                
                if line in lines:
                    pass
                else:
                    lines.append(line)
##    print(lines)
    #also stop double tree at the front of paraphyly found trees see streptophytes
    error_s = 0
    for example in lines:
        import re
        
        gitax = re.sub("(.*)(\|)(gi#?\|?)([0-9]*)(.*)", "\\4~\\1", example)
##  this bit only works if your seqID format is taxon|omic|info|gi#|whateverelse
##if switch to accession nums, replace gi(above) with acn or something.
        
        try:
            gi, tax = gitax.split("~")
        except:
            print ("ERROR WITH:"+gitax)
            print("Taxon has no gi number. rank assignation may be confused... we will see.")
            tax = gitax
        tn = tax[:-1]
        taxlist = tn.split("|")
        therankstr = taxlist[ranknum-1]
        therankstr = therankstr+"|"

        if therankstr in ranklist:
            pass
        else:
            if therankstr == "X|":
                error_s += 1
            elif therankstr == "NA|":
                error_s += 1
            elif str == "Proteobac|":
                print("PROTEOBAC WAS LOCATED WHAT WHY")
                print(example)
            else:
                ranklist.append(therankstr)
##        print(ranklist)
    print("Added "+str(len(ranklist))+" groupnames to parsing queue")
    print("There were "+str(error_s)+" NA or X missing taxonomic information seqs that were ignored.")
    return ranklist

#keep

#this is dict of all subtrees. not named or looked at taxonomy yet.
def MakeSubtreesDict(subtrees, verbose = False):
    print("Making Subtrees Dictionary")
    stdict = {}
    sttips = []
    st_trees = {}
    new_tree = "no"
    new = "no"
    first = "yes"
    #print(subtrees)
    with open(subtrees) as old:
        for line in old:
            if line == "SUBTREES_BEGIN\n":
                if verbose == True:
                    print("beginning subtree parsing")
            elif new_tree == "yes":
                st_trees[stname] = line[:-1]
                new_tree = "no"
            elif new == "yes":
                if first == "yes":
                    stname = "1"
                    first = "no"
                    new = "no"
                    new_tree = "yes"
                else:
                    #add name to dict with tips as ke
##                    print(stname+" has tips: "+str(len(sttips)))
                    stdict[stname] = sttips
                    #reset tips
                    sttips = []
                    #resetname
                    stname = line
                    new = "no"
                    if "\n" in stname:
                        stname = stname[:-1]
                    new_tree = "yes"
            elif line == "subtree#\n":
                new = "yes"
            else:
                if "\n" in line:
                    line = line[:-1]
                sttips.append(line)
        #catch the final subtre
        stdict[stname] = sttips
    print("Subtree dict finished")
##    print (stdict)
    return stdict, st_trees


#new
def verify_position(string, tip, ranknum):
    try:
        ranknum = int(ranknum)
        
    except:
        print("cant")
        print(ranknum)
    gitax = re.sub("(.*)(\|)(gi#?\|?)([0-9]*)(.*)", "\\4~\\1", tip)
    try:
        gi, tax = gitax.split("~")
    except:
        print ("ERROR WITH:"+gitax)
        print("Taxon has no gi number. rank assignation may be confused... we will see.")
        tax = gitax
    tn = tax[:-1]
    taxlist = tn.split("|")
    therankstr = taxlist[ranknum-1]
    therankstr = therankstr+"|"
    if string == therankstr:
        return True
    else:
        #print("inexact match")
        #print(ranknum)
        #print(string)
        #print(therankstr)
        return False
#keep
def CheckMonophyly(string, stdict, fewnum, ranknum, verbose = False):
    #verbose = True
    #fortesting
    #why is ranknum passing in as FALSE
    goodst = {}
    alltips = []
    for subtree in stdict:
        tips = stdict[subtree]
        goodtips = 0
        status = "good"
        for tip in tips:
            if string in tip:
                #we should verify it is in the correct position. eg: proteobacteria not deltaproteobacteria.
                #COMEBACKHERE
                a = verify_position(string, tip, ranknum)
                if a is True:
                    goodtips+=1
                    if tip in alltips:
                        pass
                    else:
                        alltips.append(tip)
                else:
                    pass
            else:
                status = "bad"
                
        if status == "good":
            goodst[subtree] = goodtips
            #break
            ##to save time - returns fewnum here rather than calculating later.
            #BUT this means we have to loop through all tips adding to alltips here
            #wheras before, we just looped until a bad status and then broke, which was
            #much faster. unclear which way is better in long run, but this works for now.
            #negligible save in >2000 seq datasets
    if len(alltips) < fewnum:
        if verbose == True:
           print(str(len(alltips))+string+"is not enough tips to reliably assign clades")
        return("few")
    if goodst == {}:
        if verbose == True:
            print(string+"is not monophyletic")
        return "NotMono"
    #get largest monophyletic subtree
    if verbose == True:
        print(goodst)
    ke, val = max(goodst.iteritems(), key=lambda x:x[1])
    if verbose == True:
        print(ke)
        print(val)
    #is this where we are failing???
            #ke = max(ex.keys(), key = lambda k: ex[k])
    #ensure all tips containing string are in largest mono clade
    mono = "yes"
    for item in alltips:
        if item in stdict[ke]:
            pass
        else:
            if verbose == True:
                print(string+" has monophyletic clades, but they do not contain ALL the tips of type str")
            mono = "no"
            break
    if mono == "yes":
        if verbose == True:
           print(string+": all sequences are in one clade and monophyletic")
        return (string, "Monophyly_Strict", ke, str(goodst[ke]), "0", str(len(alltips)))
    else:
        return "NotMono"

#    return "type, subtree, goodtips, cladetips, outsidetips"

#called from within check artifacts
def CheckMonoSingleSubtree(subtreename, subtreetipslist, alltipslist, verbose = False):
    mono = "yes"
    for item in subtreetipslist:
        if item in alltipslist:
            pass
        else:
            mono = "no"
            break
    if mono == "yes":
        return "Monophyly"
    else:
        return "NotMono"

#keep
def CheckInclusions(string, stdict, ranknum, few = "no", verbose = False):
    valst = {}
    totst = {}
    goodst = {}
    alltips = []
    for subtree in stdict:
##        print("Subtree: "+subtree)
        tips = stdict[subtree]
        
        goodtips = 0
        badtips = 0
        for tip in tips:
            
            if string in tip:
                a = verify_position(string, tip, ranknum)
                if a is True:
                    goodtips+=1
                    if tip in alltips:
                        pass
                    else:
                        alltips.append(tip)
                else:
                    badtips+=1
            else:
                badtips+=1
        #weight subtrees by number strtips * percent of tips that are strtips
        #make dictionaries of subtreename: total tips, goodtips, %goodtips, weightedvalue
        tottips = goodtips+badtips
        totst[subtree] = tottips
        goodst[subtree] = goodtips
   
        perctips = float(goodtips)/float(tottips)
       
        valtips = goodtips*perctips
        
        valst[subtree] = valtips
##        if str(goodtips) == "0":
##            pass
##        else:
##            print("good"+str(goodtips))
##            print("bad"+str(badtips))
##            print("total"+str(tottips))
##            print("perc"+str(perctips))
##            print("Val"+str(valtips))
    #get largest value subtree.
    
    ke, val = max(valst.iteritems(), key=lambda x:x[1])
    if verbose == True:
        print("Largest subtree: "+str(ke)+" "+str(val))
            #ke = max(ex.keys(), key = lambda k: ex[k])
    #ensure all tips containing string are in largest mono(w inclusions) clade
    mono = "yes"
    for item in alltips:
        if item in stdict[ke]:
            pass
        else:
            mono = "no"
            break
    if few == "few":
        return ("few", ke, goodst, totst, alltips, valst, stdict)
    if mono == "yes":
##        print ("Subtree number "+ke+" is monophyletic with inclusions with respect to "+string+" .\nIt has "+totst[ke]+" tips total, of which "+goodst[ke]+" are "+string+".")
##set a limit on how big inclusions can be?? else we get the case where it returns subtree 1: bacteria(?)

        return (string, "Paraphyly_Strict", ke, str(goodst[ke])+"/"+str(totst[ke]), "0", str(len(alltips))) 
    else:
        return ("NotParaAlone", ke, goodst, totst, alltips, valst, stdict)
    #needs a lot of testing to ensure we always want num*percent as biggest clade. what if huge inclusion?

#artcutoff  =  (str tips in clade) / (total str tips) for determining monophyly with artifacts. eg (.7 for jumpy dataset, .9 for solid)?
#duocutoff = the percent tips-in-clade to total-str-tips a subtree must have to be considered likely true clade. eg(.1 to get many, .3 to get big ones)?ArithmeticError

#keep
def CheckArtifacts(string, stdict, monoinc, artcutoff, fewval = "no", verbose = False):
    dupcutoff = 1.0 - float(artcutoff)
    bestsubtree = monoinc[1]
    goodst = monoinc[2]
    goodtipsst = goodst[bestsubtree]
    totst = monoinc[3]
    totaltipsst = totst[bestsubtree]
    allstrtips = monoinc[4]
    valst = monoinc[5]
    stdict = monoinc[6]
    bestartval = float(goodtipsst)/float(len(allstrtips))
    if fewval == "few":
        cmono = CheckMonoSingleSubtree(bestsubtree, stdict[bestsubtree], allstrtips)
        if cmono == "Monophyly":
            return(string, "Singleton", "best:"+bestsubtree+"M", str(totst[bestsubtree]), str(len(allstrtips)-int(goodst[bestsubtree])), str(len(allstrtips)))
        else:
            return(string, "Singleton", "best:"+bestsubtree+"P", str(goodst[bestsubtree])+"/"+str(totst[bestsubtree]), str(len(allstrtips)-int(goodst[bestsubtree])), str(len(allstrtips)))
    elif bestartval > float(artcutoff):
        if goodtipsst == totaltipsst:
            if verbose == True:
                print("Mono with art")
            return(string, "Mono_Few_Exc", bestsubtree, str(goodtipsst), str(len(allstrtips)-goodtipsst), str(len(allstrtips)))
        else:
            if verbose == True:
                print("Para with art")
            return(string, "Para_Few_Exc", bestsubtree, str(goodtipsst)+"/"+str(totaltipsst), str(len(allstrtips)-goodtipsst), str(len(allstrtips)))
    else:
        if verbose == True:
            print("Not mono or para with Artifacts (cutoff="+str(artcutoff)+"), best tree ="+str(bestartval))
        numberclades = 1
        if verbose == True:
            print("dc="+str(dupcutoff))
        threshold = float(len(allstrtips))*float(dupcutoff)
        if goodst[bestsubtree] > threshold:
            pass
        else:
            if verbose == True:
                print("No clade meets min threshold of "+str(dupcutoff)+" of total group tips")
            if verbose == True:
                print("Best clade is: "+str(goodst[bestsubtree])+" but "+str(threshold)+" is needed.")
            cmono = CheckMonoSingleSubtree(bestsubtree, stdict[bestsubtree], allstrtips)
            if cmono == "Monophyly":
                return(string, "NoClade_thresh", "best:"+bestsubtree+"M", str(totst[bestsubtree]), str(len(allstrtips)-int(goodst[bestsubtree])), str(len(allstrtips)))
            else:
                return(string, "NoClade_thresh", "best:"+bestsubtree+"P", str(goodst[bestsubtree])+"/"+str(totst[bestsubtree]), str(len(allstrtips)-int(goodst[bestsubtree])), str(len(allstrtips)))
        cladeslist = [bestsubtree]
        itis = "on"
        notallowed = stdict[bestsubtree]
        if verbose == True:
            print("Best st is: "+str(bestsubtree)+" with tips: "+str(len(notallowed)))
        liststdict = []
        for i in stdict:
            exi = i
            liststdict.append(i)
        if verbose == True:
            print("There are "+str(len(liststdict))+" distinct subtrees to start")
        while itis == "on":
            #remove subtrees from consideration if they contain tips from already considered subtrees.
            if verbose == True:
                print("Removing intersecting subtrees based on "+str(len(notallowed))+" used tips")
            sttotal = 0
            stremoved = 0
            lis2 = []
            for item in liststdict:
                lis2.append(item)
            for subtreel in liststdict:
                sttotal +=1
                ok = "yes"
                for tip in stdict[subtreel]: 
                    for ntip in notallowed:
                        if tip == ntip: 
                            ok = "no"
##                            print(valst[subtreel])
##                            print("removing a subtree named"+str(subtreel)+"of type"+str(type(subtreel)))
                            if subtreel in lis2:
                                stremoved +=1
                                lis2.remove(subtreel)
                                del valst[subtreel]
                            break
                    if ok == "no":
                        break
            liststdict = lis2
            if verbose == True:
                print("Start: "+str(sttotal)+" Removed: "+str(stremoved)+" Considering "+str(len(liststdict))+" exclusive subtrees")
            #if all trees contain tip, and all trees overlap, valst will be entirely removed. skipp further verification, will return one cladee tree later. later
            if len(valst) == 0:
                itis = "off"
            else:
                #get topmost non-intersecting tree
                ke, val = max(valst.iteritems(), key=lambda x:x[1])
                #see if it beats threshold to be considered a clade in its own right
                #if yes, continue search. if no close out.
                if goodst[ke] > threshold:
                    numberclades+=1
                    if verbose == True:
                        print("Found a clade. Num "+str(ke)+"Val: "+str(valst[ke])+"Good: "+str(goodst[ke])+"Thr: "+str(threshold)+"total: "+str(numberclades))
                    cladeslist.append(ke)
                    if verbose == True:
                        print("Adding "+str(len(stdict[ke]))+" tips to notallowed from subtree "+str(ke))
    ##                print (stdict[ke])
                    for eachtip in stdict[ke]:
                        valst[ke] = 0
                        if eachtip in notallowed:
                            pass
                        else:
                            notallowed.append(eachtip)
    ##                print("not allowed:"+str(len(notallowed)))
                else:
                    itis = "off"
        #create info for table
        if verbose == True:
            print("total clades: "+str(numberclades))
        goodstl = ""
        totstl = ""
        trees = "["+str(numberclades)+"]: "
        gatip = 0
        gtot = ""
        for thethings in notallowed:
            if string in thethings:
                gatip+=1
        #return if one clade found, too messy for strict para/mono with artifacts.
        if str(numberclades) == "1":
            if verbose == True:
                print("Not polyphyletic: only one clade discovered(clade cutoff="+str(dupcutoff)+")")
            cmono = CheckMonoSingleSubtree(cladeslist[0], stdict[cladeslist[0]], allstrtips)
            if cmono == "Monophyly":
                
                return(string, "Monophyly_Many_Exc", cladeslist[0], str(totst[cladeslist[0]]), str(len(allstrtips)-goodst[cladeslist[0]]), str(len(allstrtips)))
            else:
                return(string, "Paraphyly_Many_Exc", cladeslist[0], str(goodst[cladeslist[0]])+"/"+str(totst[cladeslist[0]]), str(len(allstrtips)-goodst[cladeslist[0]]), str(len(allstrtips)))
        else:
            for item in cladeslist:
                cmono = CheckMonoSingleSubtree(item, stdict[item], allstrtips)
                if cmono == "Monophyly":
                    gtot = gtot+str(goodst[item])+", "
                    trees = trees+item+"M, "
                else:
                    gtot = gtot+str(goodst[item])+"/"+str(totst[item])+", "
                    trees = trees+item+"P, "
            trees = trees[:-2]
            gtot = gtot[:-2]   
            if verbose == True:
                print("Polyphyletic")
            return(string, "Polyphyly", trees, gtot, str(len(allstrtips)-gatip), str(len(allstrtips)))

#unneeded?
def PerformScan(stringlist, treefile, artcutoff, fewnum, treesubfile, verbose = False):
    if stringlist[0] == "PerformGetRankFromSubtreesFile":
        ranknumstring = stringlist[1]
        
        ranknumlist = ranknumstring.split()
        perf = "yes"
    else:
        perf = "no"
    subtreefilename = treefile+"Subtrees.txt"
    #this bit will currently only work on my laptop (Artemis)
    #requires the following R script, which requires ape library installed:
###!/usr/bin/Rscript
###get arguments in a list, first nexus file then destination file
##myArgs <- commandArgs(trailingOnly = TRUE)
###print(myArgs)
##nexusfile<-myArgs[1]
##print("opening")
##print(nexusfile)
##subtreefile<-myArgs[2]
###create new empty destination file
##write("SUBTREES_BEGIN", file = subtreefile, append = FALSE)
###import necessary packages
##library(ape)
###open nexus file
##tree<-read.nexus(nexusfile)
##stree<-subtrees(tree)
##print("writing simple subtrees to file: ")
##print(subtreefile)
###write number of subtree followed by its tips, each on a newline.
##for(i in 1:length(stree)){
##  write(c("subtree#",i), file = subtreefile, append = TRUE)
##  write.tree(stree[[i]], file = subtreefile, append = TRUE)
##  for (n in 1:length(stree[[i]]$tip.label)){
##  write(stree[[i]]$tip.label[[n]], file = subtreefile, append = TRUE)
##  }}
##print("finished!")

    #with alias or symlink to keyword "subtreegen"
    if treesubfile == "NA":
        
        subtreefilename = treefile+"Subtrees.txt"
        subtreefilename = MakeSubtreesFile(treefile, subtreefilename, verbose)
    else:
        subtreefilename = treesubfile
    if perf == "yes":
        #preform subsampling
        stringlist = []
##        print(ranknumlist)
        for rnum in ranknumlist:
##            print(rnum)
            stringlist.append("RANKLEVEL: "+str(rnum))
            ranknum = int(rnum)
            rstringlist = GetRankFromSubtreesFile(subtreefilename, ranknum, verbose)
            for rstring in rstringlist:
                stringlist.append(rstring)
    stdict, st_trees = MakeSubtreesDict(subtreefilename, verbose)
    from prettytable import PrettyTable
    tablelist = []
    rankinfolist = []
    start = "yes"
    table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
    print("Checking Mono/Para/Polyphyly for each given group")
    #this string list is the ranks available.
    if len(ranknumlist) is 1:
        ranknum = ranknumlist[0]
    for string in stringlist:
        if "RANKLEVEL:" in string:
            table.sortby = "Topology"
            a = string.strip()
            ranknum = a[-1]
            if start == "yes":
                start = "no"
            else:
                tablelist.append(table)
                table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
            tablelist.append(string)
            continue
        
        if verbose == True:
            print("Checking group: "+string)
        CM = CheckMonophyly(string, stdict, fewnum, ranknum, verbose)
        tab = CM
        if CM == "NotMono":
            if verbose == True:
                print("Not monophyletic")
            CI = CheckInclusions(string, stdict, ranknum, "no", verbose)
            tab = CI
            if CI[0] == "NotParaAlone":
                if verbose == True:
                    print("Not paraphyletic")
                CA = CheckArtifacts(string, stdict, CI, artcutoff, verbose)
                tab = CA
        if CM == "few":
            CI = CheckInclusions(string, stdict, ranknum, "few", verbose)
            CA = CheckArtifacts(string, stdict, CI, artcutoff, "few", verbose)
            tab = CA

        table.add_row(tab)
        #whichever is good, print(?) or just return I guess if we're gunna run though lots
    table.sortby = "Topology"
    tablelist.append(table)
   
##    for titem in tablelist:
##        print titem
    return tablelist

########END OF ORIGINAL PROGRAM#######


#unneeded
def GetRankFromSubtreesFile_SubSample(subtrees, ranknum, groupstring, verbose = False):
    print("Getting ranks from subtrees file seqIDS")
    ranklist = []
    lines = []
    new = "no"
    newtree = "no"
    with open(subtrees) as old:
        for line in old:
            if line == "SUBTREES_BEGIN\n":
                print("beginning...")
            elif new == "yes":
                new = "no"
                newtree = "yes"
            elif newtree == "yes":
                newtree = "no"
            elif line == "subtree#\n":
                new = "yes"
            else:
                if "\n" in line:
                    line = line[:-1]
                
                if line in lines:
                    pass
                else:
                    lines.append(line)
    for example in lines:
        if groupstring in example:
            pass
        else:
            continue
        import re     
        gitax = re.sub("(.*)(\|)(gi#?\|?)([0-9]*)(.*)", "\\4~\\1", example)
##  this bit only works if your seqID format is taxon|omic|info|gi#|whateverelse
##if switch to accession nums, replace gi(above) with acn or something.
        
        gi, tax = gitax.split("~")
        
        tn = tax[:-1]
        taxlist = tn.split("|")
        
        therankstr = taxlist[ranknum-1]
        therankstr = therankstr+"|"
        if therankstr is "Proteobac|":
            #debugging prints
            print ("PROTEOBAC FOUND")
        if therankstr in ranklist:
            pass
        else:
            ranklist.append(therankstr)
##        print(ranklist)
##    print("Added "+str(len(ranklist))+" groupnames (subset of "+groupstring+" to parsing queue")
    print (ranklist)
    return ranklist

#unneeded
def PerformScan_SubSample(stringlist, treefile, artcutoff, fewnum, treesubfile, verbose = False):
    print("Beginning outer clade-determination")
    if stringlist[0] == "PerformGetRankFromSubtreesFile":
        ranknumstring = stringlist[1]
        
        ranknumlist = ranknumstring.split()
        perf = "yes"
    else:
        perf = "no"
    subtreefilename = treefile+"Subtrees.txt"
    #this bit will currently only work on my laptop (Artemis)
    if treesubfile == "NA":
        subtreefilename = treefile+"Subtrees.txt"
        subtreefilename = MakeSubtreesFile(treefile, subtreefilename, verbose)
    else:
        subtreefilename = treesubfile
    if perf == "yes":
        stringlist = []
##        print(ranknumlist)
        for rnum in ranknumlist:
##            print(rnum)
            stringlist.append("RANKLEVEL: "+str(rnum))
            ranknum = int(rnum)
            rstringlist = GetRankFromSubtreesFile(subtreefilename, ranknum, verbose)
            for rstring in rstringlist:
                stringlist.append(rstring)
    stdict, st_trees = MakeSubtreesDict(subtreefilename, verbose)
    from prettytable import PrettyTable
    tablelist = []
    rankinfolist = []
    start = "yes"
    table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
    subtrees_new_dict = {}
##    print("Checking Mono/Para/Polyphyly for each given group")
    for string in stringlist:
        if "RANKLEVEL:" in string:
            table.sortby = "Topology"
            
            if start == "yes":
                start = "no"
            else:
                tablelist.append(table)
                table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
            tablelist.append(string)
            continue
        
##        print("Checking group: "+string)
        CM = CheckMonophyly(string, stdict, fewnum, verbose)
        tab = CM
        if CM == "NotMono":
##            print("Not monophyletic")
            CI = CheckInclusions(string, stdict, verbose)
            tab = CI
            if CI[0] == "NotParaAlone":
##                print("Not paraphyletic")
                CA = CheckArtifacts(string, stdict, CI, artcutoff, verbose)
                tab = CA
        if CM == "few":
            CI = CheckInclusions(string, stdict, "few")
            CA = CheckArtifacts(string, stdict, CI, artcutoff, "few", verbose)
            tab = CA

        table.add_row(tab)
        subtrees_line = tab[2]
#      
##        print(subtrees_line)
        subtrees_line_list = subtrees_line.split(",")
##        print(subtrees_line_list)
        for item in subtrees_line_list:
            if "[" in item:
#                item = item[3:]
                nitem = re.sub("[A-Za-z \[\]:]*", "", item)
            if len(nitem) >4:
                print (nitem+"   "+item)
##            print(nitem)
            if nitem in subtrees_new_dict:
#                 subtrees_new_dict[nitem] =  subtrees_new_dict[nitem]+" "+string
                subtrees_new_dict[nitem] = string
            #subtrees new dict contains 1,2,3 (subtree numbers) : actinobacteria (or "dogs cat rats")?
    
        #whichever is good, print(?) or just return I guess if we're gunna run though lots
    table.sortby = "Topology"
    tablelist.append(table)
   
##    for titem in tablelist:
##        print titem
##    print("SubNewDict")
####    print(subtrees_new_dict)
    return tablelist, subtrees_new_dict, st_trees

#unneeded
def PerformScan_SubSample_Inner(stringlist, treefile, artcutoff, fewnum, treesubfile, groupstring, tree_in_text, verbose = False):
    
    subtrees_new_dict = {}
    if stringlist[0] == "PerformGetRankFromSubtreesFile":
        ranknumstring = stringlist[1]
        ranknumlist = ranknumstring.split()
        newranknumlist = []
        for item in ranknumlist:
            item = int(item) + 1
            newranknumlist.append(item)
        ranknumlist = newranknumlist
        perf = "yes"
    else:
        perf = "no"
    subtreefilename = treefile+"Subtrees.txt"
    #this bit will currently only work on my laptop (Artemis)
    if treesubfile == "NA":
        subtreefilename = treefile+"Subtrees.txt"
        subtreefilename = MakeSubtreesFile(treefile, subtreefilename, verbose)
    else:
        subtreefilename = treesubfile
    if perf == "yes":
        stringlist = []
##        print(ranknumlist)
        for rnum in ranknumlist:
##            print(rnum)
            stringlist.append("RANKLEVEL: "+str(rnum))
            ranknum = int(rnum)
            rstringlist = GetRankFromSubtreesFile_SubSample(subtreefilename, ranknum, groupstring, verbose)
            for rstring in rstringlist:
                stringlist.append(rstring)
    print(groupstring)
    #stdict is ??
    #st_trees is number:newickstr
    stdict, st_trees = MakeSubtreesDict(subtreefilename, verbose)
    from prettytable import PrettyTable
    tablelist = []
    rankinfolist = []
    start = "yes"
    table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
    for string in stringlist:
        print("checking :"+string)
        if "RANKLEVEL:" in string:
            table.sortby = "Topology"
            if start == "yes":
                start = "no"
            else:
                tablelist.append(table)
                table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
            tablelist.append(string)
            continue
##        print("Checking group: "+string)
        CM = CheckMonophyly(string, stdict, fewnum, verbose)
        tab = CM
        if CM == "NotMono":
##            print("Not monophyletic")
            CI = CheckInclusions(string, stdict, verbose)
            tab = CI
            if CI[0] == "NotParaAlone":
##                print("Not paraphyletic")
                CA = CheckArtifacts(string, stdict, CI, artcutoff, verbose)
                tab = CA
        if CM == "few":
            CI = CheckInclusions(string, stdict, "few")
            CA = CheckArtifacts(string, stdict, CI, artcutoff, "few")
            tab = CA

        table.add_row(tab)
        subtrees_line = tab[2]
    
        subtrees_line_list = subtrees_line.split(",")
        for item in subtrees_line_list:
            if "[" in item:
                item = item[3:]
            else:     
                nitem = re.sub("[A-Za-z :\[\]]*", "", item)
                if nitem in subtrees_new_dict:
                   
                    subtrees_new_dict[nitem] =  subtrees_new_dict[nitem]+" "+string
                subtrees_new_dict[nitem] = string
        print("Inner subtree determined to be:"+nitem)
        
        #whichever is good, print(?) or just return I guess if we're gunna run though lots
    table.sortby = "Topology"
    tablelist.append(table)
   
##    for titem in tablelist:
##        print titem
    #if only one order found in class-group-search (eg, all Halos are Natrialbales) we should return two halo sequences, one basal and one random.
    #else we only return basal seq of all orders represented within class-clade.
    if len(subtrees_new_dict) == 1:
        subtrees_new_dict["alone"] = string
        st_trees["alone"] = tree_in_text
    return tablelist, subtrees_new_dict, st_trees

#?? unneeded i think
def LeastIndentTip(newicktree, alone, groupstring, verbose=False):
    indent = 0
    newtip = "yes"
    tipend = "no"
    ongoingtip = "no"
    tipdict = {}
    for item in newicktree:
        if item == "(":
            indent+=1
        elif item == ")":
            indent-=1
        elif item == ",":
            pass
        elif item == ":":
            ongoingtip = "no"
            if groupstring in currenttip:
                tipdict[currenttip] = indent
            else:
                pass
        elif ongoingtip == "no":
            if item == "A":
                newtip = "yes"
                currenttip = item
                ongoingtip = "yes"
            elif item == "B":
                newtip = "yes"
                currenttip = item
                ongoingtip = "yes"
            elif item == "E":
                newtip = "yes"
                currenttip = item
                ongoingtip = "yes"
            elif item == "X":
                newtip = "yes"
                currenttip = item
                ongoingtip = "yes"
            else:
                pass
        elif ongoingtip == "yes":
            currenttip = currenttip+item
        else:
            print("Error in indentcounting")
    if len(tipdict) == 0:
        print("Zero!"+groupstring)
        print(newicktree)
        return("none")
    if alone == "yes":
        import random
##        print (newicktree)
##        print (tipdict)
        shallowtip = random.choice(tipdict.keys())
    else:  
        shallowtip, value = min(tipdict.iteritems(), key=lambda x:x[1])
    if groupstring in shallowtip:
            return shallowtip
    else:
        print("groupstring not in chosen tip")
        raise SystemExit
                
                
        
##________________
##EASYMODE: subsampling
##0. PreformScan_SubSample -> subtrees_new_list (of all subtrees found at rank (class)
##1. for each item in subtrees_new_list, get subtrees_tree and run easymode:


#main thing.
def SubSampling_Master(stringlist, treefile, artcutoff, fewnum, treesubfile, inputfasta, verbose):
    tablelist, subtrees_new_dict, st_trees = PerformScan_SubSample(stringlist, treefile, artcutoff, fewnum, treesubfile, verbose)
    chosen_seqs = []
##    print(subtrees_new_dict)
    import random
    print("Clade - dictionaries were successfully created. Moving on to inner SS... ")
    with open(treefile+"SubSamp_INFO.txt", "w") as info:
        for item in subtrees_new_dict:
            groupstring = subtrees_new_dict[item]
            a = st_trees[item]
            with open ("Temp"+item+".txt", "w") as new:
                 new.write(a)
            dire = os.getcwd()
            os.system("ConvertPhylo "+dire+" Temp"+item+".txt newick nexus")
            nexus = "Temp"+item+".nexus"
            #run the thing on the nexus file with given string(s)
            inner_treefile=nexus
            inner_treesubfile="NA"
            inner_stringlist=stringlist
            
            info.write("\nfor group "+groupstring+" subtree "+item+" the following sequences have been subsampled:")
            inner_tablelist, inner_subtrees_new_dict, inner_st_trees = PerformScan_SubSample_Inner(inner_stringlist, inner_treefile, artcutoff, fewnum, inner_treesubfile, groupstring, a, verbose)
##            if verbose == True:
            print("Created inner clade - dictionaries for group: "+groupstring+" consisting of "+str(len(inner_subtrees_new_dict))+" subtrees")
            for tree in inner_subtrees_new_dict:
                inner_string = inner_subtrees_new_dict[tree]
                if tree == "alone":
                    info.write("\nOnly one subgroup found... keeping an additional random sequence.")
                    chosen_sequence = LeastIndentTip(inner_st_trees[tree], "yes", groupstring, verbose)
                else:
                    chosen_sequence = LeastIndentTip(inner_st_trees[tree], "no", groupstring, verbose)
                if chosen_sequence in chosen_seqs:
                    
                    info.write("\n"+inner_string+": "+chosen_sequence+" is a Duplicate... avoiding")
                elif chosen_sequence == "none":
                    pass
                else:
                    chosen_seqs.append(chosen_sequence)
                    info.write("\n"+inner_string+": "+chosen_sequence)
            #remove generated files
            os.system("rm "+nexus)
            os.system("rm Temp"+item+".txt")        
            subtreefilename = nexus+"Subtrees.txt"
            os.system("rm "+subtreefilename)
        info.write("\nFinished chosing extracted sequences will be at: "+inputfasta+"Extracted.fasta\nInfo is available at: "+treefile+"SubSamp_INFO.txt")

    print("Finished choosing sequences to keep.")
    print("Extracting given sequences.... ")
    keepseqs = ""
    for seqid in chosen_seqs:
        keepseqs = keepseqs+" "+seqid
    keepseqs=keepseqs[:-1]
    print(inputfasta)
    print(keepseqs)
    os.system("feast . "+inputfasta+" -ex \""+keepseqs+"\"")
    print("Finished, extracted sequences SHOULD be at: "+inputfasta+"Extracted.fasta")
    print("Info available at: "+treefile+"SubSampInfo.txt")






#THIS IS WHAT OVERALL DTL CALLS
def PerformScan_Complex_SubSample(stringlist, treefile, artcutoff, fewnum, treesubfile, verbose = False):
    print("Beginning outer clade-determination")
    typedict = {}
    if stringlist[0] == "PerformGetRankFromSubtreesFile":
        ranknumstring = stringlist[1]
        ranknumlist = ranknumstring.split()
        perf = "yes"
    else:
        perf = "no"
    subtreefilename = treefile+"Subtrees.txt"
    print(treesubfile)
    print(ranknumlist)
    #this bit will currently only work on my laptop (Artemis)
    if treesubfile == "NA":
        subtreefilename = treefile+"Subtrees.txt"
        subtreefilename = MakeSubtreesFile(treefile, subtreefilename, verbose)
    else:
        subtreefilename = treesubfile
    if perf == "yes":
        stringlist = []
##        print(ranknumlist)
        for rnum in ranknumlist:
##            print(rnum)
            stringlist.append("RANKLEVEL: "+str(rnum))
            ranknum = int(rnum)
            rstringlist = GetRankFromSubtreesFile(subtreefilename, ranknum, verbose)
            for rstring in rstringlist:
                stringlist.append(rstring)
            #stringlist will contain all the strings found in the right position of the sequence IDS.
    stdict, st_trees = MakeSubtreesDict(subtreefilename, verbose)
    
    #this is a dict name: tips and name: tree
    from prettytable import PrettyTable
    tablelist = []
    rankinfolist = []
    start = "yes"

    #initialized the table
    
    table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
    subtrees_new_dict = {}
##    print("Checking Mono/Para/Polyphyly for each given group")
    if len(ranknumlist) == 1:
        if type(ranknumlist) == list:
            ranknum = int(ranknumlist[0])
        else:
            ranknum = int(ranknumlist)
    for string in stringlist:
        #if ranklevel is stuck in there, we are going to print multiple tables, so start a new onw. this is not implemented in the OverALLDTL_Detector analysis.
        if "RANKLEVEL:" in string:
            table.sortby = "Topology"
            string = string.strip()
            ranknum = string[-1]
            if start == "yes":
                start = "no"
            else:
                tablelist.append(table)
                table = PrettyTable(['Group', 'Topology', 'Subtree#', 'ST_Tips', 'Excluded', "Total"])
            tablelist.append(string)
            continue
##        print("Checking group: "+string)
        #check monophyly for this specific string
        CM = CheckMonophyly(string, stdict, fewnum, ranknum, verbose)
        tab = CM
        if CM == "NotMono":
##            print("Not monophyletic")
            CI = CheckInclusions(string, stdict, ranknum, "no", verbose)
            tab = CI
            if CI[0] == "NotParaAlone":
##                print("Not paraphyletic")
                CA = CheckArtifacts(string, stdict, CI, artcutoff, verbose)
                tab = CA
        if CM == "few":
            CI = CheckInclusions(string, stdict, ranknum, "few", verbose)
            CA = CheckArtifacts(string, stdict, CI, artcutoff, "few", verbose)
            tab = CA
        #add a row to the table.
        table.add_row(tab)
        subtrees_line = tab[2]
        type_line = tab[1]
        
##        print(subtrees_line)
        subtrees_line_list = subtrees_line.split(",")
##        print(subtrees_line_list)
        for item in subtrees_line_list:
            if "[" in item:
                item = item[3:]
            nitem = re.sub("[A-Za-z \[\]:]*", "", item)
            if len(nitem) >4:
                print (nitem+"   "+item)
##            print(nitem)
            if nitem in subtrees_new_dict:
                 subtrees_new_dict[nitem] =  subtrees_new_dict[nitem]+" "+string
            subtrees_new_dict[nitem] = string
            typedict[nitem] = type_line
            #subtrees new dict contains 1,2,3 (subtree numbers) : actinobacteria (or "dogs cat rats")?
    
        #whichever is good, print(?) or just return I guess if we're gunna run though lots
    table.sortby = "Topology"
    
    tablelist.append(table)
    return tablelist, subtrees_new_dict, st_trees, stdict, typedict



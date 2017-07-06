#!/usr/bin/python

#abigailc@Actaeon Jan 4 2017
#this is a module to control writing and reading progress files - for use when a run stops halfway.

class Pickle_Project_Storage:
    def __init__(self, name):
        self.name = name
        self.short_clades = []
        self.big_clades = []
        self.original_fasta = ""
        self.projectname = ""

def Pickle_Progress(list_big, list_short, original_fasta, projectname):
    import pickle
    with open("Pickle_"+projectname+".bin", 'wb') as f:
        w1 = Pickle_Project_Storage(projectname)
        w1.short_clades = list_short
        w1.big_clades = list_big
        w1.original_fasta = original_fasta
        w1.projectname = w1.name
        w1.picklefile = "Pickle_"+projectname+".bin"
        pickle.dump(w1, f)
    print("saved pickle to :"+ w1.picklefile)
    return w1.picklefile

def Pickle_Read(picklefile):
    with open(picklefile, 'rb') as f:
        w1_restore = pickle.load(f)
    return w1_restore


def Write_Progress_File(list_big, list_short, projectname):
    with open("./"+projectname+"/"+projectname+"PROGRESS.txt", "w") as new:
        new.write(">big\n")
        for subtree in list_big:
            new.write(">>tree\n")
            new.write(subtree.ret_number()+"\n")
            new.write(subtree.ret_string()+"\n")
            new.write(subtree.ret_fasta()+"\n")
            new.write(subtree.ret_gene_tree()+"\n")
            new.write(subtree.ret_species_tree()+"\n")
        new.write(">small\n")
        for subtree in list_big:
            new.write(">>tree\n")
            new.write(subtree.ret_number()+"\n")
            new.write(subtree.ret_string()+"\n")
            new.write(subtree.ret_fasta()+"\n")
            
def Read_Progress_File(input_file, projectname):
    big_list = []
    small_list = []
    with open(input_file) as old:
        big = False
        small = False
        for line in old:
            if ">big\n" == line:
                big = True
                continue
            elif ">small\n" == line:
                big = False
                small = True
                continue
            else:
                if ">>tree\n" == line:
                    count = 0
                    
                else:
                    count +=1
                    if count == 1:
                        s = Subtree(line.strip())
                        #if big is True:
                         #   list_big.append(s)
                        #if small is True:
                         #   list_small.append(s)
                    if count == 2:
                        s.set_string(line.strip())
                    if count == 3:
                        pre = line.strip()
                        pre2 = pre[:-11]
                        s.prefix = pre2
                        a = os.path.isfile(line.strip())
                        if a is True:
                            s.set_fasta(line.strip())
                            s.set_fasta_object()
                        else:
                            b = os.path.isfile(projectname+"/Gene_Trees/fasta/"+line.strip())
                            if b is True:
                                s.set_fasta(projectname+"/Gene_Trees/fasta/"+line.strip())
                                s.set_fasta_object()
                            else:
                                print("Could not load fasta for prefix: "+s.prefix)
                    if count == 4:
                        if line.strip() == "":
                            #look for the file
                            e = os.path.isfile(projectname+"/Gene_Trees/trees/RAxML_bipartitions."+s.prefix+"_gene")
                            if e is True:
                                with open(projectname+"/Gene_Trees/trees/RAxML_bipartitions."+s.prefix+"_gene") as gene_tree_fi:
                                    gene_tree_str=gene_tree_fi.read()
                                s.gene_tree = gene_tree_str
                                s.set_gene_tree_species_tips()
                            else:
                                print("Could not load gene tree for prefix:"+s.prefix)
                        else:
                            s.gene_tree = line.strip()
                            s.set_gene_tree_species_tips()
                    if count == 5:
                        if line.strip() == "":
                            #look for the file
                            c = os.path.isfile(projectname+"/Species_Trees/trees/"+s.prefix+"_CC_Rax_Bipart.newick")
                            if c is True:
                                with open(projectname+"/Species_Trees/trees/"+s.prefix+"_CC_Rax_Bipart.newick") as spec_tree_fi:
                                    spec_tree_str=spec_tree_fi.read()
                                s.species_tree = spec_tree_str
                            else:
                                print("Could not load species tree for prefix:"+s.prefix)
                        else:
                            s.species_tree = line.strip()
                        if big is True:
                            big_list.append(s)
                        elif small is True:
                            small_list.append(s)
                #this needs to set dependent on small vs big
      
                    
                #here?
    return big_list, small_list

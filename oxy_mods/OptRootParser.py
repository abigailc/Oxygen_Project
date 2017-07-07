#!/usr/bin/python

#abigailc@Actaeon Jan 23 2017

#this friend will parse the output of optroot and extract the best tree, which it will write to a file of defined name.

#useage OptRootParser.py -i INPUTNAME -o OUTPUTNAME




import argparse
#parser

if __name__ == "__main__":

 
    parser = argparse.ArgumentParser(description="All")
    parser.add_argument("-i", "--input", action = "store", help="provide input file")
    parser.add_argument("-o", "--output", action = "store", help="provide output file name")
    parser.add_argument("-s", "--species", action="store", help="provide the name of the file input to opt root (contains species tree in line 1")

    args = parser.parse_args()

    with open(args.input) as infile:
        print("Opening file: "+ args.input+" to read trees from...")
        tree = ""
        for line in infile:
            if line[0] == "(":
                if tree == "":
                    tree = line
                else:
                    break
    with open(args.species) as spinfile:
        print("opening file: "+args.species+" as the file containing a species tree in line 1")
        species_tree = ""
        i = 1
        for line in spinfile:
            if i is 1:
                species_tree = line
            else:
                pass
            i+=1
    with open(args.output, "w") as outfile:
        outfile.write(species_tree)
        outfile.write("[&R]"+tree)
    print("trees were written to the file: "+ args.output)
    

        ####This is not getting the correct species tree associated with it - just  writing the new gene tree w/ optimal rooting. why am i not combining these things?
        

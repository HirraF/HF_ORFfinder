import argparse

def comprev(genome):
    """This function takes the inputted genome sequence
	and creates the complementray reverse sequence.
	Both the orignal sequence and the complementary 
    reverse get stored in a list which is returned """
    #list for seq and its rev complement
    seq_and_comprev=[]
    seq_and_comprev.append(genome)
    #create the reverse complement
    transtab=str.maketrans('ATCG','TAGC')
    compsq=genome.translate(transtab)
    comprsq=compsq[::-1]
    seq_and_comprev.append(comprsq)
    #final list contains 2 item, sequence and its reverse complement
    return seq_and_comprev

#dictionary which allow codon to be translated into amino acid 
geneticcode={
    "AAA": "K","AAC":"N", "AAG":"K", 
    "AAT":"N","ACA":"T","ACC":"T", 
    "ACG":"T","ACT":"T","AGA":"R", 
    "AGC":"S", "AGG":"R","AGT":"S",
    "ATA":"I", "ATC":"I", "ATG":"M", 
    "ATT":"I","CAA":"Q", "CAC":"H", 
    "CAG":"Q","CAT":"H","CCA":"P", 
    "CCC":"P", "CCG":"P", "CCT":"P",
    "CGA":"R", "CGC":"R", "CGG":"R", 
    "CGT":"R","CTA":"L","CTC":"L", 
    "CTG":"L", "CTT":"L","GAA":"E", 
    "GAC":"D", "GAG":"E", "GAT":"D",
    "GCA":"A", "GCC":"A", "GCG":"A", 
    "GCT":"A","GGA":"G", "GGC":"G", 
    "GGG":"G", "GGT":"G","GTA":"V", 
    "GTC":"V", "GTG":"V", "GTT":"V",
    "TAA":"", "TAC":"Y", "TAG":"", 
    "TAT":"Y","TCA":"S", "TCC":"S", 
    "TCG":"S", "TCT":"S","TGA":"", 
    "TGC":"C", "TGG":"W", "TGT":"C",
    "TTA":"L", "TTC":"F", "TTG":"L", 
    "TTT":"F", 
}

#function for translation nt codon into amino acid 
def aatranslation(orf):
    """This function takes the ORF string and translates 
    the codons into Amino Acids according to the dictionary genetic code. """
    peptide=''
    x=0
    for i in range(0,int((len(orf)/3))): 
        codon=orf[x:x+3]
    #translate the codon from the genetic code dict
    #if codon not found in dict return ? in its place
        aa=geneticcode.get(codon, '?')
        peptide=peptide+aa
        x=x+3
    return peptide

def searchinsq(seq):
    """This function takes the sequences and searches for ORF. ORF string 
    then translated and printed to ouptut file or terminal 
    """
    global total                            #total number of orfs found  will be stored in this 
    frame= i+1                              #what frame we are in
    orfcounter=0                            #what orf of the frame is this
    w=i                                     #counts the first nt of each codon 
    adding=False
    while w<len(seq):
        currentcodon=seq[w:w+3]             #select codon
        if currentcodon =='ATG':
            ORFstring=''                     #empty string ready to recieve the orf we have found
            ORFstring=ORFstring+currentcodon #append start codon to orf 
            adding=True
            while adding==True:
                w=w+3
                nextcodon=seq[w:w+3]         # move down sequence by once codon
                if w >len(seq):
                    adding=False
                if nextcodon!='TGA'and nextcodon!='TAA' and nextcodon!='TAG':  
                        ORFstring=ORFstring+nextcodon                 #if next codon not a stop, append to orfstring
                else:
                    ORFstring=ORFstring+nextcodon
                    adding=False                                      #else if stop, stop making orfstring
                    #checks on the ORF: no N and a minimum length
                    if 'N' not in ORFstring and len(ORFstring)>= args.minlen:
                        orfcounter=orfcounter+1
                        ntlen=len(ORFstring)
                        #call translation function to translate produced orf into peptide 
                        mypeptide=aatranslation(ORFstring)
                        total=total+1                         #add 1 to total orf found so far
                        #write unique name 
                        uniquename='\n{head}_{direc}{fnumb}_{nt}_{orfn}\n'.format(head=header,direc=direction,fnumb=frame, nt=ntlen, orfn=str(orfcounter))
                        #writes to file if specified otherwise print to terminal
                        if args.outputfile != None:
                            newf.write(uniquename)
                            newf.write(mypeptide)
                        else:
                            print(uniquename,mypeptide)
                            
        else:
            w=w+3 #if codon not a start, shift down sequnece by one codon. continue searching for the next start
    return total

#parser
parser=argparse.ArgumentParser(description='This program searches an input file for ORFs')
parser.add_argument('-f','--framenumber', type= int, choices=[1,2,3,-3,-2,-1],help='Give the frame you would like to search in, default will search all 6')
parser.add_argument('-i','--inputfile', required=True ,help='Give the name of the file you wish to input', nargs='*')
parser.add_argument('-o','--outputfile', help='Give the name of the file you wish to output to')
parser.add_argument('-n', '--minlen', type=int, default= 75, help='Specify the minimum length (in nt) of orf to return. Default value is 75')
args=parser.parse_args()

#input file list- iterate through 
for files in args.inputfile:
    total=0
    #take and save users input file. then open
    genome=files
    f=open(genome,'r')

    #extract header from fasta file and concatenate rest of genome into one string 
    genomestr=''
    header=''
    for line in f:
        if line.startswith('>'):
            header=line.split()[0]
        else:                                 
            lineadd=line.rstrip("\n")
            uppercase=lineadd.upper()
            genomestr=genomestr + uppercase

    #checking format of file inputted in. raise exception if doesnt look like genome seq.
    if 'A' and 'T' and 'C' and 'G' not in genomestr:
        raise Exception ('File {} inputted does not look like a genome seq. Submit a different sequence and try again'.format(files))

    #no header on file- makes a header
    if header == '':
        header='supplied_genome'

    #open output file if given
    if args.outputfile!= None:
        newf = open(args.outputfile, "w")

    #generate list: sequence and rev sequence
    mylist=comprev(genomestr)
    directions=['+','-']   
    #direction defined as + for provided sequence and - for comp reverse  


    #if frame number specified:
        #first: select correct direction.
        #second:select correct reading frame
    #else: no frame seleted by user so all 6 searched.

    if args.framenumber != None:
        if int(args.framenumber) > 0:
            direction=directions[0]
            sequenceused = mylist[0]
        else:
            direction=directions[1]
            sequenceused=mylist[1]
        start_pos=abs(int(args.framenumber))-1
        i=start_pos
        searchinsq(sequenceused)
    else:
        for i in range(0,3):
            for s in range(0,2):
                direction=directions[s]
                searchinsq(mylist[s])

    #prints the total number of orfs foudn in the search just done 
    print('total number of orfs found: ', total)









"""
Eddie Cho
BIOL 668
SP22
Python OOP Project - Improvised/Improved 11/5/24

Using OOP concepts, it models biological sequences such as DNA, RNA, and Protein.
The script utilizes classes to represent different types of sequences and provides
methods (functions) to manipulate and analyze these sequences.


"""
##Step-1: Initialize values

#import 're' module used for regular experssions. Helpful for patter mateching and replacing strings in sequences.
import re

#Create a dictionary called 'standard_code' provided during the lecture, which maps RNA codons (three-letter sequences of nucleotides) to their corresponding amino acirds. These indicate how to translate RNA sequences into proteins.
standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

#Create a dictionary called 'kyte_doolittle', which assigns a hydrophobicity value to each amino acid. 
#Hydrophobicity is a measure of how much a molecule "hides" from water.
#Provided during the lecture
kyte_doolittle={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3}

#Create a dictionary called aa_mol_weights containing the molecular weights of different amino acids.
#It is used to calculate the weight of a protein based on its amino acid sequence.
aa_mol_weights={'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
                'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
                'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
                'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19}

#Step-2: Create Seq class
class Seq:

#Create constructor. Automatically called when a new instance (object) of a class is created.
#Purpose: to initialize the instance's attributes and set up any necessary initial state for the object.
#Note: 'self' referes to the instance. In other words, it is the instance of the object on which the methos is being called.

    def __init__(self, sequence, gene, species):
        self.sequence = sequence  #initializes the object with values for sequence (the actual genetic sequence)
        self.gene = gene #initializes gene (the name of the gene)
        self.species = species #species (the species of the organism)
        self.kmers = [] #Also initializes an empty list kmers to hold subsequences of a specified length.

#A special method that returns a string representation of the object when printing it
    def __str__(self):
        self.sequence = self.sequence.strip()  #Strip any leading/trailing spaces from the sequence.
        self.sequence = self.sequence.upper()  #Convert the sequence to uppercase
        return f"Species: {self.species}, Gene: {self.gene}, Sequence: {self.sequence}" #Returns a string

#This method prints out the species, gene, and sequence in a readable format.
    def print_record(self):
        print(self.species + " " + self.gene + ": " + self.sequence)

#This method generates 'k-mers', subseqeunces of length k from the sequence and stores them in the kmers list.
#Uses for loop to iterate over the sequence and extract every possible subsequence of length k.
#If the subsequence is not already in the kmers list, it appends it.
#In other words, 1.Iteration, 2.Extracting a k-mer 3.Checking uniqueness, 4.Storing k-mers 

    def make_kmers(self, k=3): #self refers to the current object (instance of the class) that the method is called on. 
    
    #for statement below: it will iterate through the range of indices of the sequence to extract a valid k-mer.
    #The variable i is the starting index of each k-mer to extract from the sequence.
    #'len(self.sequence)' returns the length of the sequence (# of nucleotides in it).
        for i in range(len(self.sequence) - (k - 1)): #since k= length of the k-mer, subtract k-1. Last index = i+k
           
            #if statement below: Slides the sequence starting at index i and extracts the next k
            #'not in self.kmers': checks if the k-mer that was slide is NOT already in teh list of k-mers (self.kmers).
            #The idea is that we only want to store unique k-mers in the list. If in the list, skip adding it again.
            if self.sequence[i:i + k] not in self.kmers:
             	
             	#Append if not list (since we want only unique k-mers)
                self.kmers.append(self.sequence[i:i + k]) 

#Generates a FASTA format string with a header line (>) followed by the sequence on the next line
    def fasta(self):
        self.fasta_str = '> ' + self.species + ' ' + self.gene + '\n' + self.sequence
        return self.fasta_str

#Step-3: Create DNA class (must inherit Seq class)

#Inheriting from Seq class

class DNA(Seq):

#**kwargs stands for keyword arguments to pass # of keyword arguments as a dictionary
#The keys of this dictionary are the argument names, and the values are the corresponding argument values
    def __init__(self, sequence, gene, species, geneid, **kwargs): 
        
        #call the parent class's __init__ method with super()
        super().__init__(sequence, gene, species)
        
        #Initialize the geneid attribute 
        self.geneid = geneid  

#Overrides the __str__ method for DNA objects to include geneid as well
    def __str__(self):
        return f"Gene: {self.gene}, Species: {self.species}, Gene ID: {self.geneid}, Sequence: {self.sequence}"

#Calculates the GC content of the DNA sequence by counting the occurances of 'G' and 'C' using regular experssions
    def analysis(self):
        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence))
        return gc

#Returns the reverse complement of the DNA sequence. It replaces nucleotide with its complement then reverses seq. order.
#In oter words, 1.Reverse sequence 2.Complement each base 3.Build the reverse complement 4.Return
    def reverse_complement(self):
        
        #Reverse complement. Unknown remains as N.
        complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', "N":"N"} 
        
        #An empty string is created to hold seequence as it processes
        self.reverse_complement_sequence = '' 
        
        #Iterates over each nucleotide in the reversed manner
        for j in reversed(self.sequence): 
            
            #for each reversed j, looks up complement of J from the dictionary
            #Then, complement is appended using '+='
            self.reverse_complement_sequence+= complement_dict[j] 
            
            #After the loop finishes, it will contain reverse complment of the original sequence and returns
        return self.reverse_complement_sequence

#This method prints out the species, gene, and sequence in a readable format.
    def print_info(self):
        print(self.geneid + " " +self.species+' '+self.gene+' ' +self.sequence)

#Generates the six reading frames of the RNA sequence, which includes three frames from the forward sequence
#and three from the reverse complement sequence. Commonly used to analyze possible coding sequences
#useful for analyzing different possible open reading frames (ORFs) in a sequence
#Since some coding sequences may only be apparent when read in a specific reading frame
#either in the forward and reverse direction

    def six_frames(self):
        self.SixFrames = [] #Initalize an empty list to hold the six reading frames: 3 from the foward and 3 from reverse
	
	#Starts the first loop to generate the frames from the foward sequence
	#range(3) means the loop will iterate 3 times, i taking values 0, 1, 2
	#representing the different starting positions for the frames in the forward sequence
	#by reading the sequence from the left to right       
        for i in range(3): 
        
        #with each iteration, append a substring of the original sequence starting at position i to the 'SixFrames' list
        #For example, when i=0, self.sequence[0:] takes the whole sequence from the first nucleotide.
        #when i=1, self.sequence[1:] starts from the second nucleotide onward.
        #when i=2, self.sequence[2:] starts from the third nucleotide onward. 
            self.SixFrames.append(self.sequence[i:]) 
            
            #Used to generate the next three frames from the reverse complement
            self.reverse_complement_sequence = self.reverse_complement()
        
        #Starts the second loop to generate frames from the reverse complement sequence
        #Just like the first loop, range(3) will iterate 3 times (i= 0, 1, 2)
        #Append substrings of the reverse complement sequence to the 'SixFrames' list starting i=0, 1, 2
        #For example, when i=0, self.reverse_complement_sequence[0:] takes the whole reverse complement sequence.
        #when i=1, self.reverse_complement_sequence[1:] starts from the second nucleotide of the reverse complement 
        #when i=2, self.reverse_complement_sequence[2:] starts from the third nucleotide of the reverse complement     
        for i in range(3):
            self.SixFrames.append(self.reverse_complement_sequence[i:])
        
        #Finally, after both loops are finished, the method returns the 'SixFrames' list containing six sequences:
        #first three frames are derived from the forward sequence.
        #next three frames are derived from the reverse complement sequence.  
        return self.SixFrames

#Step-4: Create RNA class(Must inherit DNA class)
class RNA(DNA):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species,geneid)
        self.__str__()
        self.codons = []

#Replaces the sequence by replacing all 'T's with 'U's to represent RNA seq
    def __str__(self):
        self.sequence = re.sub('T', 'U', self.sequence)
        return self.sequence

#Divides the RNA sequence into codon (three letters of nucleotides) and stores them in the codons list.
    def make_codons(self):
    
    	#Starting at index 0, the beginning of the RNA seq then stops at the length of the sequence -2
    	#In order to ensure enough characters left to make a complete codon and avoid going out of bounds (3).
    	#The loop increments by 3 on each iteration, so each will process a block of 3 letters at a time.
        for i in range (0,len(self.sequence) - 2, 3):
        
            #Verify if the substring of the sequence from i to i+3 is 3 letters long.
            if len(self.sequence[i:i+3])==3:
           
            #If three letters are correct, then it will be appended to self.codons.
                self.codons.append(self.sequence[i:i+3])

#Translates the RNA sequence into a protein sequence by looking up each codon in the 'standard_code' dictionary.
#If the codon does not exist, it adds 'X' to denote unknown amino acid
#Essentially, 1.initialize, 2.iterate, 3.verify valid codon. 4.Return the final translated protein sequence 
    def translate(self):
    	
    	#Initialize an empty string to hold the translated protein sequence
        self.translated_sequence = ''

	#To iterate over each codon in the list self.codons
        for codon in self.codons:
        
	    #Verify if the codon is not present in the 'standard_code', if not found in dictionary, then true
            if codon not in list(standard_code.keys()):
                
                #If unknown, then appends 'X' to the 'self.translated_sequence'
                self.translated_sequence += 'X' 
             
            #If found, then proceed; appends the corresponding amino acid to 'self.translated_sequence'   
            else:
                self.translated_sequence += standard_code[codon]
        
        #After processing all the codons, the method returns the final       
        return self.translated_sequence

#Step-5: Create Protenin class (must inherit Seq class)
class Protein(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species)
        self.__str__()

#replaces any characters in the protein sequence that are not amino acids based on 'standard_code' dictionary 
#with N to denote unknown amino acid
    def __str__(self):
        LETTER = ''.join(list(standard_code.values()))
        self.sequence = re.sub('[^'+LETTER+']','N', self.sequence) 
        return self.sequence

#Calculates the total hydrophobicity of the protin sequence by summing the hydrophobicity values from
#'kyte_dolittle' dictionary then values are summed up to give an overall hydrophobicity score for the protein
#1.initalize. 2.Iterate through each amino acid 3.look up its hydrophobicity value 4.return total hydrophobicity value 
    def total_hydro(self):
    	
    	#initializes an attribute and set it to 0.
	#will accumulate the hydrophobicity values of each amino acid in the protein sequence as it iterates
        self.hydrophobicity = 0
        
        #looks up the hydrophobicity value for that amino acid in the kyte_doolittle dictionary
        #adds it to the total self.hydrophobicity
        for let in self.sequence:
            self.hydrophobicity += kyte_doolittle[let]
        
        #returns the final total hydrophobicity score
        #represents the overall hydrophobicity of the entire protein sequence.
        return self.hydrophobicity

#Calculates the molecular weight of the protein sequence by summiung the molecular weights of each amino acid from
#'aa_mol_weights' dictionary
    def mol_weight(self):
    
    	#initializes an attribute and set it to 0.
        self.MW = 0
        
        #Looks up molecular weight for that amino acid in the "aa_mol_weights" dictionary then add to "self.MW"
        for let in self.sequence:
            self.MW += aa_mol_weights[let]
        
        #Return molecular weight
        return self.MW
        
        
        
####EXAMPLE OUTPUT Starts below


# Create a DNA object
x = DNA("ATGCGTAA", "Gene_X", "Species_Y", "G001")

# Print the DNA object 
#calls __str__ method of DNA class, which inherits from the Seq class
#__init__ method of Seq is also called to initialize those attributes
print(x)

# Perform GC content analysis
#calculates the GC content by counting the occurrences of 'G' and 'C' in the DNA sequence
#using regular expressions (re.findall()).
print("GC Content:", x.analysis())

# Get the reverse complement of the DNA sequence
#calculates the reverse complement of the DNA sequence by first reversing the sequence
#Then, replaces  each base with its complement
print("Reverse Complement:", x.reverse_complement())

# Get six-frame translation of the DNA sequence
#Returns six different reading frames—three from the forward sequence and three from the reverse complement
print("Six Frame Translations:", x.six_frames())

# Print detailed DNA information
x.print_info()

# Generate a FASTA format string for the DNA sequence
print(x.fasta())

####EXAMPLE OUTPUT Ends



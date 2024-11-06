"""
Eddie Cho
BIOL 668
SP22
Pythn OOP Project
This program analyzes given seqeunce(s) and yields outputs
in terms of DNA, RNA,Protein using Objected Orientated Programming(OOP) concepts
"""


#Step-1: Initialize values
import re #Imports the regular expression (re) module for pattern matching

#Dictionary to map RNA codons to amino acids or stop codons
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

#Dictionary with Kyte-Doolittle hydropathy index for amino acids
kyte_doolittle={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3}
#Dictionary with molecular weights of amino acids
aa_mol_weights={'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
                'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
                'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
                'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19}

#Step-2: Create Seq class
class Seq:
    def __init__(self,sequence,gene,species): #Constructor to initialize the sequence, gene, and species
        self.sequence=sequence #Assigns the sequence to the instance
        self.__str__()  #Calls the __str__ method to process the sequence
        self.gene=gene #Assigns the gene name to the instance
        self.species=species #Assigns the species name to the instance
        self.kmers = [] #Initializes an empty list to store kmers

    def __str__(self): #Method to clean and return the sequence as a string
        self.sequence = self.sequence.strip() #Strips any leading or trailing whitespace
        self.sequence = self.sequence.upper() #Converts the sequence to uppercase
        return self.sequence #Returns the cleaned sequence

    def print_record(self): #Method to print the species and gene with sequence
        print(self.species + " " + self.gene + ": " + self.sequence)

    def make_kmers(self, k=3): #Method to create k-mers (subsequences of length k) from the sequence
        for i in range (len(self.sequence) - (k-1)): #Iterates through the sequence to generate kmers
            if self.sequence[i:i+k] not in self.kmers: #If the k-mer isn't already in the list,
              self.kmers.append(self.sequence[i:i+k]) #Adds the k-mer to the list

    def fasta(self): #Method to return the sequence in FASTA format
        self.fasta_str = '> ' + self.species + ' ' + self.gene + '\n'+self.sequence
        return self.fasta_str #Returns the formatted FASTA string

#Step-3: Create DNA class (must inherit Seq class)
class DNA(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):  #Constructor to initialize DNA attributes
        super().__init__(sequence,gene,species) #Calls the constructor of the parent Seq class
        self.sequence=sequence #Assigns the sequence to the instance
        self.__str__() #Processes the sequence (uppercase and strip) through __str__    
        self.geneid=geneid #Assigns the gene ID to the instance
        self.species = species #Assigns the species name to the instance
        self.gene = gene #Assigns the gene name to the instance

    def __str__(self): #Method to clean the sequence (remove non-ATGC characters and replace with 'N')
        super().__str__() #Calls the parent class's __str__ method
        self.sequence = re.sub('[^ATGCU]','N', self.sequence) #Replaces non-ATGC characters with 'N'
        return self.sequence #Returns the cleaned sequence

    def analysis(self): #Method to calculate the GC content of the DNA sequence
        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence)) #Counts the G and C bases
        return gc #Returns the GC content

    def reverse_complement(self): #Method to generate the reverse complement of the DNA sequence
        complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', "N":"N"}  #Dictionary for base pairing
        self.reverse_complement_sequence = '' #Initialize the reverse complement sequence
        for j in reversed(self.sequence): #Iterate over the sequence in reverse
            self.reverse_complement_sequence+= complement_dict[j] #Add the complement of each base
        return self.reverse_complement_sequence  #Returns the reverse complement sequence

    def print_info(self): #Method to print gene ID, species, gene, and sequence
        print(self.geneid + " " +self.species+' '+self.gene+' ' +self.sequence)

    def six_frames(self): #Method to generate six reading frames of the DNA sequence
        self.SixFrames = [] #Initialize an empty list for the six frames
        
        for i in range(3): #Three forward reading frames
            self.SixFrames.append(self.sequence[i:]) #Adds the sequence starting at positions 0, 1, and 2
            self.reverse_complement_sequence = self.reverse_complement() #Get the reverse complement sequence
        for i in range(3): #Three reverse reading frames
            self.SixFrames.append(self.reverse_complement_sequence[i:]) #Adds the reverse complement sequence
        return self.SixFrames  #Returns the six frames

#Step-4: Create RNA class(Must inherit DNA class)
class RNA(DNA):

    def __init__(self,sequence,gene,species,geneid,**kwargs): #Constructor to initialize RNA attributes
        super().__init__(sequence,gene,species,geneid)  #Calls the constructor of the parent DNA class
        self.__str__()  #Processes the sequence through __str__ method
        self.codons = [] # Initializes an empty list to store codons

    def __str__(self): #Method to clean the sequence (replace 'T' with 'U' for RNA)
        super().__str__()  #Calls the parent class's __str__ method
        self.sequence = re.sub('T', 'U', self.sequence) #Replaces all 'T's with 'U's
        return self.sequence #Returns the RNA sequence

    def make_codons(self): #Method to generate codons from the RNA sequence
        for i in range (0,len(self.sequence) - 2, 3): #Iterates through the sequence in steps of 3
            if len(self.sequence[i:i+3])==3: #Ensures the codon is of length 3
                self.codons.append(self.sequence[i:i+3]) #Adds the codon to the list

def translate(self):  #Method to translate the RNA sequence to a protein sequence
        self.translated_sequence = ''  #Initialize an empty string for the translated sequence

        for codon in self.codons:  #Iterate through each codon
            if codon not in list(standard_code.keys()):  # If the codon is not valid, add 'X'
                self.translated_sequence += 'X'
            else:  #Otherwise, add the corresponding amino acid
                self.translated_sequence += standard_code[codon]
        return self.translated_sequence  #Returns the translated protein sequence

#Step-5: Create Protenin class (must inherit Seq class)
class Protein(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs): #Constructor to initialize Protein attributes
        super().__init__(sequence,gene,species) #Calls the constructor of the parent Seq class
        self.__str__() #Processes the sequence through the __str__ method

    def __str__(self): #Method to clean the protein sequence (remove invalid characters)
        super().__str__() #Calls the parent class's __str__ method
        LETTER = ''.join(list(standard_code.values())) #concatenates all valid amino acid letters
        self.sequence = re.sub('[^'+LETTER+']','N', self.sequence)   #Replace invalid amino acids with 'N'
        return self.sequence #Returns the cleaned sequence

    def total_hydro(self): #Method to calculate the hydrophobicity of the protein sequence
        self.hydrophobicity = 0 #Initialize hydrophobicity score
        for let in self.sequence: #Iterates over each amino acid in the sequence
            self.hydrophobicity += kyte_doolittle[let] #Adds the hydropathy value for each amino acid
        return self.hydrophobicity #Returns the total hydrophobicity

    def mol_weight(self): #Method to calculate the molecular weight of the protein sequence
        self.MW = 0  #Initialize molecular weight
        for let in self.sequence: #Iterates over each amino acid in the sequence
            self.MW += aa_mol_weights[let] #Adds the molecular weight for each amino acid
        return self.MW #Returns the total molecular weight

x=DNA("G","tmp","m",000) #Creates a new DNA object with a sequence "G", gene "tmp", species "m", and geneid 000 for output

print(x)

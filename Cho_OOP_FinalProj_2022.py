"""
Eddie Cho
BIOL 668
SP22
Pythn OOP Project
This program analyzes given seqeunce(s) and yields outputs
in terms of DNA, RNA,Protein using Objected Orientated Programming(OOP) concepts
"""


#Step-1: Initialize values
import re

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

kyte_doolittle={'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
                'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3}

aa_mol_weights={'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
                'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
                'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
                'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19}

#Step-2: Create Seq class
class Seq:
    def __init__(self,sequence,gene,species):
        self.sequence=sequence
        self.__str__() 
        self.gene=gene
        self.species=species
        self.kmers = []

    def __str__(self):
        self.sequence = self.sequence.strip()
        self.sequence = self.sequence.upper()
        return self.sequence

    def print_record(self):
        print(self.species + " " + self.gene + ": " + self.sequence)

    def make_kmers(self, k=3):
        for i in range (len(self.sequence) - (k-1)):
            if self.sequence[i:i+k] not in self.kmers:
              self.kmers.append(self.sequence[i:i+k])

    def fasta(self):
        self.fasta_str = '> ' + self.species + ' ' + self.gene + '\n'+self.sequence
        return self.fasta_str

#Step-3: Create DNA class (must inherit Seq class)
class DNA(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species)
        self.sequence=sequence
        self.__str__()    
        self.geneid=geneid
        self.species = species
        self.gene = gene

    def __str__(self):
        super().__str__()
        self.sequence = re.sub('[^ATGCU]','N', self.sequence) 
        return self.sequence

    def analysis(self):
        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence))
        return gc

    def reverse_complement(self):
        complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', "N":"N"}
        self.reverse_complement_sequence = ''
        for j in reversed(self.sequence):
            self.reverse_complement_sequence+= complement_dict[j]
        return self.reverse_complement_sequence

    def print_info(self):
        print(self.geneid + " " +self.species+' '+self.gene+' ' +self.sequence)

    def six_frames(self):
        self.SixFrames = []
        
        for i in range(3):
            self.SixFrames.append(self.sequence[i:])
            self.reverse_complement_sequence = self.reverse_complement()
        for i in range(3):
            self.SixFrames.append(self.reverse_complement_sequence[i:])
        return self.SixFrames

#Step-4: Create RNA class(Must inherit DNA class)
class RNA(DNA):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species,geneid)
        self.__str__()
        self.codons = []

    def __str__(self):
        super().__str__()
        self.sequence = re.sub('T', 'U', self.sequence)
        return self.sequence

    def make_codons(self):
        for i in range (0,len(self.sequence) - 2, 3):
            if len(self.sequence[i:i+3])==3:
                self.codons.append(self.sequence[i:i+3])
    def translate(self):
        self.translated_sequence = ''

        for codon in self.codons:
            if codon not in list(standard_code.keys()):
                self.translated_sequence += 'X'
            else:
                self.translated_sequence += standard_code[codon]
        return self.translated_sequence

#Step-5: Create Protenin class (must inherit Seq class)
class Protein(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species)
        self.__str__()

    def __str__(self):
        super().__str__()
        LETTER = ''.join(list(standard_code.values()))
        self.sequence = re.sub('[^'+LETTER+']','N', self.sequence) 
        return self.sequence

    def total_hydro(self):
        self.hydrophobicity = 0
        for let in self.sequence:
            self.hydrophobicity += kyte_doolittle[let]
        return self.hydrophobicity

    def mol_weight(self):
        self.MW = 0
        for let in self.sequence:
            self.MW += aa_mol_weights[let]
        return self.MW

x=DNA("G","tmp","m",000)


#!/usr/bin/env python
# coding: utf-8

# In[ ]:


class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        """ 
        The init method initializes an instance of the class with the attribute "protein".
        The string is converted to uppercase and a new sequence with only valid amino acids joined is stored.
        """
        
        up_protein = str.upper(protein) #Converts string to uppercase
        split_aminos = [] #stores split characters in empty list
        valid_aminos = self.aa2mw.keys() #Obtains keys from the dictionary that has only valid amino acids
        
        #Iterates over each character in the uppercase "protein" string
        for char in up_protein:
            #If character is in valid amino acids, it is appended to the split_aminos list
            if char in valid_aminos:
                split_aminos.append(char)
        
        p = ''.join(split_aminos).split() #joins spli_aminos list and splits into one string
        self.protein = ''.join(p).upper() #sets final uppercased string as value for self.protein
        #Overrides previous value with input protein
        self.protein = protein 
        
        #Creates empty dictionary with 20 amino acids as keys and initial count as 0
        self.aaComp = {
            'A': 0, 'G': 0, 'M': 0, 'S': 0, 'C': 0,
            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
            'W': 0, 'F': 0, 'L': 0, 'R': 0, 'Y': 0
            }
        
        #Iterates over each key in dictionary of valid aminos
        for aa in self.aa2mw.keys():
            self.aaComp[aa] = self.protein.count(aa) #Each valid amino is counted and stored in the aaComp dictionary
        

    def aaCount (self):
        """ 
        aaCount returns the number of total amino acids in the protein sequence.
        This is after uppercased, split, and valid aminos joined.
        """
        
        #Iterates over each character in protein sequence
        for aa in self.protein:
            #If character is a valid amino, according to keys in self.aaComp
            if aa in self.aaComp.keys():
                return sum(self.aaComp.values()) #returns sum of all values in the dictionary
        

    def pI (self):
        """ 
        pI calculates the isolectric point of the protein by referring to another function _charge_.
        The only input is just the initial protein seqeunce, and outputs theoretical pI.
        """
        if sum(self.aaComp.values()) == 0:
            return 0.0
        else:
            
            #adapted from Srikar's code, information found online
            temp_pH = 0.0 #initialize
            best_pH = 42.0 #set to an arbitrarily high value so that first pH value in iteration is guranateed to be less
            best_charge = 999 #same reasoning as above, except with charge
            
            #iterate over pH values until the pH value is less than 14
            while temp_pH < 14:
                charge = abs(self._charge_(temp_pH)) #absolute value of the charge at a given pH
                #If current charge is less than or equal to best charge seen so far, store the current pH as best pH and current charge as best charge
                if charge <= best_charge: 
                    best_charge = charge
                    best_pH = temp_pH
                temp_pH += 0.01 #increment by 0.01 after each iteration
            return round (best_pH, 2)
    

    def aaComposition (self) :
        """ 
        aaComposition returns the composition percentages of amino acids given in the protein.
        Composition is stored as a dictionary where each key represents an amino acid, and corresponding value 
        is the number of occurences of that amino in protein. 
        """
        return self.aaComp
        

    def _charge_ (self, pH):
        """ 
        aaComposition returns the composition percentages of amino acids given in the protein and uses pH.
        Composition is stored as a dictionary where each key represents an amino acid, and corresponding value 
        is the number of occurences of that amino in protein. 
        """
        
        pH = round(pH, 2) #initialize pH to 2 decimal places
        #initialize to store each
        pos = 0.0
        neg = 0.0
        num = 0.0
        den = 0.0
        
        #Iterate tho each amino acid in the protein
        for aa in self.protein:   
            #If the amino acid is present in the positive charge dictionary
            if aa in self.aa2chargePos.keys():   
                num = 10 ** self.aa2chargePos[aa]  #calculate numerator for positive charge calculation 
                den = 10 ** self.aa2chargePos[aa] + (10 ** pH)  #calculate denominator for positive charge calculation 
                pos += num / den #update
            #If the amino acid is present in the negative charge dictionary    
            elif aa in self.aa2chargeNeg.keys():
                num = 10** pH #calculate numerator for negative charge calculation 
                den = 10 ** self.aa2chargeNeg[aa] + (10 ** pH) #calculate numerator for negative charge calculation 
                neg += num / den #update
            
        pos += (10**self.aaNterm)/(10**self.aaNterm + 10**pH) # add N-terminal contribution to positive charge value
        neg += (10**pH)/(10**self.aaCterm+10**pH) # add C-terminal contribution to negative charge value
        return pos - neg # return the difference of positive and negative charge value
        
        


    def molarExtinction (self):
        tyro = self.aaComp['Y'] * self.aa2abs280['Y']
        trypto = self.aaComp['W'] * self.aa2abs280['W']
        cyst = self.aaComp['C'] * self.aa2abs280['C']
        
        mol_extinct = tyro + trypto + cyst
        
        return mol_extinct
        

    def massExtinction (self, Cysteine = True):
        """
        This function returns the mass extinction coefficient of the protein based on its molecular weight and molar extinction coefficient. 
        Cysteine is Boolean value indicating whether to include cysteine in the calculation. Default is True.
        """
        myMW =  self.molecularWeight() #calculates molecular weight by reference molecularWeight function
        #If not equal to zero, it returns the mass extinction coefficient
        if myMW != 0:
            return self.molarExtinction() / myMW #formula is to divide these two values, already given
        #If molecular weight is zero, so it molar extinction
        else: 
            return 0.0
        
    
    def molecularWeight (self):
        """
        This function calculates the molecular weight of the protein. 
        It adds up the molecular weight of each amino acid in the protein and subtracts the mass of water lost during peptide bond formation.
        Returns the molecular weight as a float.
        """
        
        # checks if the amino acid composition of the protein is available
        if sum(self.aaComp.values()) == 0:
            return 0.0 #if not yes to above check, return zero
        
        else:
            total = 0 #initialize a variable to store total molecular weight
            water_loss = (len(self.protein) - 1) # calculate the number of water molecules lost during peptide bond formation (in hint)
            # Iterate through each amino acid in the protein
            for aa in self.protein:
                total += self.aa2mw.get(aa) # add the molecular weight of the current amino acid to the total
            return total - (water_loss * self.mwH2O) # subtract the mass of water lost during peptide bond formation

# Please do not modify any of the following.  This will produce a standard output that can be parsed
import sys 
def main():
    """
     Prompts user to input a protein sequence, calculates protein parameters, and outputs results.
     Creates an instance of ProteinParam class, and calls the various methods to calculate the parameters. The results are printed.
    """
    inString = input('protein sequence?')
    while inString : 
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        
        if myAAnumber == 0 : myAAnumber = 1 #in the case where no amino acids are present
       
        for aa,n in sorted(myParamMaker.aaComposition().items(), 
                           key= lambda item:item[0]):
            print ("\t{} = {:.2%}".format(aa, n/myAAnumber))
    
        inString = input('protein sequence?') #keeps asking for protein sequence input
        
    

if __name__ == "__main__":
    main()
    
#VLSPADKTNVKAAW


# In[ ]:


class NucParams:
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        """
        __init__ initializes the object by creating dictionaries to store the nucleotide, amino acid, codon composotion of each sequence.
        A string containing the seqeunce to be analyzed is the parameter.
        """
        self.codonComp = {codon: 0 for codon in self.rnaCodonTable} #stores frequency of each codon in the iput sequence
        self.aaComp = {aa: 0 for aa in self.rnaCodonTable.values()} #stores frequency of each amino acid in the input sequence
        self.nucComp = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0} #stores the frequency of each nucleotide
        
        #Add given sequence to the object
        self.addSequence(inString) 
    
    def addSequence(self, inSeq):
        """
        addSequence adds a sequence to the existing dictionaries and updates values accordingly.
        A string of nucleotides to be added to the existing sequence is the paramter. 
        """
        
        #convert sequence to uppercase and remove spaces
        inSeq = inSeq.upper().replace(' ','')
        
        #Iterate over sequence in codons, groups of 3
        for i in range(0, len(inSeq), 3):
            codon = inSeq[i : i +3]
           
        #Iterate over nucleotides in codon and update
            for nuc in codon:
                self.nucComp[nuc] += 1
            
            codon = codon.replace('T', 'U') #to get RNA codons
            
            #update amino acid and codon dictionaries
            if codon in self.rnaCodonTable:
                aa = self.rnaCodonTable[codon]
                self.aaComp[aa] += 1
                self.codonComp[codon] += 1
    #returns amino acid composition dictionary   
    def aaComposition(self):
        return self.aaComp
    
    #returns the nucleotide composition dictionary
    def nucComposition(self):
        return self.nucComp
    
    #retruns codon composition dictionary
    def codonComposition(self):
        return self.codonComp
    
    #returns total count of nucleotides in seqeunce, sum of calues in nuc comp
    def nucCount(self):
        return sum(self.nucComp.values())


# In[ ]:


import sys
class FastAreader :
    ''' 
    Define objects to read FastA files.
    
    instantiation: 
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=None):
        '''contructor: saves attribute fname '''
        self.fname = fname
            
    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is None:
            return sys.stdin
        else:
            return open(self.fname)
        
    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        
        with self.doOpen() as fileH:
            
            header = ''
            sequence = ''
            
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence


# In[ ]:


class OrfFinder():
    """
    Takes in a sequence of a FASTA file and stores ORFs in a list within a list.
    Attributes:
         List of valid stop codons.
         List of valid start codon.
         Dictionary of valid nucleotides.
    """
   
    nucPair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} #for later

    def __init__(self, seq):
        """Creates a list of open reading frames (ORFs) found in the input DNA sequence in the FASTA format. Each ORF is represented as a list containing the frame, 
        start position, stop position, and length 
    Parameters:
    sequence: A string representing the DNA sequence in the FASTA format.
        """
        self.seq = seq
        self.Orfs = []  # The lists where the ORFs frame, start, stop, and length are stored
        
        
    #start with forward ORFs and the do reverse, essentially split class into 2
    def findOrfs(self):
        """
        Find Orfs on 3'-5' strand and return list of Orfs.
        
        Returns:
        list: A list of dictionaries, where each dictionary represents an ORF and has the keys:
              'frame': the reading frame of the ORF (-1, -2, or -3)
              'start': the 1-based starting position of the ORF
              'stop': the 1-based ending position of the ORF
              'length': the length of the ORF in nucleotides
        
        Iterates over all possible frames and scans for start and stop codons to identify each ORF.
        Saves ORF found in each list. 
        
        Each iteration is seprated based on which list found on, used Boolean variables to keep track of certain conditions.
        After each condition is True/False, determines whether ORF should be saved or not. 
        """
        #startPos = []
        #foundStart = False
        #foundCodon = False

        for frame in range(3):  # Check first 3 frames
            foundStart = False  # After finding codons and start codons, set flag for start codon
            foundCodon = False #Set flag for any found codon
            startPos = []  # Starts/clears list for start position for each frame
            
            #-1 * ((frame % 3) + x) formar: 
            #the frame in which the ORF was found, where frame values -1, -2, and -3 denote ORFs on the reverse strand and frame values 1, 2, and 3 denote ORFs on the forward strand.
            
            for i in range(frame, len(self.seq), 3):
                codon = self.seq[i : i + 3] #get a codon, length of 3 nucleotides
                
                #check if codon is start codon, and if so continue
                if codon == 'ATG':  
                    startPos.append(i) #add current index to start position list
                    
                    #set flags to true because both found
                    foundStart = True 
                    foundCodon = True
                    
                #check if current codon is a stop codon, conintue if so
                if codon in ['TGA', 'TAG', 'TAA']:
                    
                    #also if start codon was previously found
                    if foundStart:
                        
                        #calculate start, stop, and length of the ORF found
                        start = startPos[0] + 1 - frame
                        stop = i + 3
                        length = stop - start + 1
                        self.saveOrf((frame % 3) + 1, start, stop, length)
                        
                        #save the ORF and set conditions accordingly
                        startPos = []
                        foundStart = False
                        foundCodon = True
                        
                        
                #if neither start nor stop codons were found
                if not foundCodon:
                    
                    #check if current codon is stop codon again
                    if codon in ['TGA', 'TAG', 'TAA']:  # If no start codon was found but stop codon found
                        
                        #create an ORF with a dangling stop codon
                        start = 1
                        stop = i + 3
                        length = stop - start + 1
                        
                        #save the ORF
                        self.saveOrf((frame % 3) + 1, start, stop, length)
                        startPos = [] #clear start position list and set flag True
                        foundCodon = True

            if foundStart:  # If no stop codon was found but start codon was found. - dangling start
                
                #create an ORF with a dangling start codon
                start = startPos[0] + 1
                stop = len(self.seq)
                length = stop - start + 1
                self.saveOrf((frame % 3) + 1, start, stop, length)

        return self.Orfs #return list of ORFs found in the sequence
    
    #went through test file after this was written to check for errors and functionality.
    #worked so will be using same format for the reverse ORF finder

    def reverseSequence(self):
        """ Returns the reverse complement of the DNA sequence.
    
            Uses a dictionary to match each nucleotide to its complement.
            Reverses the DNA sequence and replaces each nucleotide with complement.
    
            Returns: The reverse complement of the DNA sequence.
        """
        #nucPair = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} #dictionary to match each respective nucleotide (top)
        
        # Reverse the DNA sequence and get the complementary nucleotides using the dictionary
        # Returns the reversed and complementary strand of DNA
        return ''.join([self.nucPair[base] for base in self.seq[::-1]]) #uses dictionary value to put together
    
    def findRevOrfs(self): 
        """ Find Orfs on 5'-3' strand and returns that list of Orfs.
        Searches for start/ stop in all three frames of the reverse complement.
        
        When start codon is found,add to list.
        When stop codon is found, checked against each start position to see if in frame with the start codon before save.
        
        Returns:
        A list of all ORFs found on the reverse complement.
        """
        reverseSeq = self.reverseSequence()
        #startPos = []
        #foundStart = False
        #foundCodon = False

        for frame in range(3):  # Check  all frames, first 3
            foundStart = False  # Flag when finding codons and start codons.
            foundCodon = False
            startPos = []  # Clears the list to save position for each frame
            
            #iterates over indices of reverse sequence starting from frame to end
            for i in range(frame, len(reverseSeq), 3):
                codon = reverseSeq[i : i + 3]  # The codon is 3 nucleotides.

                if codon == 'ATG':  # When start codon is found, same as previous prgram
                    startPos.append(i)
                    foundStart = True
                    foundCodon = True

                if codon in ['TGA', 'TAG', 'TAA']:
                    if foundStart:  #if both start and stop codons are present
                        stop = len(reverseSeq) - startPos[0] #distance from start of reverse to start codon
                        start = len(reverseSeq) - (i + 2) #distance from end of reverse to position of last nuc, subtract from total of reverse
                        
                        #adjust stop position depending on the frame
                        if frame == 1: 
                            stop += 1
                        
                        elif frame == 2: 
                            stop += 2
                        
                        length = stop - start + 1 #calculate length after adjusting
                        
                        #save ORF and reset variables
                        #computes the reading frame of the ORF relative to the 3'-5' strand of the DNA sequence, but with a negative sign because reverse complement
                        self.saveOrf(-1 * ((frame%3) + 1), start, stop, length)
                        
                        startPos = []
                        foundStart = False
                        foundCodon = True
                #dangling stop = stop codon was found without a corresponding start codon
                if not foundCodon:
                    if codon in ['TGA', 'TAG', 'TAA']:  # If no start codon was found but stop codon found, make dangling stop
                        start = len(reverseSeq) - i - 2 #calculate start position
                        stop = len(reverseSeq)
                        length = stop - start + 1
                        
                        self.saveOrf(-1 * ((frame % 3) + 2), start, stop, length)
                        startPos = []
                        foundCodon = True

            if foundStart:  # If no stop codon was found but start codon was found, make dangling start
                start =  startPos[0] + 1
                stop = 1
                length = stop - start + 1 
                self.saveOrf(-1 * ((frame % 3) + 1), start, stop, length)

        return self.Orfs

    def saveOrf(self, frame, start, stop, length):
        """ Saves ORF info in following format:
        
        Adds list to the Orfs attribute, each list contains:
        - Frame: the reading frame of the ORF (0, 1, or 2)
        - Start: the index of the start codon for the ORF
        - Stop: the index of the stop codon for the ORF
        - Length: the length of the ORF in nucleotides
        """
        self.Orfs.append([frame, start, stop, length]) #puts all together


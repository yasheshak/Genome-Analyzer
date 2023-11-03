#!/usr/bin/env python
# coding: utf-8

# In[1]:


from SequenceAnalysis import FastAreader
from SequenceAnalysis import NucParams
#adds seqeunce to NucParams object
def main (fileName = None):
    myReader = FastAreader('testGenome.fa') 
    myNuc=NucParams()
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)
    
    #calculate AT and GC counts and content
    A_count = myNuc.nucComp.get('A')
    T_count = myNuc.nucComp.get('T')
    AT_count = A_count + T_count
    
    G_count = myNuc.nucComp.get('G')
    C_count = myNuc.nucComp.get('C')
    GC_count = G_count + C_count
    GC_content = (GC_count / myNuc.nucCount())*100
    
    #print sequence length and GC content
    print("Sequence Length = {0:0.02f}Mb\n".format(myNuc.nucCount()/1000000))
    print("GC Content = {0:0.01f}%\n".format(GC_content))
    
    # sort codons in alpha order, by Amino Acid
    aaDict = myNuc.aaComposition()
    RNAList = []
    for codon, aa in myNuc.rnaCodonTable.items():
        RNAList.append((aa, codon))
    RNAList = sorted(RNAList)
    
    # calculate relative codon usage for each codon and print
    for aa, codon in RNAList:
        val = myNuc.codonComp[codon] / myNuc.aaComp[aa]
        
        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(codon, aa, val*100, myNuc.codonComp[codon]))

if __name__ == "__main__":
    main()
    


# In[ ]:





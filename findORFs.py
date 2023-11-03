#!/usr/bin/env python3
# Name: Yashesha Kothari (ykothari)
# Group Members: Srikar Bevara, Donya Mirzazadeh, Cameron Ahyahee, Tia Abrahama


########################################################################
# CommandLine
########################################################################

from SequenceAnalysis import OrfFinder, FastAreader
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond. 
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
 
    '''
    
    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''
        
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - a brief description of what this thing does', 
                                             epilog = 'Program epilog - some other stuff you feel compelled to say', 
                                             add_help = True, #default is True 
                                             prefix_chars = '-', 
                                             usage = '%(prog)s [options] -option1[default] <input >output'
                                             )
        self.parser.add_argument('-lG', '--longestGene', action = 'store', nargs='?', const=True, default=False, help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices= (100,200,300,500,1000), default=100, action = 'store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action = 'append', default = ['ATG'],nargs='?', 
                                 help='start Codon') #allows multiple list options
        self.parser.add_argument('-t', '--stop', action = 'append', default = ['TAG','TGA','TAA'],nargs='?', help='stop Codon') #allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')  
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

# Main
# Here is the main program
   
def main(inFile = None):
    '''
    Find some genes. 
    
    Can take an input file in FASTA format containing one or more DNA sequences or read the input sequences 
    from the command line using the CommandLine class.
    '''
    
    #Check if input file is provides
    if inFile is None:
        thisCommandLine = CommandLine() #if not use CommanLine class for reading input sequences
        
        #Check is longest Gene argument passed
        if thisCommandLine.args.longestGene:
            fastaFile = FastAreader() #create an instance of the FastAreader class to read from FASTA
            
            #Iterate over each sequence and print header of current
            for header, sequence in fastaFile.readFasta():
                print(header)
                
                #Create instance of the OrfFinder class for current sequence and find ORFs
                OrfData = OrfFinder(sequence)
                OrfData.findOrfs()
                OrfData.findRevOrfs()
                
                #Filter ORFs based on length with minGene
                cleanedList = filter(lambda Orf: Orf[3] > thisCommandLine.args.minGene, OrfData.Orfs)
                
                #Iterate to sort ORFs by length and print each info piece in wanted format
                for frame, start, stop, length in sorted(cleanedList, key = lambda Orf: Orf[3], reverse = True):  # key specified function to be called on each element, returns length
                    print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(frame, start, stop, length)) #formars and prints in that order
                
                #if longest Gene not passed, use input file specified
                else:
                    thisCommandLine = CommandLine(inFile)
                print(thisCommandLine.args)#print arguments pass to the object
                
if __name__ == "__main__":
    main() #call main function, deleted given code because calling from command line
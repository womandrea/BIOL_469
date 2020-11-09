import argparse
import pandas as pd
import subprocess
import sys
import os
from Bio.Seq import Seq

'''
The purpose of this script is to run annotate on the 
required fasta file(s), ideally given an assembled phage genome. 

The script will store the PHANOTATE output, store all the saved positions, 
and extract the expected ORFs from the given fasta file for subsequent analysis. 

Expected analysis 
- GO term 
- BLASTing 
- other stuff???
'''
phanotate_file_path = "/home/andrea_gabriela_alexei_gmail_com/PHANOTATE/phanotate.py"

class ORF(object):
    """
    A class representing each potential open reading frame detected by Phanotate.
    """
    def __init__(self):
        self.frame = None # True if positive, False if negative
        self.start = 0 # beginning position
        self.end = 0 # end position (+1 to the val for substring stealing!
        self.contig = "" # for later idenfitication if assembly isn't resolved
        self.score = 0 # Learn more about the score, might use for initial filtering? Might be better to subproc it to awk


class Identification(object):
    """
    A class whose operations contribute to the identification of the 
    open reading frames within the fasta file input to Phanotate.
    """
    def __init__(self, args):
        self.input_file = args.input 
        if args.output_dir:
            self.out_file = args.output_dir + "/" + "phanotate_output.tsv"
        else: 
            self.out_file = os.getcwd() + "/" + "phanotate_output.tsv"
        self.orf_file = None
        self.get_output()
        self.write_ORF(self.retrieve_ORF())

    def write_ORF(self, dic):
        # Write a fasta file using the dictionary 
        #get the lists to make our fasta file!
        data = []
        counter = 1
        for i in dic:
            new_line = ">{} ORF{} Score: {} Frame: {} \n {}\n".format(i.contig, counter, i.score, i.frame, dic[i])
            data.append(new_line)
            counter += 1
        
        orf_file_output = "/".join(self.out_file.split("/")[:-1]) + "/"  + "phanotate_output_ORF.fasta"
        with open(orf_file_output, "w") as fasta:
            for i in dic:
                fasta.writelines(data)
        return None

    def get_output(self):
        """
        This code is meant to only save the output of Phanotate, and store 
        the file path as a variable of the class object

        get_output: self -> None 
        Effects: Updates self to include ORF file path from 
                 Phanotate. Writes the output of Phanotate to a file.
        """
        cmd = ["python", phanotate_file_path,  self.input_file]
        output_reading = subprocess.Popen(cmd,stdout=subprocess.PIPE)
        output_file_path = self.out_file
        output_file = open(output_file_path, "w")
        out, err = output_reading.communicate()
        out = out.decode(encoding="utf-8")
        output_file.write(out)
        output_file.close()

    def retrieve_ORF(self):
        """
        Calls for the storage of our ORFs as separate classes, 
        with parameters that will help us rewrite out FASTA. This function 
        only updates and returns the dictionary with the substring that 
        is the ORF.

        retrieve_ORF: -> Dict
        Effects: Updates dictionary from store_ORF.
        """
        orf_dict = self.store_ORF()
        with open(self.input_file, "r") as fasta_file:
            fasta_input = fasta_file.read().split(">")
            for orf in orf_dict:
                
                contig_name = orf.contig
                i = 0
                
                while contig_name not in fasta_input[i]:
                    i += 1

                contig_sequence = fasta_input[i].split("\n")[1:]
                contig_sequence = "".join(contig_sequence)
                 
                if orf.frame == "+":
                    orf_sequence = contig_sequence[orf.start-1:orf.end]
                    orf_dict[orf] += orf_sequence
                
                else:
                    orf_sequence = contig_sequence[orf.end-1:orf.start]
                    seq = Seq(orf_sequence)
                    reverse_seq = seq.reverse_complement()
                    orf_dict[orf] = reverse_seq
        return orf_dict

    
    def store_ORF(self):
        """
        Stores the information from the output into a dictionary. Will be called by retreive_ORFs to get 
        the necessary information. 

        store_ORF: self -> Dict
        """
        # STEP 6: Write all the information from the dictionary as a fasta output for later analysis.
        
        with open(self.out_file, "r") as fin:
            data = fin.read().splitlines(True)
            data = [data[1].split("#")[1]] + data[2:]
        with open(self.out_file, "w") as fout:
            fout.writelines(data)


        with open(self.out_file, "r"):
            tsv_output = pd.read_csv(self.out_file, sep="\t")
        tsv_output.reset_index(level=0, inplace=True)
        row_count = len(tsv_output.index)

        i = 0
        orf_dict = {}
        while i < row_count:
            """
                For each row, we are going to store the information of
                each hit in the ORF class. The key to this will map back to a
                string, which will be updated to be the substring 
                mapping to the ORF. 
            """
            new_orf = ORF()
            new_orf.start = int(tsv_output.iloc[i,0])
            new_orf.end = int(tsv_output.iloc[i,1])
            new_orf.frame = (tsv_output.iloc[i,2])
            new_orf.contig = tsv_output.iloc[i,3]
            new_orf.score = float(tsv_output.iloc[i,4])
            orf_dict[new_orf] = ""
            i += 1
        return orf_dict
        

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="For BIOL 469. Enjoy! Annotate your phage genomes for half of the efficiency of the original pipeline!")
    parser.add_argument("-in", help="The fasta files of the sequence data. May be a resolved assembly or separate contigs", required=True, dest="input")
    parser.add_argument("-out", help="Provide the full directory of where the output is to be stored. Otherwise saves file in the directory of the input.", dest=
"output_dir")
    args = parser.parse_args()
    Identification(args)

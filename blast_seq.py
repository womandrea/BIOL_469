import argparse
import pandas as pd
import subprocess
import sys
import os
from Bio.Seq import Seq
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from collections import Counter 
from concurrent import futures
import multiprocessing
import time

class prot_seq(object):
    def __init__(self):
        self.fasta_line = ""
        self.label = ""
        self.name = ""
        self.len = 0

class BLAST(object):
    def __init__(self, args):
        self.input_file =  args.input
        self.blast_call = args.alg
        
        if args.output_dir:
            if os.path.isdir(args.output_dir):
                self.output_dir = args.output_dir
            else:
                os.mkdir(args.output_dir)
                self.output_dir = args.output_dir
        else: 
            self.output_dir = os.getcwd()
        
        self.sequence = {}  # An empty dictionary that gets updated by the following call. 
        self.find_sequence()
        self.identify_sequences()
    
    def identify_sequences(self):
        print("So anyway, I just started BLASTing...")
        if self.blast_call == "blastx":
            file_path = self.output_dir + "/blastx_summary.tsv"
            self.blastx_call_all()
            row_input = [["Contig Name", "length_bp", "Blast Result"]]
               
            for value in self.sequence.values():
                new_row = [value.name, value.len, value.label]
                row_input.append(new_row)
                tsv = pd.DataFrame(row_input)
                tsv.to_csv(file_path, sep="\t")

             
        else: 
            file_path = self.output_dir + "/blastn_summary.tsv"
            blast_results_all = self.blastn_call_all()
            row_input = [["Contig Name", "length_bp", "Blast Result"]]
                
            for value in blast_results_all.values():
                new_row = [value.name, value.len, value.label]
                row_input.append(new_row)
                tsv = pd.DataFrame(row_input)
                tsv.to_csv(file_path, sep='\t')
        return None


    def find_sequence(self):
        with open(self.input_file, "r") as fin:
            all_lines = fin.read().split("\n")
            
            new_lines = []
        
            new_string = ""
            for sti in all_lines:
                if ">" in sti:
                    new_lines.append(new_string)
                    new_string = sti + "\n"
                else:
                    new_string += sti
            new_lines.append(new_string)
            new_lines = new_lines[1:]
            
            for line in new_lines:
                new_protein_sequence = prot_seq()
                new_protein_sequence.fasta_line= line
                new_protein_sequence.len = len(line.split("\n")[1])
                fasta_title = line.split("\n")[0]
                try: 
                    #new_protein_sequence.name = fasta_title[0] + "pos:{}-{}".format(fasta_title[1], fasta_title[2])
                    # DEFINE YOUR NAME HERE!
                    self.sequence[new_protein_sequence] = new_protein_sequence
                except: 
                    print("Oh oh oh! You didn't update how you wanted this fasta header parsed. Get on over to line 90 and fix me up!")
                
        return None
        
    def blastx_call_all(self):
        
        with futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor: 
            orf = tuple(self.sequence.values()
            for results in executor.map(lambda x: BLAST.blastx_call_each(self.sequence, x),orf):
                pass
        
        return None
    
    def blastn_call_all(self):
        dic = self.sequence
        
        with futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor: 
            orf = tuple(self.sequence.values())
            for results in executor.map(lambda x: BLAST.blastn_call_each(dic, x),orf):
                pass
        return None

    @staticmethod
    def blastx_call_each(dic, key):
        """
        Calls blast on each individual sequence. If a time-out (ValueError) is encountered
        the command is tried again, after a 5 second rest period. Command is re-attempted for a max of 10 times. 

        blast_call_each: Dic Str -> None
        Effects: Remote blating of sequencing, updates dictionary with the most frequently 
        occuring output. 

        """
        fasta_sequence = dic[key].fasta_line
        
        for attempt in range(10):
            try:
                top_hits = []
                result_handle = NCBIWWW.qblast("blastx", "nr", fasta_sequence)
                for record in NCBIXML.parse(result_handle):
                    if record.alignments:
                        for align in record.alignments:
                            for hsp in align.hsps:
                                if hsp.expect < 1e-20:
                                    top_hits.append(align.title[:100])
            
                top_genes = list(map(lambda x: x.split("|"), top_hits))
                top_genes = list(map (lambda x: "".join(x), top_genes))
                top_genes = list(map (lambda x: x.split(","), top_genes))
                top_genes = list(map (lambda x: "".join(x), top_genes))
            
                result = [item for items, c in Counter(top_genes).most_common() 
                                      for item in [items] * c]
                try:
                    dic[key].label = result[0]
                    print("Houston, we found a gene: {}".format(result[0]))
                except IndexError:
                    # If there are no genes identified 
                    dic[key].label = "NA"
                    print("Blast could not find a protein for one of these genes! Find me in the output tsv.")
            
            except ValueError:
                # ValueError occurs when too many sequences are given simulatenously to NCBI. This is meant to compensate for this!

                time.sleep(5)
                # print("Error! We're going to try this again! Strap in your seatbelt kids.")
            
            else: 
                #  When the ValueError does not occur, the for loop is broken
                break

        return None
    
    @staticmethod
    def blastn_call_each(dic, key):
        """
        Calls blast on each individual sequence. If a time-out (ValueError) is encountered
        the command is tried again, after a 5 second rest period. Command is re-attempted for a max of 10 times. 

        blast_call_each: Dic Str -> None
        Effects: Remote blating of sequencing, updates dictionary with the most frequently 
        occuring output. 

        """
        
        untranscribed = dic[key].fasta_line.split("\n")[1]
        transcribed = Seq(unstranscribed).transcribe()
        fasta_sequence =  dic[key].fasta_line.split("\n")[0] "\n" + transcribed
        for attempt in range(10):
            try:
                result_handle = NCBIWWW.qblast("blastn", "refseq_rna", fasta_sequence)
                for record in NCBIXML.parse(result_handle):
                   if record.alignments:
                        for align in record.alignments:
                            for hsp in align.hsps:
                                if hsp.expect < 1e-20:
                                    top_hits.append(align.title[:100])
            
                top_genes = list(map(lambda x: x.split("|"), top_hits))
                top_genes = list(map (lambda x: "".join(x), top_genes))
                top_genes = list(map (lambda x: x.split(","), top_genes))
                top_genes = list(map (lambda x: "".join(x), top_genes))
            
                result = [item for items, c in Counter(top_genes).most_common() 
                        for item in [items] * c]
                try:
                    dic[key].label = result[0]
                    print("Houston, we found a gene: {}".format(result[0]))
                except IndexError:
                    # If there are no genes identified 
                    dic[key].label = "NA"
                    print("Blast could not find a sequence for one of these guys! Oops!")

            
            except ValueError:
                # ValueError occurs when too many sequences are given simulatenously to NCBI. This is meant to compensate for this!
                time.sleep(5)
                # essentially repeats the process above if it gets a ValueError.   
            else: 
                 break
        return None 
      
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="WARNING: You have to update how to parse the file header for this to work! Make it based on your fasta file input.")
    parser.add_argument("-in", help="Input your nucleotide fasta file here!", required=True, dest="input")
    parser.add_argument("-out", help="Full file path for your blast summary!!", dest="output_dir")
    parser.add_argument("-alg [blastx, blastn]", help="Choose which one, this is required! Blastn will go for rna_sequences, blastx will go for proteins. Don't know the utility of this yet...",
    	dest="alg")
    args = parser.parse_args()
    BLAST(args)

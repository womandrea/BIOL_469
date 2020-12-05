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


class ORF(object):
    def __init__(self): 
        """
        Object representing important information from each respctive open reading 
        frame identified by Phanotate. 
        """

        self.frame = None # True if positive, False if negative
        self.start = 0 # beginning position
        self.end = 0 # end position (+1 to the val for substring stealing!
        self.contig = "" # for later idenfitication if assembly isn't resolved
        self.score = 0 # Learn more about the score, might use for initial filtering? Might be better to subproc it to awk
        
        self.prot = ""  # the protein sequence alone
        self.nuc = ""  # the nucleotide sequence alone
        self.prot_fas = ""  # the fasta file line corresponging to the protein sequence
        self.nuc_fas = ""  # the fasta file line corresponding to the nucleotide sequence
        self.label = ""  # the gene label for this open reading frame
  

class Identification(object):
    """
    Class encapsulating all the methods involing calling 
    Phanotate, parsing the output, and annotating the open reading frames. 
    """
    def __init__(self, args):
        self.input_file = args.input  # input fasta file which is the assembled phage genome
        
        if args.sample_name:  # for naming purposes alone
            self.sample_name = "{}_".format(args.sample_name)
        else:
            self.sample_name = "phanotate_"

        if args.output_dir:  # determining the output of all the processes
            self.output_dir = args.output_dir
            
            if not os.path.isdir(args.output_dir):
                os.mkdir(args.output_dir)
            else:
                None

            self.out_file = args.output_dir + "/" + "{}output.tsv".format(self.sample_name) # the output file for the initial phanotate calls

        else: 
            self.out_file = os.getcwd() + "/" + "{}output.tsv".format(self.sample_name)  # store the output in the current directory if none spec.
            self.output_dir = os.getcwd()
        
        self.orf_file = None  # the fna file
        self.all_orfs = None  # the ORF dictionary 
        self.prot_file = None  # the faa file 

        self.get_output() 
        self.write_ORF(self.retrieve_ORF())
        # Consider not making this 
        updated_orfs = Identification.blast_call_all(self.all_orfs)
        self.blast_output(updated_orfs)


    def write_ORF(self, dic):
        """
        Taking self and an ORF dictionary as an input, 
        writes the protein and nucleotide fasta files. 

        write_ORF: Identification Dic -> None
        Effects: Writes two new files, a protein fasta file and a nucleotide fasta 
                fasta files respectively. Updates protein_fas and nuc_fas parameters in the 
                dictionary. 
        """
        
        orf_data = []
        prot_data = []
        counter = 1
        
        for i in dic:

            # for each ORF in the dictionary, add the respective fasta file line to each list and update the ORF

            orf_new_line = ">{} ORF{} Score: {} Frame: {} \n {}\n".format(i.contig, counter, i.score, i.frame, dic[i].nuc)
            orf_data.append(orf_new_line)
            
            prot_new_line = ">{} ORF{} Score: {} Frame: {} \n {}\n".format(i.contig, counter, i.score, i.frame, dic[i].prot)
            prot_data.append(prot_new_line)


            dic[i].prot_fas = prot_new_line 
            dic[i].nuc_fas = orf_new_line

            counter += 1
        
        orf_file_output = "/".join(self.out_file.split("/")[:-1]) + "/"  + "{}ORF.fna".format(self.sample_name)
        
        with open(orf_file_output, "w") as fasta:  # write the fna file 
            fasta.writelines(orf_data)
        
        protein_file_output = "/".join(self.out_file.split("/")[:-1]) + "/"  + "{}ORF.faa".format(self.sample_name)
        
        with open(protein_file_output, "w") as fa:  # write the faa file 
            fa.writelines(prot_data)
            
        self.orf_file = orf_file_output
        self.prot_file = protein_file_output
        return None

    def get_output(self):
        """
        Saves the output of Phanotate, and stores 
        the file path as a variable of the class object

        get_output: self -> None 
        Effects: Updates self to include ORF file path from 
                 Phanotate. Writes the output of Phanotate to a file.
        """
        cmd = ["python","/home/andrea_gabriela_alexei_gmail_com/PHANOTATE/phanotate.py",  self.input_file]
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
                    orf_sequence = contig_sequence[orf.start-1:orf.end-1]
                    orf_dict[orf].nuc += Seq(orf_sequence)
                    rna_seq = orf_dict[orf].nuc.transcribe()
                    orf_dict[orf].prot = rna_seq.translate()

                
                else:
                    orf_sequence = contig_sequence[orf.end-1:orf.start-1]
                    seq = Seq(orf_sequence)
                    reverse_seq = seq.reverse_complement()
                    orf_dict[orf].nuc = reverse_seq
                    rna_seq = orf_dict[orf].nuc.transcribe()
                    orf_dict[orf].prot = rna_seq.translate()
        
        self.all_orfs = orf_dict
        return orf_dict

    
    def store_ORF(self):
        
        """
        Stores the information from the output into a dictionary. Will be called by retreive_ORFs to get 
        the necessary information. 

        store_ORF: self -> Dict
        Effects: Initialize ORF objects, store in an orf dictionary. 
        """
        
        with open(self.out_file, "r") as fin: 
            data = fin.read().splitlines(True)
            data = [data[1].split("#")[1]] + data[2:]
        with open(self.out_file, "w") as fout:
            fout.writelines(data)

        # create a dataframe from the Phanotate output to initialize each ORF object
        with open(self.out_file, "r"):
            tsv_output = pd.read_csv(self.out_file, sep="\t")
        tsv_output.reset_index(level=0, inplace=True)
        row_count = len(tsv_output.index)

        i = 0
        orf_dict = {}  # a little redundant not to make this automatically the parameter in self, this was an oversight from earlier!
        
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
            orf_dict[new_orf]  = new_orf
            i += 1
        return orf_dict
    
    def blast_output(self, dic):
        """
        Creates a summary tsv for the phage annotation. 

        blast_output: self Dic -> Dic 
        Effects: Creates and saves the tsv file storing blast identification 

        """

        new_output = self.output_dir + "/blast_summary.tsv"
        
        column_titles = ["ORF_Tag", "length_bp", "Product"]
        rows = [column_titles]

        for i in dic:
            orf = dic[i].nuc_fas
            orf_title = ("".join(orf.split(">")[1:])).split(" ")[1]  # save the ORF corresponding to each gene
            length_bp = len(dic[i].nuc)  # the length of each nucleotide sequence
            
            if dic[i].label:
                product_parsed = dic[i].label.split("[")[0] 
            else: 
                product_parsed = "NA"

            row = [orf_title, length_bp, product_parsed]
            rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(new_output, sep="\t")
        
        return None
    
    
    @staticmethod
    def blast_call_all(dic):
        """
        Multi-threads the blast call to avoid timing out errors. 

        blast_call_all: Dic -> Dic 

        TO UPDATE: Stop it from being a static method
        """
        with futures.ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor: 
            orf = tuple((dic.keys()))
            
            for results in executor.map(lambda x: Identification.blast_call_each(dic, x),
                    orf):
                pass
        
        return dic

    @staticmethod
    def blast_call_each(dic, key):
        """
        Calls blast on each individual sequence. If a time-out (ValueError) is encountered
        the command is tried again, after a 5 second rest period. Command is re-attempted for a max of 10 times. 

        blast_call_each: Dic Str -> None
        Effects: Remote blating of sequencing, updates dictionary with the most frequently 
        occuring output. 

        """
        fasta_sequence = dic[key].nuc_fas
        header = "".join(fasta_sequence.split("\n")[0])
        top_hits = []
        print("Starting on {}".format(header))
        
        for attempt in range(10):
            try:
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
                    dic[key].label = None
                    print("Blast could not find a gene for {}. Oops!".format(header))
            
            except ValueError:
                # ValueError occurs when too many sequences are given simulatenously to NCBI. This is meant to compensate for this!

                time.sleep(5)
                print("Error! We're going to try this again! Strap in your seatbelt kids.")
            
            else: 
                #  When the ValueError does not occur, the for loop is broken
                break

        return None 
      
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="For BIOL 469. Enjoy! Annotate your phage genomes for half of the efficiency of the original pipeline!")
    parser.add_argument("-in", help="The fasta files of the sequence data. May be a resolved assembly or separate contigs", required=True, dest="input")
    parser.add_argument("-sn", help="Sample name. Do not put any spaces, ex. E_coli_O157H7.", dest=
"sample_name")
    parser.add_argument("-out", help="Not required. Add the full filepath to the output directory. Stores files in current dir if unspecified.", dest="output_dir")
    args = parser.parse_args()
    Identification(args)

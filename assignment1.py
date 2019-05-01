import mysql.connector
import pysam
from pathlib import Path
import numpy as np
import subprocess

__author__ = 'Alexander Benkö'

class Assignment1:
    '''
    Class will calculate all relevant parameters on init.
    Print summary with print_summary() or call parameters directly.
    '''

    def __init__(self, gene, genome_reference, uscs_file_name, bamfile):
        self.gene = gene
        self.genome_reference = genome_reference
        self.uscs_file = uscs_file_name
        self.bamfile = bamfile
        self.alignfile = pysam.AlignmentFile(self.bamfile, "rb")

        ### FETCH DATA ###
        if not Path.cwd().joinpath(self.uscs_file).exists():
            print("Connecting to UCSC to fetch data")

            ## Open connection
            cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=self.genome_reference)

            ## Get cursor
            cursor = cnx.cursor()

            ## Build query fields
            query_fields = ["refGene.name2",
                            "refGene.name",
                            "refGene.chrom",
                            "refGene.txStart",
                            "refGene.txEnd",
                            "refGene.strand",
                            "refGene.exonCount",
                            "refGene.exonStarts",
                            "refGene.exonEnds"]

            ## Build query
            query = "SELECT DISTINCT {} from refGene WHERE name2 = '{}'".format(",".join(query_fields), self.gene)

            ## Execute query
            cursor.execute(query)

            ## Write to file
            ## TODO this may need some work
            with open(self.uscs_file, "w") as fh:
                for row in cursor:
                    fh.write(str(row) + "\n")


            ## Close cursor & connection
            cursor.close()
            cnx.close()

            print("Done fetching data")


        ### GENE COORDINATES ###
        # SELECT FIRST ENTRY ONLY
        with open(self.uscs_file, "r") as fh:
            for row in fh:
                if "KCNE1" in row:
                    row_split = row.replace(")","").replace("(","").replace("'","")
                    row_split = row_split.split(", ")
                    break
        self.chromosome = row_split[2]
        self.start = int(row_split[3])
        self.stop = int(row_split[4])

        ### GET EXONS ###
        self.exons = row_split[6]

        ### GET READS ###
        self.reads = list(self.alignfile.fetch(self.chromosome, self.start, self.stop))

        ### GET PROPERLY PAIRED READS ###
        self.proper_reads = len([i for i in self.reads if i.is_proper_pair])

        ### GET READS WITH INDELS ###
        # Cigar = Compact Idiosyncratic Gapped Alignment Report
        # cigar operation 1 = insertion
        # cigar operation 2 = deletion
        rd_indel = []
        for i in self.reads:
            if not i.is_unmapped:
                cig = i.cigartuples
                for (operation, length) in cig:
                    if (operation == 1) or (operation == 2):
                        rd_indel.append(i)
        self.reads_indel = len(rd_indel)
        
        ### GET NUMBER OF MAPPED READS ###
        self.reads_mapped = 0
        for i in self.reads:
            if not i.is_unmapped:
                self.reads_mapped += 1

        ### GET TOTAL AVERAGE COVERAGE ###
        # SQ = reference sequence dictionary
        # LN = length of reference sequence
        # SN = name of reference sequence
        chr_len = [i["LN"] for i in self.alignfile.header["SQ"] if i["SN"] == self.chromosome][0]
        self.chr_avg_cov = round(np.mean(self.alignfile.count_coverage(self.chromosome, start=0, stop=chr_len)),2)

        ### GET GENE AVERAGE COVERAGE ###
        self.gene_avg_cov = round(np.mean(self.alignfile.count_coverage(self.chromosome, start=self.start,stop=self.stop)),2)


    def print_summary(self):
        ### GENE SYMBOL ###
        # The official gene symbol approved by the HGNC, which is typically a short form of the gene name.
        # Symbols are approved in accordance with the Guidelines for Human Gene Nomenclature
        print(f"\nGene Symbol: {self.gene}")
        print(f"\nGenomic Region:\n\tLocated on: {self.chromosome}\n\tStart: {self.start}\n\tStop: {self.stop}")

        ### HEADER ###
        print("\nSam Header:")
        for k, v in self.alignfile.header["HD"].items():
            if k == "SO":
                print(f"\tSO (Sorting order of Alignments): {v}")
            if k == "VN":
                print(f"\tVN (Format version): {v}")
            if k == "GO":
                print(f"\tGO: (Grouping of alignments): {v}")

        ### EXONS ###
        print(f"\nNumber of Exons: {self.exons}")

        ### PROPERLY PAIRED READS ###
        print(f"\nProperly Paired Reads: {self.proper_reads}")

        ### READS WITH INDELS ###
        print(f"\nAmount of Reads with Indels: {self.reads_indel}")

        ### MAPPED READS ###
        print(f"\nNumber of Mapped Reads: {self.reads_mapped}")

        ### TOTAL AVERAGE COVERAGE ###
        print(f"\nTotal Average Coverage: {self.chr_avg_cov}")

        ### AVERAGE GENE COVERAGE ###
        print(f"\nAverage Gene Coverage: {self.gene_avg_cov}")

def main():
    print("Assignment 1")
    assignment1 = Assignment1("KCNE1","hg38","USCS_data","chr21.bam")
    print("\nAvailable Attributes of Class 'Assignment1':\n" + [i for i in dir(assignment1) if "__" not in i])
    assignment1.print_summary()
    print("\nDone with assignment 1")

if __name__ == '__main__':
    main()
    
    

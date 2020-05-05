### using bioPython

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
import pandas as pd
import re


# handling sequences
my_dna = Seq("AGTACACTGGTC", generic_dna)
my_dna.complement()
my_rna = my_dna.transcribe()

# will throw errors if there aren't complete codons
my_protein = my_rna.translate()

# similarly, this will throw error
my_protein.complement()

# reading in a GFF file
# file is from wget ftp://ftp.ensembl.org/pub/release-98/gff3/homo_sapiens/Homo_sapiens.GRCh38.98.gff3.gz
col_names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
df = pd.read_csv('Homo_sapiens.GRCh38.98.gff3.gz', compression='gzip',
                         sep='\t', comment='#', low_memory=False,
                         header=None, names=col_names)

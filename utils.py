import pandas as pd
def read_fasta( fasta_path, split_char="!", id_field=0):
    seqs = dict()
    with open( fasta_path, 'r' ) as fasta_f:
        for line in fasta_f:
            if line.startswith('>'):
                uniprot_id = line.replace('>', '').strip().split(split_char)[id_field]
                uniprot_id = uniprot_id.replace("/","_").replace(".","_")
                seqs[ uniprot_id ] = ''
            else:
                seq= ''.join( line.split() ).upper().replace("-","")
                seq = seq.replace('U','X').replace('Z','X').replace('O','X')
                seqs[ uniprot_id ] += seq
    print("\nRead {} sequences.".format(len(seqs)))
    return seqs

def read_smi(smi_path):
    smiles = pd.read_csv(smi_path, sep="\t")
    return smiles["smiles"], smiles["id"]
        
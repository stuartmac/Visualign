
import ast
from Bio import AlignIO
from operator import itemgetter
import pandas as pd
import requests
from varalign import alignments
import argparse


def parse_pdb_xrefs(seq):
    """
    Parse Pfam PDB dbxref annotations from sequence.

    Example: ['PDB; 4K7D B; 399-458;', 'PDB; 4K95 J; 399-458;']
    """
    pdb_xrefs = [x.split()[1:] for x in seq.dbxrefs if x.startswith('PDB;')]
    pdb_mappings = []
    for pdb_id, chain_id, res_range in pdb_xrefs:
        chain_id = chain_id.replace(';', '')
        start, end = map(int, res_range.replace(';', '').split('-'))  # Split start/end
        length = len(range(start, end + 1))  # Calculate length
        pdb_mappings.append((pdb_id, start, end, chain_id, length))
    return pdb_mappings


def chimera_command(pdb_id, start, end, chain_id, marked=None, name=None, template_n=0):
    """
    Construct Chimera command to open a PDB and isolat or mark a specific region.
    """
    # Select command template
    templates = ['open {pdb_id}; select #MODEL_ID:{start}-{end}.{chain_id}; select invert sel; delete sel',
                 'open {pdb_id}; setattr r domain true #MODEL_ID:{start}-{end}.{chain_id}; select #MODEL_ID:.{chain_id}; select invert sel; delete sel']
    template = templates[template_n]

    # Construct command from template and variables
    command = template.format(pdb_id=pdb_id, start=start, end=end, chain_id=chain_id)

    # Add marked attributes; format multi- or single-residue selection
    if isinstance(marked, (list, tuple)):
        marked = [str(x) + '.' + chain_id for x in marked]
        marked = ','.join(marked)
    elif isinstance(marked, str):
        marked += '.' + chain_id
    # Add to command
    if marked:
        command += '; setattr r marked true #MODEL_ID:{}'.format(marked)

    # Add name attribute setting to command
    if isinstance(name, str):
        command += '; setattr M name {} #MODEL_ID'.format(name)

    return command


def uniprot_pdb_query(uniprot_id):
    "Query SIFTS 'best_structures' endpoint."
    url = ''.join([sifts_best, uniprot_id])
    result = requests.get(url)
    return result.json()


def find_overlap(mapping, seq_range):
    "Calculate overlap between Pfam sequence and SIFTS mapped PDB."
    uniprot_resnums = range(*[mapping[k] for k in ('unp_start', 'unp_end')])
    covered = set(seq_range).intersection(set(uniprot_resnums))
    return (mapping['pdb_id'], mapping['chain_id'], len(covered) / float(len(seq_range)))


def uniprot_to_pdb(mapping):
    "Map UniProt residue numbers to PDB residue numbers."
    pdb_resnums = range(*[mapping[k] for k in ('start', 'end')])
    uniprot_resnums = range(*[mapping[k] for k in ('unp_start', 'unp_end')])
    return dict(zip(uniprot_resnums, pdb_resnums))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("alignment", help="Input the stockholm format alignment path to be analyzed.",
                        type=str)
    parser.add_argument("residue", help= "Input the residue you want to analyze.", type=int)
    parser.add_argument("magnification", help="Enter how many armstrongs the resolution of the images are to be.", type=int)
    #parser.add_argument("--modelcount", help="Enter your desired number of models.", type=int)
    args = parser.parse_args()

    alignment_path = args.alignment

    aln = AlignIO.read(alignment_path, 'stockholm')
    print aln



    # pdb_mappings = parse_pdb_xrefs(seq)
    #
    # # Read in columns from file
    # umd_family_info = '/homes/smacgowan/projects/umd_families/columns.csv'
    # column_table = pd.read_csv(umd_family_info,
    #                            index_col=0,
    #                            converters={'columns_pandas':ast.literal_eval})
    # column_table.head(2)
    #
    # alignment_name = alignment_path.split('/')[-1][:7]
    # umd_entry = column_table[column_table['AC'].str.contains(alignment_name)]
    # umd_entry
    #
    # umd_columns = umd_entry.get_value(122, 'columns_pandas')
    # umd_columns
    #
    # itemgetter(*umd_columns)(dict(alignments.index_seq_to_alignment(seq)))
    #
    #
    #
    #
    #
    # # Select example sequence mapped PDB (sort so that choice is first)
    # pdb_mappings.sort(key=lambda x: x[4], reverse=True)  # Sort by length
    # pdb_mappings.sort(key=lambda x: x[0], reverse=True)  # Sort by PDB
    #
    # # Lookup marked columns in example seq
    # marked = itemgetter(*umd_columns)(dict(alignments.index_seq_to_alignment(seq)))
    #
    # # Build command
    # command = chimera_command(*pdb_mappings[0][:-1], marked=marked, template_n=1)
    # command

    umd_columns = [args.residue]
    model = 0  # Initialise model ID
    chimera_script = []
    for seq in aln:
        if any([x.startswith('PDB') for x in seq.dbxrefs]):
            # Get PDB xrefs
            pdb_mappings = parse_pdb_xrefs(seq)
            pdb_mappings.sort(key=lambda x: x[4], reverse=True)  # Sort by length
            pdb_mappings.sort(key=lambda x: x[0], reverse=True)  # Sort by PDB

            # Identify marked residue
            index_dict = dict(alignments.index_seq_to_alignment(seq))
            #marked = itemgetter(*umd_columns)(index_dict)
            marked = index_dict[umd_columns[0]]

            # Write command
            command = chimera_command(*pdb_mappings[0][:-1], marked=marked)
            command = command.replace('MODEL_ID', str(model))  # Substitue model ID placeholder

            chimera_script.append((seq.id, command))
            model += 1

    chimera_script

    seqs_known_structure = zip(*chimera_script)[0]
    seqs_known_structure

    # Get SIFTS "best structure" for a sequence.
    sifts_best = 'http://www.ebi.ac.uk/pdbe/api/mappings/best_structures/'
    example = 'Q9JK66'




    result = uniprot_pdb_query(example)
    mapping = result[example][0]  # Extract first result from SIFTS query
    mapping




    # Find overlaps for all retrieved SIFTS mappings
    seq_range = range(313, 378)  # Seq start/end for example sequence
    overlaps = [find_overlap(mapping, seq_range) for mapping in result[example]]
    overlaps.sort(key=lambda x: x[2], reverse=True)  # Reorder
    overlaps




    uniprot_to_pdb(mapping)

    model = 0  # Initialise model ID
    commands = []
    for seq in aln:
        if any([x.startswith('PDB') for x in seq.dbxrefs]):
            # if seq.id in seqs_known_structure:

            # Lookup PDBe mappings
            uniprot_id = seq.annotations['accession'].split('.')[0]
            try:
                pdb_mappings = uniprot_pdb_query(uniprot_id)[uniprot_id]
            except KeyError:
                print 'skipping {}'.format(uniprot_id)
                continue
            pdb_mappings = [x for x in pdb_mappings if x['experimental_method'] == 'X-ray diffraction']
            if len(pdb_mappings) == 0:
                continue

            # Identify pdb/chain with most overlap
            start, end = [seq.annotations[k] for k in ['start', 'end']]
            seq_range = range(start, end)
            overlaps = [find_overlap(mapping, seq_range) for mapping in pdb_mappings]
            overlaps.sort(key=lambda x: x[2], reverse=True)
            pdb_id, chain_id = overlaps[0][:2]
            # Find corresponding mapping
            for mapping in pdb_mappings:
                if mapping['pdb_id'] == pdb_id and mapping['chain_id'] == chain_id:
                    break

            # Identify start and end residues of PDB coverage
            res_map = uniprot_to_pdb(mapping)
            if start not in res_map.keys():
                start = min(res_map.keys())
            if end not in res_map.keys():
                end = max(res_map.keys())

            # PDB resnums
            pdb_start, pdb_end = [res_map[k] for k in [start, end]]

            # Marked columns
            try:
                marked = itemgetter(*umd_columns)(dict(alignments.index_seq_to_alignment(seq)))  ## umd_columns global
            except KeyError:
                marked = None
            if isinstance(marked, list):
                marked = [x for x in marked if x in res_map.keys()]
            if isinstance(marked, str):
                marked = [marked] if marked in res_map.keys() else []
            # Write command
            model_name = seq.id + '({})'.format(pdb_id)
            command = chimera_command(pdb_id, pdb_start, pdb_end, chain_id, marked, model_name, template_n=1) + '; wait'
            command = command.replace('MODEL_ID', str(model))
            commands.append((seq.id, command))
            model += 1
    commands

    # Add match maker
    chimera_script = list(zip(*commands)[1])
    #mm_script = ['mm #{} #{}; wait'.format('0', str(i+1)) for i in range(len(chimera_script)-1)]
    mm_script = ['mm #{}:/domain #{}:/domain; wait'.format('0', str(i+1)) for i in range(len(chimera_script)-1)]
    #mm_script = ['mm #0 #1-{}; wait'.format(len(chimera_script)-1)]
    chimera_script = chimera_script + mm_script
    chimera_script

    #Add quality of life commands
    chimera_script.write = ("display :/marked")
    chimera_script.write = ("focus :/marked z < {}").format(args.magnification)
    chimera_script.write = ("center :/marked")
    chimera_script.write = ("cofr :/marked")
    chimera_script.write = ("select :/marked; namesel marked")
    #\n might work

    #chimera_script.write = ("display :/marked \n focus :/marked z < {} \n center :/marked \n cofr :/marked \n select :/marked \n namesel marked").format()
    chimera_script.write = ("findhbond selRestrict \"marked & without CA/C1'\"reveal true intermodel false")
    chimera_script.write = ("~modeldisp #1-2")
    chimera_script.write = ("copy file {}.png png").format(seq.id + _ + pdb_id)
    chimera_script.write = ("~modeldisp #0")
    chimera_script.write = ("modeldisp #1")
    chimera_script.write = ("copy file {}.png png").format(seq.id + _ + pdb_id)
    chimera_script.write = ("~modeldisp #1")
    chimera_script.write = ("modeldisp #2")
    chimera_script.write = ("copy file {}.png png").format(seq.id + _ + pdb_id)

    # com_file_name = 'PF01485_chimera_alignment.com'
    com_file_name = 'PF00001_chimera_alignment.com'
    com_file_name = 'PF00104_chimera_alignment.com'

    with open(com_file_name, 'w') as output:
        for line in chimera_script:
            output.write(line)
            output.write('\n')
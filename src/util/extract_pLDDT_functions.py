import gzip
from pathlib import Path
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder


def extract_pLDDT(structure_file):
    """
    Extracts the per-residue pLDDT scores from a PDB or mmCIF file.

    Args:
        structure_file (str): Path to the PDB or mmCIF file.

    Returns:
        list: A list of pLDDT scores, one for each residue with a 'CA' atom.
    """

    structure_file = Path(structure_file)

    # compressed or not
    if structure_file.name.endswith(".gz"):
        f = gzip.open(structure_file, "rt")
    else:
        f = open(structure_file, "r")

    # PDB or mmCIF
    if structure_file.name.rstrip(".gz").lower().endswith(".pdb"):
        parser = PDBParser()
        structure = parser.get_structure(structure_file.stem, f)
    elif structure_file.name.rstrip(".gz").lower().endswith(".cif"):
        parser = MMCIFParser()
        structure = parser.get_structure(structure_file.stem, f)
    else:
        raise ValueError(f"Unknown file format: {structure_file.name}")

    # extract pLDDT
    pLDDT = []
    for residue in structure.get_residues():
        if residue.has_id("CA"):
            pLDDT.append(residue["CA"].bfactor)

    f.close()

    return pLDDT


def extract_pLDDT_pdb_gz(pdb_gz_structure_file):
    """
    Extracts the per-residue pLDDT scores from a PDB or mmCIF file.

    Args:
        structure_file (str): Path to the PDB or mmCIF file.

    Returns:
        list: A list of pLDDT scores, one for each residue with a 'CA' atom.
    """

    structure_file = Path(pdb_gz_structure_file)
    f = gzip.open(structure_file, "rt")
    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, f)

    # extract pLDDT
    pLDDT = []
    for residue in structure.get_residues():
        pLDDT.append(residue["CA"].bfactor)

    f.close()

    return pLDDT


def extract_pLDDT_pdb(pdb_structure_file):
    """
    Extracts the per-residue pLDDT scores from a PDB or mmCIF file.

    Args:
        structure_file (str): Path to the PDB or mmCIF file.

    Returns:
        list: A list of pLDDT scores, one for each residue with a 'CA' atom.
    """

    structure_file = Path(pdb_structure_file)
    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, structure_file)

    # extract pLDDT
    pLDDT = []
    for residue in structure.get_residues():
        pLDDT.append(residue["CA"].bfactor)

    return pLDDT


def extract_protein_seq_pdb_gz(pdb_gz_structure_file):
    """
    Extracts the residue sequence from a pdb.gz file.
    """

    structure_file = Path(pdb_gz_structure_file)
    f = gzip.open(structure_file, "rt")
    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, f)

    ppb = PPBuilder()
    seq = ppb.build_peptides(structure)[0].get_sequence()

    f.close()

    return seq

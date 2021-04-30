# S5296

## About

<font size=4>

**S5296.txt**: the training set for parameterizing PremPS model and it contains 5,296 single-point mutations.

**Mutant structures produced by FoldX.tar.gz** contains the initial protein 3D structures for reverse mutations produced by FoldX. 

</font> 

## Terms in the text files

<font size=4>

PDB Id: The PDB entry of protein.

Mutated Chain: Protein chain with mutation.

Mutation_PDB: The mutation corresponding to the residue numbering found in the Protein Data Bank. The first character is one letter amino acid code for the wild-type residue, the second to penultimate characters indicate residue number, and the final character indicates the mutant amino acid.

UniProt: The UniProt ID of protein.

Mutation_UNP: The mutation corresponding to the residue numbering found in the protein sequence.

Label: 'forward' indicates the forward mutations; 'reverse' indicates the reverse mutations.

DDGexp: Experimental changes of unfolding Gibbs free energy upon mutations (in kcal/mol).

Location: Location of the mutated site in the protein structure (core or surface).

similar proteins: The “similar proteins” of each protein. MMseqs2 software was used to find the “similar proteins”; the sequence identity is set to 25% and the alignment covers at least 50% of query and target sequences.

**The rest of columns are the features described below （see more details in the paper）:**

DCS: the change of conservation after mutation calculated by PROVEAN method. 

DOMH: the difference of hydrophobicity scale between mutant and wild-type residue type.

PSSM: the Position-Specific Scoring Matrix created by PSI-BLAST.

P\_L, P\_FWY and P\_RKDE: the fraction of aromatic residues (F, W and Y), charged residues (R, K, D and E) and leucine (L) buried in the protein core, respectively.

N\_Hydro and N\_Charg: the number of hydrophobic (V, I, L, F, M, W, Y or C) and charged amino acids (R, K, D or E) at 23 sites in the protein sequence, respectively, the center of which is at the mutated site. 

SASA\_pro and SASA\_sol: the solvent accessible surface area (SASA) of the mutated residue in the protein and in the extended tripeptide respectively. 

<font>

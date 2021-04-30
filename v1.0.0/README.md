# PremPS-1.0.0
## About
<font size=4> 
  
Before running PremPS, you need to create a folder including two input files (see the example of 2020100417132606935696574). 
  
</font>

## Two input files in the folder of 2020100417132606935696574
<font size=4> 

1. The 3D structure of a protein (PREMPS.pdb1), which can be obtained from the Protein Data Bank (PDB) or created by the user.

2. The file containing mutation information (2020100417132606935696574.input), who's name must be consistent with the input folder name.

- PDBfile: coordinate file of a protein structure.
- Chains: the selected chains that will be taken into account during the calculation.
- MutChain: the protein chain where the mutation occurs.
- Mutation_PDB: the mutation: the first character is one letter amino acid code for the wild-type residue, the second to penultimate characters indicate   residue number in MutChain, and the final character indicates the mutant amino acid.
- Result_Id: a number defined by the user.
- isPI: equal to one or zero if mutant structure is produced or not for each mutation.

  The columns are separated by tabs.

</font>

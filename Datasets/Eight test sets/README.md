# Eight test sets

## About

<font size=4>

The eight datasets that were taken from the previous studies are used to assess the predictive performance of PremPS and compare with other computational methods (see more details in the paper). 

</font> 

## Terms in the text files

<font size=4>

PDB Id: The PDB entry of protein.

Mutated Chain: Protein chain with mutation.

Mutation_PDB: The mutation corresponding to the residue numbering found in the Protein Data Bank. The first character is one letter amino acid code for the wild-type residue, the second to penultimate characters indicate residue number, and the final character indicates the mutant amino acid.

UniProt: The UniProt ID of protein.

Mutation_UNP: The mutation corresponding to the residue numbering found in the protein sequence.

DDGexp: Experimental changes of unfolding Gibbs free energy upon mutations (in kcal/mol).

Location: Location of the mutated site in the protein structure (core or surface).

Label (in Ssym.txt, S250.txt and S2000.txt): 'forward' indicates the forward mutations; 'reverse' indicates the reverse mutations.

Label (in S350.txt): '1' indicates that the mutations for which the values in change of stability are available for all methods.

PremPS: Predicted changes of unfolding free energy upon mutations by PremPS (in kcal/mol). Positive and negative sign corresponds to destabilizing and stabilizing mutations, respectively. 

PremPS_M: The PremPS model was retrained after removing the overlapped mutations and their corresponding reverse mutations with each test set from the training dataset.

PremPS_P: The PremPS model was retrained after removing all mutations in the “similar proteins” from the training dataset. 

Twenty-fold (in S1925.txt): PremPS score when the model was retrained on S1925 and then preforming 20-fold cross-validation (in kcal/mol).

<font>

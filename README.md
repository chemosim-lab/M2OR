# M2OR : Database of Olfactory Receptors-Molecule Measurments

The sense of smell, is initiated by the activation specific receptors expressed at the
surface of olfactory neuron, induced by the binding of an odorant. This first step activates signaling pathways leading to the perception
of distinct odors. To gain a deeper understanding of the molecular mechanisms involved in
olfaction, we have gathered a database of OR-odorant pairs, drawn from published literature and
public databases. It includes over 51,483 unique pairs of ORs-molecules. The database contains
information about the responses, receptors, molecules, and experiments that have been performed.

## Database Overview

### Molecule
Provided identifier in the source are used to find InChI Keys on PubChem (*Name*, *CID*, 
*CAS*, *InChI Key*). We use InchI Keys as the standard identifier for the database. In case only the structure is available, 
the canonical SMILES *canonicalSMILES* is deduced using openbabel. If the isomerism is not specified, the molecule 
is considered as “sum of isomers” as long as it has at least one chiral center. The count of chiral centers is determined with RDKit. 
Mixtures of odorants are identified as "mixture" and written as a space-separated list of InChI Keys (*Mixture*).

### Receptor
Uniprot IDs are obtained from either the sequence, given identifiers, or when necessary, the name (*Gene Name*, *Uniprot ID*). 
If the Uniprot ID is not available, only the sequence is included (*Sequence*). In cases where the receptor is mutated, 
the wild-type identifier or sequence is retrieved, and the mutation is indicated separately using the AApositionAA format (*Mutation*). 
The taxon name of the family from which the OR come from is also included (*species*).

### Response
Response is described by a binary code: 0 non-agonist, 1 agonist (*Responsive*). All decisions regarding responsiveness 
are made by the original authors, except for Mainland et al. In the case of primary or secondary screening, 
the concentration used (*Value_Screen*) and its unit (*Unit_Screen*) is specified as well as, when available, the raw response value 
(*Value*) and its unit (*Unit*). For dose-response measurements, EC50 value is added in *Value*, *Unit* columns. 
In case the ligand is a non agonist, "n.d" is used to describe the pair as non-responsive and the maximum concentration used in the test is also included.
Moreover, the number of independent experiment is included (*Nbr_measurements*). 

### Bio-Assay
The biological test is described by the quantity measured (*Type*), the OR expression system (*Cell_line, *Gprotein*, *Co_transfection*, *Tag*), 
the experimental conditions (*Delivery*, *Assay*) and finally the tooling (*Assay system*). 

### Source
The reference (*Reference*, *DOI*) which the pairs originated from is given as well as the location (*Reference Position*) (Table, Figure, Supp Mat) 
where the information was found. 

![image](https://user-images.githubusercontent.com/73403769/235442338-85d09711-8130-4c8f-974b-108e759975d3.png)

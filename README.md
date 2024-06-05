# :airplane: TLDR: 
The maintained and up-to-date data can be searched through and downloaded here: [m2or.chemsensim.fr](https://m2or.chemsensim.fr/)

# M2OR : Database of Olfactory Receptors-Molecule Measurments [1]

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

### To cite us
```
@article{10.1093/nar/gkad886,
author = {Lalis, Maxence and Hladiš, Matej and Khalil, Samar Abi and Briand, Loïc and Fiorucci, Sébastien and Topin, Jérémie},
title = "{M2OR: a database of olfactory receptor–odorant pairs for understanding the molecular mechanisms of olfaction}",
journal = {Nucleic Acids Research},
pages = {gkad886},
year = {2023},
month = {10},
abstract = "{Mammalian sense of smell is triggered by interaction between odorant molecules and a class of proteins, called olfactory receptors (ORs). These receptors, expressed at the surface of olfactory sensory neurons, encode myriad of distinct odors via a sophisticated activation pattern. However, determining the molecular recognition spectrum of ORs remains a major challenge. The Molecule to Olfactory Receptor database (M2OR, https://m2or.chemsensim.fr/) provides curated data that allows an easy exploration of the current state of the research on OR-molecule interaction. We have gathered a database of 75,050 bioassay experiments for 51 395 distinct OR-molecule pairs. Drawn from published literature and public databases, M2OR contains information about OR responses to molecules and their mixtures, receptor sequences and experimental details. Users can obtain information on the activity of a chosen molecule or a group of molecules, or search for agonists for a specific OR or a group of ORs. Advanced search allows for fine-grained queries using various metadata such as species or experimental assay system, and the database can be queried by multiple inputs via a batch search. Finally, for a given search query, users can access and download a curated aggregation of the experimental data into a binarized combinatorial code of olfaction.}",
issn = {0305-1048},
doi = {10.1093/nar/gkad886},
url = {https://doi.org/10.1093/nar/gkad886},
eprint = {https://academic.oup.com/nar/advance-article-pdf/doi/10.1093/nar/gkad886/52409607/gkad886.pdf},
}
```

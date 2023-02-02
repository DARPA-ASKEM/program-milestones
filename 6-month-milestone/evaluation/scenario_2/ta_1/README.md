Directory for TA1 delivarables.

## Py-ascets

SIDARTHE-py-ascet_from_ta2.json: is the py-ascet from ta_2 as the gold standard

SIDARTHE-MIT_a.json: is the MIT un-annotated py-ascet from the code version a

SIDARTHE-MIT_a.json: is the MIT un-annotated py-ascet from the code version b

SKEMA_SIDARTHE_PN.json: This is the extraction of the equations to the PN done by SKEMA, it is a full PN. 

SKEMA_SIDARTHE_PN_V.json: This is the extration from the equations to the PN done by SKEMA, for the SIDARTHE_V model. 

## Initial Conditions and Parameters

SIDARTHE-ic-unit1.json: is the initial concentrations for the unit tests

SIDARTHE-params-t#.json: is the parameters for unit test 1 and unit test 2. t1 is  for unit test 1 and the first of the parameter value for unit test 2. The remaining ones are for the time dependence of the parameters for unit test 2. 

## Connection between Artifacts

SIDARTHE-extracted-vars.txt: is the list of variables and their definitions extracted by MIT from the [SIDARTHE paper text](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7175834/pdf/41591_2020_Article_883.pdf) (paper PDF was pre-processed through COSMOS to get raw text).

SIDARTHE-params-dkg.txt: is the connection between the variables extracted from paper text and the DKG terms.

SIDARTHE-latex-params.txt is the connection between the latex extract from paper formula and the variables extracted from paper text.

## Enriched Gromets
The Gromet FNs for the SIDARTHE sources with text extraction metadata attached inside the `gromet_fn` directory

The skema text extractions linked to those gromet fn are found in `scenario2_skema_extractions.xlsx`


## Equation and code alignment
Equation Alignment in the SIDARTHE and SIDARTHE+V papers.xlsx is the alignment result for the core equations in the SIDARTHE and SIDARTHE+V papers.

Equation Alignment in the SIDARTHE paper and code version A.xlsx is the alignment result for the core equations in the SIDARTHE and the implementation in the code version A.

Equation Alignment in the SIDARTHE paper and code version B.xlsx is the alignment result for the core equations in the SIDARTHE and the implementation in the code version B.


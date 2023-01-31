HMS notes:
The scenario calls for finding models that represent hospitalization and deaths. 
Terarium already has ~27 models loaded that were obtained from BioModels via MIRA. These are searchable
through MIRA DKG concepts. We used the Terarium model explorer to identify two models
that are relevant in this case: Ndairou202 (BIOMD0...0958) and Paiva2020 (BIOMD0...0960).

We found that the curators of these models grounded equivalent concepts somewhat inconsistently.
We therefore remapped some groundings to use consistent DKG terms, thereby improving consistency.

These models can be used donwstream using the following exports in this folder:
- scenario3_biomd958.json: Petri net JSON export of the Ndairou2020 model from MIRA
- scenario3_biomd958_mira.json: MIRA template model JSON export of the Ndairou2020 model from MIRA
- scenario3_biomd958_rate_law_subs.json: Petri net JSON export of Ndairou2020 after MIRA's template and rate law simplification and parameter formula simplification (substitute formulas over multiple parameters and introduce a new parameter with the substituted value) to make clean mass-action representations.
- scenario3_biomd960.json: Petri net JSON export of the Paiva2020 model from MIRA
- scenario3_biomd960_mira.json: MIRA template model JSON export of the Paiva2020 model from MIRA
- scenario3_biomd960_rate_law_subs.json: Petri net JSON export of Paiva2020 after MIRA's template and rate law simplification and parameter formula simplification (substitute formulas over multiple parameters and introduce a new parameter with the substituted value) to make clean mass-action representations.

The MIRA structural comparison of the two models (using MIRA's internal
visualization, not the more sophisticated HMI) is available in
sc3_model_958_960_delta_corrected.png.

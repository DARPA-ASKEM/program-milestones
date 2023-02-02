# MIT Products Overview, Feb 1, 2023

Our work comprised several outputs:
1. *Data-generating REST APIs* that generate integrated extractions from descriptive model inputs
2. A user-facing *tool for accepting and editing TA-1 Petri Net approximations* that was written by TA-4 (thanks Manfred and Daniel) and powered by MIT and UA data
3. *An initial experiment* to evaluate our data inside this tool
4. *Various reachback and data services* for helping to manage flawed data sources

---
## Data-Generating REST APIs

We have several REST APIs that start with an XDD query for a paper and yield a range of extractions. The extractions include:
1. Variables extracted from program source
2. Variables and useful values from the paper text
3. Formulas and variables from the formulas (which are parsed from images)

These are tied together, and also grounded against elements in the DKG. The resulting integrated file gives a comprhensive multi-information-channel 
background story about an input model.

This metadata model is intended to be used by a human being inside TA-4 tooling to repair and better understand the recovered models.

Note that we have functionality that links the above data items to structured data columns and was shown during the hackathon, but 
for narrow software engineering reasons is not included in the integrated metadata set here.

This directory includes a PDF of Jupyter output that summarizes the full set of data we can generate (see file `mit-feb1-extraction-summary.pdf`). You can also 
get the Jupyter notebook at https://github.com/mikecafarella/mitaskem/blob/0ddfbbd17732ed9424df135aeb0c1a9dc83fcd2c/demos/2023-02-01/mit-feb1-demo.ipynb

This version of the data includes some data elements from Arizona, but we need to do a bit more integration work to fully exploit what they have to offer.

---
## User-facing tool

To be clear, we didn't write this tool! The tool was written by Manfred and Daniel from TA-4. But the tool is meant to display metadata produced by MIT and UA processes.

You can see a video of the tool in action at this URL:
https://drive.google.com/drive/folders/1VI-iDNR3Km9RP7mjKrBVoV5mT-OGK8V8

It is an extended attempt by James of Florida to use the tool to recover the Petri Net.  More about this shortly.

---
## Initial Experiment

One of the goals of TA-1 is to produce a model that can be used by downstream teams. 


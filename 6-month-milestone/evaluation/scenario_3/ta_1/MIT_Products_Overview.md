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
https://drive.google.com/file/d/1uNLWvTCSKpenkbTnn0jOAuBKvB-D2yBz/view?usp=share_link

It is an extended attempt by James of Florida to use the tool to recover the Petri Net.  More about this shortly.

---
## Initial Experiment

One of the goals of TA-1 is to produce a petri net model that can be used by downstream teams. If we do a better job of extracting information from code, paper, formula, and data sources, then the resulting incomplete-but-annotated petri net should be straightforward for a user to correct.

We can now ship extractions and metadata to the user-facing editing tool described above, so we ran a simple (and flawed) initial experiment.

Justin from UA attempted to build the SIDDARTHE Petri Net using an editor that showed all of the correct formulas, as well as good editing functionality that showed the formula described by the current-Petri-Net-in-progress. It took him about 15 minutes.

We sent our extractions for the same model to a novel editing tool. This tool lacked almost all of the nice features described in the paragraph above. It didn't show the correct formulas, and it didn't show the formulas that describe the Petri-Net-in-progress. This tool didn't even allow the user to delete an edge that was incorrectly placed! However, it did show off all the upstream metadata extractions.  James used this tool to do the same task as Justin. His experience was recorded in the video above.

The bad news is that it took James 25 minutes before he ran out of gas: although he made a lot of progress, the network was still not done. So we failed to "beat" the naive approach in terms of time needed to arrive at a correct Petri Net.

The good news is that James gave a lot of useful feedback, much of it centered around a few small changes we can make to the tool itself rather than the extractions.  One important suggestion was to show extracted species and transitions annotated with actual rendered portions of the original paper  (rather than unrendered extracted LaTeX). Another was to reorganize the way we present information; right now the "per-node-metadata" display made it hard to have a TODO list that he worked his way down.  It was a constant struggle to remember where he was in the Petri Net process.  There are a few moments in the video above where improved extraction quality might have made a differnce, but right now it seems we can improve a lot by more pure tool based edits.

---
## Reachback and Data Services

Finally, we implemented a few services that are not integrated in the objects listed above.

The first is the extraction-to-datacolumn integrator. As mentioned above, this REST service works but is not integrated into the overall metadata extraction flow. It will be shortly.

The second is a "reachback service" for variables that have seemingly-unreliable data. A TA-3 modeler may decide a particular variable-with-data is possibly not reliable --- perhaps the data failed a statistical test or perhaps the user read an article in the news about how certain COVID stats are becoming less reliable with the growth of home testing.  The user should be able to "reach back" for more support data from TA-1.

Right now the system works by taking as input a particular variable name.  This is a set of data the user views as flawed (say, COVID infection numbers). The system then emits as output a list of other known datasets that might be relevant or correlated to the original questionable dataset. For example, the COVID infection data might be a query input, and the system might return datasets that are fairly closely correlated, such as wastewater data. The goal is to provide the user with data that might be useful as a proxy, or in predicting the variable's true value.

You can see an example of the output of this system in the local file `mit-feb1-dataquality-reachback-example.json`

This service works as of the end of Wednesday but has not yet been integrated into any user-facing system



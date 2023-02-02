# TA3 Evaluation Results

Included is a reproducible markdown script which includes the analysis and inline evalution code recreating
the full analysis of Scenario 2 from start to finish. For full reproducibility, we have ensured that it
works on a remote server that runs from scratch with strict unit testing enabled. The results of the build
can be found at:

https://chrisrackauckas.github.io/ASKEM_Evaluation_Staging/dev/Scenario2/Evaluation_Scenario_2/

Full reproducibility instructions are given in:

https://chrisrackauckas.github.io/ASKEM_Evaluation_Staging/dev/

A static render of a Pluto notebook of the results can be loaded via loading the following into a web browser:

https://github.com/DARPA-ASKEM/program-milestones/blob/main/6-month-milestone/evaluation/scenario_2/ta_3/Scenario2.html

A cloud runnable version which allows for changing any plots interactively can be found via signing in at the link:

https://juliahub.com/pluto/editor.html?id=69609761-ea6f-4e63-b3d3-94a7dd8e989b

A locally runnable Pluto notebook can be found as the Scenario2.jl and can be opened using the instructions 
from [Pluto.jl](https://github.com/fonsp/Pluto.jl)

## Using the Pluto Notebooks

The Pluto notebooks are fully relocatable and download the dependencies and data on-demand, making them easy
install and use. To do this, first install Julia v1.8.5 from https://julialang.org/downloads/. Then
open up Julia and run:

```julia
using Pkg
Pkg.add("Pluto")
import Pluto
Pluto.run()
```

When this window opens in your browser, navigate to the Pluto notebook. Opening it will start the download
of all data and dependencies, along with the automatic installation and compilation of all dependencies.
This may take awhile due to the complexity of the downloads.
# TA3 Evaluation Results

Included is a reproducible markdown script which includes the analysis and inline evalution code recreating
the full analysis of Scenario 3 from start to finish. For full reproducibility, we have ensured that it
works on a remote server that runs from scratch with strict unit testing enabled. The results of the build
can be found at:

https://chrisrackauckas.github.io/ASKEM_Evaluation_Staging/dev/Scenario3/Evaluation_Scenario_3/

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
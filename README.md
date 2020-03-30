# Covid-19-Policy-Simulator

This is the R code underlying our preprint ``How will this continue? Modelling interactions between the COVID-19 pandemic and policy responses'' (Axel G. Rossberg and Robert J. Knell, March 2020).

The main code is in the file Covid-19-Policy-Simulator.R. If you simply source/run the file, it will generate a figure showing one sample run of the model with the published standard parameters. You can generate multiple runs at once by increasing the values of the parameter n_runs in the script to e.g. 200 or 2000 (may take a few minutes). The graph generated then overlays the results of all runs in pales colour, and the result of the last run in clear lines.

To do multiple runs without the need to wait for the graphics to rebuild, execute the command 

  plotting <- FALSE
  
once.

The code is heavily annotated. It should be easily to read and modify. Suggestions for improvements are welcome.

The script MakeGraphs.R contains command lines that we used to generate the graphs in Figs. 2-4 in the preprint.


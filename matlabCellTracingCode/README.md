# README #

This repository contains the analyses for the single cell tracing experiment. All code specific to those project should go here. More general code should be kept in ARA_Tools or even stitchIt. It does nice things like this:

![allcells_overlay.png](https://bitbucket.org/repo/n6A9az/images/604005974-allcells_overlay.png)

For details on a selection of specific topics see [the wiki](https://bitbucket.org/lasermouse/single_cell_trace_analysis/wiki/Home)

## Installation

1. Either [download the zip](https://bitbucket.org/lasermouse/single_cell_trace_analysis/downloads) or (better) clone the repository
in a Git client like [SmartGit](http://www.syntevo.com/smartgit/) or [SourceTree](https://www.sourcetreeapp.com). 

2. Add ```./code``` and all its sub-directories to the MATLAB path. Other directories, such as ```./plots``` are intended to contain routines that do things with the data but these routines shouldn't be added to the MATLAB path. So to run these specific plotting 
functions you will ```cd ``` to the appropriate directory and run the code locally from within that. 



## Useful functions
* ```listAllTracedCells``` -- list the paths of all traced cells (either registered to ARA or in sample space)
* You can make an HTML page with lots of cool stuff by going to ```single_cell_code/plots/clustering``` and running  ```publish clusterBasedOnProjectionPatterns_PUB.m```
* You should also browse the code too -- there's more in there.


## Generating plots of single neurons
* `overlayCellOnOutlines` - Produces three image files each showing the three brain orientations with brain area outlines and an traced cell overlaid.  The function saves eps, png, and fig. These can be used to make webpages, etc. 
* `overlayCellOnARA` - A function associated with `overlayCellOnOutlines`. 
This is used by changing to the *experiment* directory of one sample and just calling the function. It will produce plots in an `overlays` sub-directory in the `sample2ARA` directory.


## Generating plots of multiple neurons
* ```overlayCellsOnThreeProjections``` -- makes publication ready plots with the three orientations side by side


## Examples
* ```help(overlayCellsOnThreeProjections)```


## Legacy code
* ```areaOutline``` -- Build outlines that was used for the above plotting functions. 
For new outline system see `help generateProjections`



## Dependencies
If you're installing yourself, note that the toolbox also requires:

* [ARA_tools](https://bitbucket.org/lasermouse/ara_tools) and its dependencies
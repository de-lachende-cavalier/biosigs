##STRUCTURE

The code has been divided in multiple files for modularity and ease of perusal. The MATLAB UI should automatically signal which files are functions and which ones are scripts (I made sure to create the former with the New > Function functionality).

In any case, and for extra clarity, the scripts are "compare_mses.m", "hrv.m", "main.m", "misc.m", and "demo.m":

- compare_mses.m: Compares the mean squared errors from the various sources, used to make sure that we're using the best possible source.
- hrv.m: Implements HRV detection.
- main.m: The "main" file to explore the heart rate and breathing rate estimation.
- misc.m: A file to explore miscellaneous experiments.
- demo.m: A file for playing around with the code.

The "spec.txt" file containes the specification as written on Ariel.

##USAGE

The way to use the code should be straight forward: all the scripts can be run with the classical "Run" button and have also been divided in pieces for ease of use, so that we can run certain sections only with "Run Section".

I have tried to mantain the code as readable and obvious as possible, adding comments where it seemed necessary.


##DATA

The data was obtained using phyphox (https://phyphox.org/), so if you want to try out new data I would suggest using the same application and exporting the resuling data in CSV format (with commas, which should be the second option one gets when one chooses to export data).

For all these experiments, I have created a "blueprint" demo.m file and added an empty subdirectory to data/vert_xiphoid/ for storage, which I suggest you utilise.

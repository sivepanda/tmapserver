# Tissuumaps Setup with H5AD!

## What you'll need
1. Installed ( Tissuumaps )[https://tissuumaps.github.io/installation/]. We'll go over installing Conda on the DCC later for this step.
2. You need an H5AD file on hand to view!


## Running
Running Tissuumaps is fairly simple. In this case, we are setting up for DCC usage, which means that we'll need to both install dependencies, conda, and pip.

To install Conda, follow the guide (here)[https://oit-rc.pages.oit.duke.edu/rcsupportdocs/software/user/#miniconda-installation-sample] under "Miniconda intallation sample". Then, you can follow the directions on the Tissuumaps installation website for Conda installation. 

You then install Tisuumaps with pip using a pipenv in your working directory.

In order to run Tisuumaps on the DCC, you'll need port forwarding to access the core's localhost on your computer:

Run this script in order to Port Forward the locally hosted Tissuumaps instance to your PC
`ssh  -t -t -i .ssh/id  [netid]@dcc-login.oit.duke.edu -L 5000:localhost:5000 ssh -i .ssh/id dcc-core-[corenumber] -L 5000:localhost:5000`

You can also locally set some of these variables if you do not forsee them changing. The NetID is likely the more likely choice. You can do this by:
`export netid=abc123`

And optionally the corenumber the same way, though this changes every time you make a request. Then you can run:

`ssh -t -t -i .ssh/id ${netid}@dcc-login.oit.duke.edu -L 5000:localhost:5000 ssh -i .ssh/id dcc-core-${corenumber} -L 5000:localhost:5000`


Then, running the Tissuumaps instance in the DCC is as simple as tissuumaps_server "/path/to/image"


### Common issues
The #1 issue I faced is what the script in this repo addresses. Tissuumaps is somewhat opinionated in the input it recieves, and certain non escaped characters cause issues with visualizing all data properly. If you see that not all genes are viewable or suspect you have an issue, you can run the sanitization script in this repo by cloning and cd-ing into it, then:

1. `pip install .`
2. `adata-sanitize [input].h5ad [output].h5ad`

#### Other Common Quirks
- Path to image must be a *direct path*
- Firefox (& other Gecko based browser like Zen, Waterfox, etc) appear to have a memory leak and tend to crash when running Tissuumaps

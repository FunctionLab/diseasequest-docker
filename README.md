# diseaseQUEST-docker


## diseaseQUEST setup

First, install docker: https://docs.docker.com/install/. (Make sure docker is running after you install it.)

Then, build and initialize the dq (diseaseQUEST) docker image. Initialization will trigger a download of a smaller ~5GB test data file (all data, etc., would result in a >500GB file), which includes some example tissue networks for disease prediction; all tissue networks can be downloaded from http://wisp.princeton.edu/download/.

```bash
$ git clone git@github.com:FunctionLab/diseasequest-docker.git
$ cd diseasequest-docker; bash setup.sh 2>&1 | tee setup.log
```
Note that if you are running docker on a mac or windows machine, a VM with a 2GB memory limit is created for docker. Since running network integration needs ~2GB of memory, and running disease prediction requires ~2GB of memory, this 2GB allocation is easily exceeded (e.g., running more than one process at a time), and processes will be killed. To avoid this, the memory allocation for docker should be increased (preferably to at least 4GB).

For details on how to increase the memory allocation, see the "Advanced" sections for your corresponding OS: 
* Mac: https://docs.docker.com/docker-for-mac/#advanced
* Windows: https://docs.docker.com/docker-for-windows/#advanced

We have also provided precompiled binaries (`bin/`) for relevant sleipnir functions (code in `src/sleipnir`).

## Run commands for specific diseaseQUEST results

### Predict tissue-specific networks using semi-supervised regularized Bayesian integration

This will take ~1hr:

```bash
docker run -v `pwd`/data/:/dq/data -v `pwd`/outputs:/dq/outputs dq networks [tissue name]
```

Example:

```bash
docker run -v `pwd`/data/:/dq/data -v `pwd`/outputs:/dq/outputs dq networks dopaminergic_neuron
``` 
Will generate the dopaminergic neuron network located at ./output/all/predictions.

The network is a .dab file (gene ids are entrez gene ids) that can be queried and explored using `Dat2Dab` (in `bin/`); more documentation for `Dat2Dab` is available at https://libsleipnir.bitbucket.io/Dat2Dab.html.

You can start a bash docker to interact with the predicted network interactively:
Example for finding the entire network neighborhood of bcat-1 (entrez gene id 181423) or querying the edge weight between bcat-1 and ifc-2 (worm ortholog to human NEFL gene, one of bcat-1's top neighbors, which has been associated with neurodegenerative disorders):
```bash
docker run -v `pwd`/bin:/home/bin/ -v `pwd`/outputs/:/home/outputs/ -it --rm bash
$ cd home
$ bin/Dat2Dab -i outputs/all/predictions/dopaminergic_neuron_weights.dab -l 181423
$ bin/Dat2Dab -i outputs/all/predictions/dopaminergic_neuron_weights.dab -l 181423 -L 180414
```

Alternatively, to output the entire network as a .dat file (non-binarized file with a list of all edge weights):
```bash
docker run -v `pwd`/bin:/home/bin/ -v `pwd`/outputs/:/home/outputs/ -it --rm bash -c '/home/bin/Dat2Dab -i /home/outputs/all/predictions/dopaminergic_neuron_weights.dab -o /home/outputs/all/predictions/dopaminergic_neuron_weights.dat'
```

### Predict disease-gene associations

```bash
docker run -v `pwd`/data/:/dq/data -v `pwd`/outputs:/dq/outputs dq predictions [disease name] [tissue name]
```

Example:

```bash
docker run -v `pwd`/data/:/dq/data -v `pwd`/outputs:/dq/outputs dq predictions parkinsons_disease dopaminergic_neuron
```

Will generate the predictions for Parkinson's disease using the dopaminergic neuron network located at ./output/parkinsons_disease.predictions.


```bash
docker run -v `pwd`/data/:/dq/data -v `pwd`/outputs:/dq/outputs dq predictions longevity intestine
```

Will generate the predictions for longevity using the intestine network located at ./output/longevity.predictions.

### Help menu

```bash
docker run dq -h
```

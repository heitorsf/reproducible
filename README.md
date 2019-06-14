# A model of the neuromuscular system

This repository contains a computational model of the human neuromuscular system implemented using the NEURON simulator (https://neuron.yale.edu) with Python and the module NetPyNE.

## Environment

These are the requirements for running the model:
1. Python 2.7 (https://www.python.org/)
2. NEURON 7.4 (https://neuron.yale.edu)
3. NetPyNE 0.7.0 (https://netpyne.org)
4. The files in this repository.

You can either install the requirements manually or **run in docker container**.

## Using the docker container

First, install Docker following the instructions on https://docs.docker.com/install/.

You will use the *docker image* [heitorsf/nerlab:reproduce] (https://hub.docker.com/r/heitorsf/nerlab) to create a *container*. You can name the container for further use. For that, type in the command line to create a container using the docker image

```
docker run -it -p 8888:8888 --name my_container heitorsf/nerlab:reproduce
```

A URL will be prompted in your screen, copy and paste it to your browser.

Once in the jupyter notebook, open `/work/reproducible/deliver/Artigo_Executavel.ipnb`.

To close everything, go to the command line screen and type Ctrl+C.

If you want to run your container again, use:

```
docker start -i my_container
```

For support, please e-mail heitorsf@gmail.com.
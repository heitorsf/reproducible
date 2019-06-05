## A model of the neuromuscular system

This repository contains a computational model of the human neuromuscular system implemented using the NEURON simulator (https://neuron.yale.edu) with Python and the module NetPyNE.

## Environment

These are the requirements for running the model:
1. Python 2.7 (https://www.python.org/)
2. NEURON 7.4 (https://neuron.yale.edu)
3. NetPyNE 0.7.0 (https://netpyne.org)
4. The files in this reposotory.

You can either install the requirements manually or **run in docker container**.

In either case, download or clone this repository first.

# Using the docker container

First, install Docker following the instructions on https://docs.docker.com/install/.

You will use the container [heitorsf/nerlab:netpyne](https://hub.docker.com/r/heitorsf/nerlab). For that, type in the command line:

```
docker run -it -p 8888:8888 -v /YOUR_PATH_TO_THE_CLONED_REPO/reproducible:work/reproducible heitorsd/nerlab:netpyne
```
A URL will be prompted in your screen, copy and paste it to your browser.

Once in the jupyter notebook, open `/work/reproducible/deliver/Executable_Paper.ipnb`.


For support, please e-mail heitor@gmail.com.

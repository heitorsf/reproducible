FROM heitorsf/nerlab:reproduce
WORKDIR /work
USER root
EXPOSE 8888
USER neuron
CMD jupyter notebook

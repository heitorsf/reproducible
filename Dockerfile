FROM heitorsf/nerlab:reproduce

# create user with a home directory
ARG NB_USER=neuron
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /work/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}

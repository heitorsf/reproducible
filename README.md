# Modelo do sistema neuromuscular

Este repositório contém um modelo computacional simplificado do sistema neuromuscular humano. Foi implementado usando o simulador NEURON com Python e com o módulo NetPyNE.

## Ambiente

Estes são as dependências para executar o modelo:
1. Python 2.7 (https://www.python.org/)
2. NEURON 7.4 (https://neuron.yale.edu)
3. NetPyNE 0.7.0 (https://netpyne.org)
4. Os arquivos neste repositório.

Você pode executar o modelo usando um **container docker**, descrito abaixo, ou instalar todas as dependências manualmente.

## Usando o container docker

Primeiro, instale o Docker seguindo as instruções em https://docs.docker.com/install/, de acordo com seu sistema operacional (na página, escolha o sistema pelo menu lateral esquerdo).

Você irá usar a imagem docker chamada [heitorsf/nerlab:reproduce] (https://hub.docker.com/r/heitorsf/nerlab) para criar um *container* (o link é apenas para referência ao Docker Hub). Você pode dar um nome a esse container para utilizá-lo novamente depois. Para isso, use o seguinte comando em um console ou terminal:

```
docker run -it -p 8888:8888 --name my_container heitorsf/nerlab:reproduce
```

*Em Linux, pode haver problemas de permissão, neste caso use `sudo docker (...)`.*

Uma URL irá aparecer em sua tela, copie-a em um *web browser* para acessar o *jupyter notebook*.

Uma vez no jupyter notebook, acesse o diretório `reproducible/deliver` e ali abra o arquivo: `Artigo_Executavel.ipnb`.

Para fechar o container, pressione Ctrl+c no terminal ou console a partir do quel você executou o container.

Se você quiser rodar o mesmo container novamente, use:

```
docker start -i my_container
```

Para suporte, envie um e-mail para heitorsf@gmail.com.

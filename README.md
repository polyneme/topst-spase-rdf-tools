# SPASE RDF Tools

Toolset to produce SPASE RDF and explore teh resulting Knowledge Graph.

## The SPASE Knowledge Graph

The SPASE Knowledge is composed of two main parts:

1. The SPASE ontology, which is an automatically generated OWL Ontology using the SPASE Base Model XSD file available [here](https://spase-group.org/data/schema/spase-2.6.0.xsd). The ontology generation algorithm takes every entity on the SPASE XSD file and turns it into an OWL Class, all the relationships between entities get mapped to owl:ObjectProperties and every literal property of each entity gets mapped to owl:DataTypeProperty, all properties get assigned their corresponding domain and range.
2. SPASE RDF Individuals Data, which is an automatically generate TTL file containing RDF that represents the different SPASE resources on the XML files provided by [hpde](https://github.com/hpde/hpde.io/). This RDF complies with the SPASE Ontology.


## RDF Exploration
### Running on Docker with Docker compose
#### Requirements

- [Git](https://git-scm.com/download/mac)
- [Docker](https://docs.docker.com/get-docker/)

#### Get the code

- Clone this repo and its submodules:
  ```shell 
    git clone --recurse-submodules -j8 https://github.com/alexgarciac/spase-rdf-tools.git
    cd spase-rdf-tools
  ```

#### Decompress pre-processed data

Decompress the pre-processed data under:
`spase-rdf-tools/data/spase.ttl.zip`:

  ```shell 
    cd data
    unzip spase.ttl.zip
    cd ..
  ```

#### Setup

```shell
docker-compose build
```

#### Execution
```shell
docker-compose up
```

Open:

-  Graph Explorer: https://localhost/explorer
-  Jupyter notebook: Check `docker-compose`out for a link like http://127.0.0.1:8888/?token=<token\>. Then Open the browser on the given link and open the `SPASE RDF Exploration.ipynb` file.
-  Fuseki: http://localhost:3030

### Running (almost) without Docker

#### Requirements

- Java
- Docker
- Python 3.8+

#### Setup

1. Install and run Fuseki:
   - Download the latest version of Jena Fuseki from [here](https://jena.apache.org/download/). You can use:
      ```
      cd ~/Applications/
      curl https://archive.apache.org/dist/jena/binaries/apache-jena-fuseki-<fuseki_version>.zip -o apache-jena-fuseki-<fuseki_version>.zip
      ```
     
      Just replace `<fuseki_version>` with the latest version available:
      ```
      curl https://archive.apache.org/dist/jena/binaries/apache-jena-fuseki-4.9.0.zip -o apache-jena-fuseki-4.9.0.zip
      ```
   - Unzip the Jena Fuseki package:
      ````shell
      unzip apache-jena-fuseki-4.9.0.zip
      ````
   - Run Fuseki Server:
      ```shell
      cd apache-jena-fuseki-4.9.0
      ./fuseki-server
      ```
   - Open your browser to check your Fuseki is up and running: http://localhost:3030
   

2. Load the RDF data into a Fuseki dataset:
   - Create a dataset by opening Fuseki on a browser and click on the `add one link`:
   ![Fuseki dataset creation step 1](docs/img/Screenshot%202024-02-17%20at%2011.52.22.png)
   - Name the dataset `spase` and select the dataset type (choose persistent if you plan to re-use this dataset on future runs) and then click create dataset:
   ![Fuseki dataset creation step 2](docs/img/Screenshot%202024-02-17%20at%2011.54.16.png)
   - Upload the pre-processed data, click on add data > select files and select the spase.owl and spase.ttl files under [data](data), then click on upload all:
   ![Fuseki dataset creation step 3](docs/img/Screenshot%202024-02-17%20at%2011.57.09.png)
   - Your new Fuseki dataset should be available under  http://localhost:3030/spase

3. Install and run the RDF Exploration Jupyter notebook:
   - Go to this repo directory:
      ```shell
      cd ~/git/spase-rdf-tools/ # replace with the right location
      ```
   - Get into the python package for the RDF Tools:
      ```shell
     cd spase_rdf_tools
      ```
   - [Optional] Create and activate a virtual environment:
      ```shell
     python3 -m venv venv
     source venv/bin/activate
      ```
   - Install python requirements:
     ```shell
        pip install -r requirements.txt
     ```
   - Setup Jupyter extensions for [KG Exploration](https://github.com/aws/graph-notebook/):
     ```shell
      jupyter nbextension enable  --py --sys-prefix graph_notebook.widgets
      python -m graph_notebook.static_resources.install
      python -m graph_notebook.nbextensions.install
      python -m graph_notebook.ipython_profile.configure_ipython_profile
     ```
   - Run the Jupyter notebook with the extensions:
     ```shell
        jupyter notebook --NotebookApp.kernel_manager_class=notebook.services.kernels.kernelmanager.AsyncMappingKernelManager --ip 0.0.0.0 ./
     ```
   - Open the Jupyter notebook URL and navigate to the [SPASE RDF Exploration](./spase_rdf_tools/SPASE%20RDF%20Exploration.ipynb) notebook.
   - For more information on the graph-notebook please check [their repository](https://github.com/aws/graph-notebook/).

4. Install and run [graph-explorer](https://github.com/aws/graph-explorer):
   - Clone the **graph-explorer** repo (a copy of the repo is included here as a git submodule):
      ```shell
        git clone https://github.com/aws/graph-explorer/
      ```
   - Navigate to the **graph-explorer** directory:
     ```shell
       cd  graph-explorer
     ```
   - Build the Docker image:
     ```shell
      docker build -t graph-explorer .
     ```
   - Start the Docker container:
    ```shell
      docker run -p 80:80 -p 443:443 --env HOST=localhost graph-explorer
    ````
   - Go to **graph-explorer** in your browser by opening https://localhost/explorer (Click on `Advanced > Proceed to localhost` if prompted):
   ![Fuseki dataset creation step 3](docs/img/Screenshot%202024-02-17%20at%2012.20.51.png)
   - Add a new connection by clicking on plus sign in the top right corner, name the connection `spase`, choose `RDF (Resource Description Framework) - SPARQL` as Graph type, and set the endpoint value to `http://localhost:3030/spase`:
   ![Fuseki dataset creation step 3](docs/img/Screenshot%202024-02-17%20at%2012.24.12.png)
   - Synchronise your connection and navigate your graph.
   - For more information on [graph-explorer](https://github.com/aws/graph-explorer), please check their repository.
   
   

## RDF Generation

### Requirements

- Python 3.9.0+

### Install dependencies

- Install python dependencies:
    ```shell
    pip install -r requirements.txt
    ```
### Commands available
```shell
python3 spase_rdf_tools.py --help

Usage: spase_rdf_tools.py [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  create-owl             Creates OWL Ontology using python module
  create-python-model    Creates Python model from XSD file using xsdata
  download-hpde          Downloads and decompress HPDE files from GitHub...
  download-spase-schema  Downloads SPASE XSD schema file from spase-group...
  generate-rdf           Creates TTL RDf File using python module to lo...
```

This is also available as a Jupyter notebook under `spase_rdf_tools/SPASE RDF Generation.ipynb`.
# SPASE RDF Tools

Toolset to produce SPASE RDF and explore teh resulting Knowledge Graph.

## The SPASE Knowledge Graph

The SPASE Knowledge is composed of two main parts:

1. The SPASE ontology, which is an automatically generated OWL Ontology using the SPASE Base Model XSD file available [here](https://spase-group.org/data/schema/spase-2.6.0.xsd). The ontology generation algorithm takes every entity on the SPASE XSD file and turns it into an OWL Class, all the relationships between entities get mapped to owl:ObjectProperties and every literal property of each entity gets mapped to owl:DataTypeProperty, all properties get assigned their corresponding domain and range.
2. SPASE RDF Individuals Data, which is an automatically generate TTL file containing RDF that represents the different SPASE resources on the XML files provided by [hpde](https://github.com/hpde/hpde.io/). This RDF complies with the SPASE Ontology.


The code to download and process the data is available as a Command Line Client `spase-rdf-tools/spase_rdf_tools.py`. You can run:

```
python spase_rdf_tools/spase_rdf_tools.py --help
```

To see the available commands.

This is also available as a Jupyter notebook under `spase_rdf_tools/SPASE RDF Generation.ipynb`.


## Explore RDF

### Requirements

- [Docker](https://docs.docker.com/get-docker/)
- [Git](https://git-scm.com/download/mac)

### Clone repositories

- Clone this repo and its submodules:
  ```bash 
    git clone --recurse-submodules -j8 https://github.com/alexgarciac/spase-rdf-tools.git
    cd spase-rdf-tools
  ```

### Decompress pre-processed data

Decompress the pre-processed data under:
`spase-rdf-tools/data/spase.ttl.zip`:

  ```bash 
    cd data
    unzip spase.ttl.zip
    cd ..
  ```

### Setup

```bash
docker-compose build
```

### Execution
```bash
docker-compose up
```

Open:

-  Graph Explorer: https://localhost/explorer
-  Jupyter notebook: Check `docker-compose`out for a link like http://127.0.0.1:8888/?token=<token\>. Then Open the browser on the given link and open the `SPASE RDF Exploration.ipynb` file.
-  Fuseki: http://localhost:3030


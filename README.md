# SPASE RDF Tools

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

### Download pre-processed data

Decompress the pre-processed data under:
`spase-rdf-tools/data`


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
-  Jupyter notebook: Check `docker-compose`out for a link like http://127.0.0.1:8888/?token=<token>ç
-  Fuseki: http://localhost:3030


# SPASE RDF Tools

## Explore RDF

### Requirements

- [Docker](https://docs.docker.com/get-docker/)
- [Git](https://git-scm.com/download/mac)

### Clone repositories

- Clone this repo:
  ```bash 
    git clone https://github.com/alexgarciac/spase-rdf-tools.git && cd spase-rdf-tools
  ```
- Clone the Graph explorer repo:
  ```bash 
    git clone https://github.com/aws/graph-explorer.git
  ```

### Download pre-processed data




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


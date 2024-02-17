import click

from utils.spase_to_rdf import create_python_model_from_xsd, create_owl_from_python_module, xml_to_rdf
from utils.data_download import download_hpde_files, download_spase_schema_xsd


@click.command()
@click.option('--branch', '-b', default='master', help='Branch name of the HPDE repository')
@click.option('--output-dir', '-o', default='.', help='Output folder to save the downloaded HPDE files')
def download_hpde(branch: str, output_dir: str):
    """Downloads and decompress HPDE files from GitHub repository"""
    download_hpde_files(branch, output_dir)


@click.command()
@click.argument('version')
@click.option('--output-dir', '-o', default='.', help='Output directory to save the downloaded SPASE schema XSD file')
def download_spase_schema(version: str, output_dir: str):
    """Downloads SPASE XSD schema file from spase-group website"""
    download_spase_schema_xsd(version, output_dir)


@click.command()
@click.argument('input-file')
@click.option('--output-mod', '-o', default='spase_model',
              help='Output python module to save the generated Python Model')
def create_python_model(input_file: str, output_mod: str):
    """Creates Python model from XSD file using xsdata"""
    create_python_model_from_xsd(input_file, output_mod)


@click.command()
@click.argument('python_module')
@click.option('--output-file', '-o', default='spase.owl', help='Output OWL file to save the generated Ontology')
def create_owl(python_module: str, output_file: str):
    """Creates OWL Ontology using python module"""
    create_owl_from_python_module(python_module, output_file)


@click.command()
@click.argument('input_dir')
@click.option('--output-dir', '-o', default='data', help='Output directory to save the generated RDF triples')
@click.option('--model-module', '-m', default='spasE_model', help='Model to use on the XML load')
@click.option('--partition', '-p', default='1', help='Number of output files')
def generate_rdf(input_dir: str, model_module: str, output_dir: str, partition: int):
    """Creates TTL RDf File using python module to lo XML files"""
    xml_to_rdf(input_dir, model_module, output_dir, partition)


@click.group()
def cli():
    pass


cli.add_command(download_hpde)
cli.add_command(download_spase_schema)
cli.add_command(create_python_model)
cli.add_command(create_owl)
cli.add_command(generate_rdf)

if __name__ == '__main__':
    cli()

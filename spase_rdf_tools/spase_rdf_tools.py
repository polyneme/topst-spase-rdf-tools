import click

from utils.data_download import download_hpde_files, download_spase_schema_xsd

@click.command()
@click.option('--branch', '-b', default='master', help='Branch name of the HPDE repository')
@click.option('--output-dir', '-o', default='.', help='Output folder to save the downloaded HPDE files')
def download_hpde(branch: str, output_dir: str):
    download_hpde_files(branch, output_dir)


@click.command()
@click.argument('version')
@click.option('--output-dir', '-o', default='.', help='Output directory to save the downloaded SPASE schema XSD file')
def download_spase_schema(branch: str, output_dir: str):
    download_hpde_files(branch, output_dir)
@click.group()
def cli():
    pass

cli.add_command(download_hpde)
cli.add_command(download_spase_schema)

if __name__ == '__main__':
    cli()
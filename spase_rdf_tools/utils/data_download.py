import click
import requests
import os
from zipfile import ZipFile
import io

def download_hpde_files(branch, output):
    """Downloads HPDE files from GitHub repo"""
    url = f"https://github.com/hpde/hpde.io/archive/refs/heads/{branch}.zip"

    # Send a GET request to download the zip file
    response = requests.get(url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Extract the zip file content
        with ZipFile(io.BytesIO(response.content)) as zip_file:
            # Extract all contents to the output folder
            zip_file.extractall(output)

        print(f"Downloaded HPDE files from branch '{branch}' to: {output}")
    else:
        print(f"Failed to download HPDE files. Status code: {response.status_code}")

def download_spase_schema_xsd(version: str, output_dir: str):
    """Downloads SPASE XSD schema file from spase-group website"""
    schema_url = f"https://spase-group.org/data/schema/spase-{version}.xsd"

    # Send a GET request to download the schema file
    response = requests.get(schema_url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Get the file name from the URL
        file_name = os.path.basename(schema_url)

        # Save the schema file to the specified output directory
        output_path = os.path.join(output_dir, file_name)
        with open(output_path, "wb") as f:
            f.write(response.content)

        print(f"Downloaded latest SPASE schema XSD to: {output_path}")
    else:
        print(f"Failed to download SPASE schema XSD. Status code: {response.status_code}")


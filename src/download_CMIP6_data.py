'''
Below is the draft of script to download CMIP6 precipitation data (historical, RCP4.5, 7.0, 8.5, and historical)
under different climate models from NASA NCCS server.'''
import os
import requests
from bs4 import BeautifulSoup

# Base URL for the data
base_url = "https://ds.nccs.nasa.gov/thredds/fileServer/AMES/NEX/GDDP-CMIP6"

# Models and scenarios to download data from
models = [
    "ACCESS-CM2", "ACCESS-ESM1-5", "BCC-CSM2-MR", "CESM2", "CESM2-WACCM", 
    "CMCC-CM2-SR5", "CMCC-ESM2", "CNRM-CM6-1", "CNRM-ESM2-1", "CanESM5", 
    "EC-Earth3", "EC-Earth3-Veg-LR", "FGOALS-g3", "GFDL-CM4", "GFDL-CM4_gr2", 
    "GFDL-ESM4", "GISS-E2-1-G", "HadGEM3-GC31-LL", "HadGEM3-GC31-MM", "IITM-ESM", 
    "INM-CM4-8", "INM-CM5-0", "IPSL-CM6A-LR", "KACE-1-0-G", "KIOST-ESM", 
    "MIROC-ES2L", "MIROC6", "MPI-ESM1-2-HR", "MPI-ESM1-2-LR", "MRI-ESM2-0", 
    "NESM3", "NorESM2-LM", "NorESM2-MM", "TaiESM1", "UKESM1-0-LL"
]
scenarios = ["historical", "ssp245", "ssp370", "ssp585"]
run = "r1i1p1f1"
variable = "pr"

# Function to download a file from a URL
def download_file(url, save_path):
    response = requests.get(url)
    if response.status_code == 200:
        with open(save_path, 'wb') as file:
            file.write(response.content)
        print(f"Downloaded {url} to {save_path}")
    else:
        print(f"Failed to download {url}")

# Function to get file list from HTML page
def get_file_list(html):
    soup = BeautifulSoup(html, 'html.parser')
    files = [a['href'] for a in soup.find_all('a') if a['href'].endswith('.nc')]
    return files

# Create directories and download files
for model in models:
    for scenario in scenarios:
        path = f"{base_url}/{model}/{scenario}/{run}/{variable}/"
        file_list_url = f"{base_url}/{model}/{scenario}/{run}/{variable}/catalog.html"
        
        # Get the HTML content of the file list page
        response = requests.get(file_list_url)
        if response.status_code == 200:
            files = get_file_list(response.text)
            for file_name in files:
                file_url = f"{path}{file_name}"
                save_path = os.path.join("data", model, scenario, run, variable, file_name)
                os.makedirs(os.path.dirname(save_path), exist_ok=True)
                download_file(file_url, save_path)
        else:
            print(f"Failed to access {file_list_url}")

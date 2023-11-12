import urllib.request
import os
import zipfile
dataset_link = 'https://figshare.com/ndownloader/files/38597144'
dataset_file = 'datasets.zip'
output_dir = 'simulation'
download_path = os.path.join(output_dir, dataset_file)

print("Downloading dataset file")
urllib.request.urlretrieve(dataset_link, download_path)
print("Downloaded dataset file")
print("Extracting dataset file")
with zipfile.ZipFile(download_path, 'r') as zip_ref:
    zip_ref.extractall(output_dir)
os.remove(download_path)
print("Exctracted dataset file")
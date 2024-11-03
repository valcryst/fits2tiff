## fits2tiff
Simple fits to tiff converter  

### Requirements:
```
pip install astropy
pip install pillow
```
### Usage:
```
usage: fits2tiff.py [-h] [-i INPUT] [-o OUTPUT] [-d DIRECTORY] -m {rgb,yuyv,mono}
```
  
Use -m to define the input fits format {rgb,yuyv,mono}  
  
Single file convert:  
```
python fits2tiff.py -m rgb -i input.fits -o output.tiff
```
Batch convert:
```
python fits2tiff.py -m rgb -d F:\Capture\M31
```

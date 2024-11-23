## fits2tiff
Simple fits to tiff converter  

### Requirements:
```
pip install astropy
pip install pillow
```
### Usage:
```
usage: fits2tiff.py [-h] [-i File] [-o File] [-d Path] -m {mono8,mono16,rgb,rggb,lrgb,yuyv422,auto}

Konvertiert FITS-Dateien in TIFF-Dateien.

options:
  -h, --help            show this help message and exit
  -i File, --input File
                        Pfad zur FITS-Eingabedatei
  -o File, --output File
                        Pfad zur TIFF-Ausgabedatei
  -d Path, --directory Path
                        Pfad zum Verzeichnis mit FITS-Dateien
  -m {mono8,mono16,rgb,rggb,lrgb,yuyv422,auto}, --mode {mono8,mono16,rgb,rggb,lrgb,yuyv422,auto}
                          mono8    (8-Bit-Graustufen)
                          mono16   (16-Bit-Graustufen)
                          rgb      (RGB)
                          rggb     (Bayer-Pattern)
                          lrgb     (Luminance + RGB)
                          yuyv422  (YUYV 4:2:2)
                          auto     (automatische Erkennung)
```
Single file convert:  
```
python fits2tiff.py -m rgb -i input.fits -o output.tiff
```
Batch convert:
```
python fits2tiff.py -m rgb -d F:\Capture\M31
```

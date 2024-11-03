import argparse
import os
from astropy.io import fits
from PIL import Image
import numpy as np

def yuyv_to_rgb(yuyv_data):
    """Konvertiert YUYV 4:2:2 in RGB."""
    h, w = yuyv_data.shape[-2:]
    rgb_image = np.zeros((h, w, 3), dtype=np.uint8)

    for i in range(0, h * w, 2):
        y1 = yuyv_data[i]
        u = yuyv_data[i + 1] - 128
        y2 = yuyv_data[i + 2]
        v = yuyv_data[i + 3] - 128

        rgb_image[i] = yuv_to_rgb(y1, u, v)
        rgb_image[i + 1] = yuv_to_rgb(y2, u, v)

    return rgb_image

def yuv_to_rgb(y, u, v):
    """Wandelt YUV-Pixelwert in RGB um."""
    r = y + 1.402 * v
    g = y - 0.344136 * u - 0.714136 * v
    b = y + 1.772 * u
    return np.clip([r, g, b], 0, 255).astype(np.uint8)

def fits_to_tiff(fits_path, output_path, mode):
    with fits.open(fits_path) as hdul:
        image_data = hdul[0].data

    if mode == "rgb":
        if image_data.ndim == 3 and image_data.shape[0] == 3:
            r = (255 * (image_data[0] - np.min(image_data[0])) / np.ptp(image_data[0])).astype(np.uint8)
            g = (255 * (image_data[1] - np.min(image_data[1])) / np.ptp(image_data[1])).astype(np.uint8)
            b = (255 * (image_data[2] - np.min(image_data[2])) / np.ptp(image_data[2])).astype(np.uint8)
            rgb_image = np.stack((r, g, b), axis=-1)
            image = Image.fromarray(rgb_image, mode="RGB")
        else:
            raise ValueError("Das Bildformat passt nicht zu RGB (erwartet 3 Kanäle).")
    elif mode == "yuyv":
        if image_data.ndim == 3 and image_data.shape[0] == 3:
            rgb_image = yuyv_to_rgb(image_data)
            image = Image.fromarray(rgb_image, mode="RGB")
        else:
            raise ValueError("Das Bildformat passt nicht zu YUYV.")
    elif mode == "mono":
        # Falls nur ein Kanal verwendet werden soll, nehme die erste Ebene oder das Graustufenbild
        if image_data.ndim == 3:
            image_data = image_data[0]
        image_data = (255 * (image_data - np.min(image_data)) / np.ptp(image_data)).astype(np.uint8)
        image = Image.fromarray(image_data, mode="L")  # "L" für Graustufen
    else:
        raise ValueError("Ungültiger Modus. Verwende -rgb, -yuyv oder -mono.")

    image.save(output_path, format="TIFF")
    print(f"Konvertiert: {fits_path} -> {output_path}")

def convert_directory(directory, mode):
    """Konvertiert alle FITS-Dateien in einem Verzeichnis in TIFF-Dateien."""
    for filename in os.listdir(directory):
        if filename.lower().endswith(".fits"):
            fits_path = os.path.join(directory, filename)
            output_path = os.path.join(directory, filename.rsplit(".", 1)[0] + ".tiff")
            fits_to_tiff(fits_path, output_path, mode)

def main():
    parser = argparse.ArgumentParser(description="Konvertiert FITS-Dateien in TIFF-Dateien.")
    parser.add_argument("-i", "--input", type=str, help="Pfad zur FITS-Eingabedatei")
    parser.add_argument("-o", "--output", type=str, help="Pfad zur TIFF-Ausgabedatei")
    parser.add_argument("-d", "--directory", type=str, help="Pfad zum Verzeichnis mit FITS-Dateien")
    parser.add_argument("-m", "--mode", choices=["rgb", "yuyv", "mono"], required=True,
                        help="Gibt das Format der FITS-Datei an: rgb, yuyv oder mono (Graustufen)")

    args = parser.parse_args()

    if (args.input or args.output) and args.directory:
        parser.error("Die Optionen -i und -o können nicht zusammen mit -d verwendet werden.")
    
    if args.input and args.output:
        fits_to_tiff(args.input, args.output, args.mode)
    elif args.directory:
        convert_directory(args.directory, args.mode)
    else:
        parser.error("Bitte entweder -i und -o für Einzeldateien oder -d für ein Verzeichnis angeben.")

if __name__ == "__main__":
    main()

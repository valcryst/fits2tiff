import argparse
import os
import time
from astropy.io import fits
from PIL import Image
import numpy as np

def preprocess_image_data(image_data, header):
    """Verarbeitet die Bilddaten basierend auf Header-Informationen."""
    bscale = header.get('BSCALE', 1.0)  # Skalierungsfaktor (Standardwert: 1.0)
    bzero = header.get('BZERO', 0.0)    # Offset-Wert (Standardwert: 0.0)

    # Wende Skalierung und Offset an
    image_data = (image_data * bscale) + bzero

    # Normiere Daten auf den Bereich [0, 255]
    if np.ptp(image_data) > 0:  # Prüfen, ob es Pixelwerte mit Abweichungen gibt
        image_data = (255 * (image_data - np.min(image_data)) / np.ptp(image_data)).astype(np.uint8)
    else:
        image_data = np.zeros_like(image_data, dtype=np.uint8)  # Schwarzes Bild für konstante Werte

    return image_data

def detect_mode_from_header(fits_header):
    """Erkennt den Modus basierend auf dem FITS-Header."""
    if fits_header.get('COLRCCD', 0) == 0 and fits_header.get('NAXIS') == 2:
        return "mono"
    elif fits_header.get('COLRCCD', 1) == 1 and fits_header.get('NAXIS') == 3 and fits_header.get('NAXIS3') == 3:
        return "rgb"
    elif fits_header.get('NAXIS') == 3 and fits_header.get('NAXIS3') == 4:
        return "lrgb"
    elif fits_header.get('NAXIS') == 2:
        return "rggb"
    elif fits_header.get('NAXIS') == 3 and fits_header.get('NAXIS3') == 3:
        return "yuyv422"
    else:
        raise ValueError("Das Datenformat konnte nicht erkannt werden. Bitte Modus manuell angeben.")

def fits_to_tiff(fits_path, output_path, mode, summary):
    """Konvertiert eine FITS-Datei in eine TIFF-Datei."""
    with fits.open(fits_path) as hdul:
        image_data = hdul[0].data
        header = hdul[0].header

    if mode == "auto":
        detected_mode = detect_mode_from_header(header)
        mode_mapping = {
            "mono": "mono8",
            "rgb": "rgb",
            "lrgb": "lrgb",
            "rggb": "rggb",
            "yuyv422": "yuyv422",
        }
        mode = mode_mapping.get(detected_mode, None)
        if not mode:
            raise ValueError(f"DEBUG: Automatisch erkannter Modus '{detected_mode}' ist nicht unterstützt.")

    summary[mode] += 1

    if mode == "mono8":
        if image_data.ndim > 2:
            image_data = image_data[0]
        image_data = preprocess_image_data(image_data, header)
        image = Image.fromarray(image_data, mode="L")

    elif mode == "mono16":
        if image_data.ndim > 2:
            image_data = image_data[0]
        image_data = (65535 * (image_data - np.min(image_data)) / np.ptp(image_data)).astype(np.uint16)
        image = Image.fromarray(image_data, mode="I;16")

    elif mode == "rgb":
        if image_data.ndim == 3 and image_data.shape[0] == 3:
            r = preprocess_image_data(image_data[0], header)
            g = preprocess_image_data(image_data[1], header)
            b = preprocess_image_data(image_data[2], header)
            rgb_image = np.stack((r, g, b), axis=-1)
            image = Image.fromarray(rgb_image, mode="RGB")
        else:
            raise ValueError("Das Bildformat passt nicht zu RGB (erwartet 3 Kanäle).")

    elif mode == "lrgb":
        if image_data.ndim == 3 and image_data.shape[0] == 4:
            l = preprocess_image_data(image_data[0], header)
            r = preprocess_image_data(image_data[1], header)
            g = preprocess_image_data(image_data[2], header)
            b = preprocess_image_data(image_data[3], header)
            rgb_image = np.stack((r, g, b), axis=-1)
            rgb_image = (rgb_image * l[:, :, None]).astype(np.uint8)
            image = Image.fromarray(rgb_image, mode="RGB")
        else:
            raise ValueError("Das Bildformat passt nicht zu LRGB (erwartet 4 Kanäle).")

    elif mode == "rggb":
        if image_data.ndim == 2:
            h, w = image_data.shape
            r = image_data[0:h:2, 0:w:2]
            g1 = image_data[0:h:2, 1:w:2]
            g2 = image_data[1:h:2, 0:w:2]
            b = image_data[1:h:2, 1:w:2]
            g = (g1 + g2) // 2
            rgb_image = np.stack((r, g, b), axis=-1)
            rgb_image = preprocess_image_data(rgb_image, header)
            image = Image.fromarray(rgb_image, mode="RGB")
        else:
            raise ValueError("Das Bildformat passt nicht zu RGGB (erwartet 2D-Bild).")

    elif mode == "yuyv422":
        if image_data.ndim == 3 and image_data.shape[0] == 3:
            rgb_image = yuyv_to_rgb(image_data)
            image = Image.fromarray(rgb_image, mode="RGB")
        else:
            raise ValueError("Das Bildformat passt nicht zu YUYV 4:2:2.")

    else:
        raise ValueError(f"Ungültiger Modus '{mode}'. Unterstützte Modi sind: mono8, mono16, rgb, lrgb, rggb, yuyv422.")

    image.save(output_path, format="TIFF")

def print_progress(current, total):
    """Druckt einen Fortschrittsbalken in einer Zeile."""
    progress = int((current / total) * 30)
    bar = "[" + "#" * progress + "-" * (30 - progress) + "]"
    print(f"\r{bar} {current}/{total}", end="", flush=True)

def main():
    parser = argparse.ArgumentParser(
        description="Konvertiert FITS-Dateien in TIFF-Dateien.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("-i", "--input", metavar="File", type=str, help="Pfad zur FITS-Eingabedatei")
    parser.add_argument("-o", "--output", metavar="File", type=str, help="Pfad zur TIFF-Ausgabedatei")
    parser.add_argument("-d", "--directory", metavar="Path", type=str, help="Pfad zum Verzeichnis mit FITS-Dateien")
    parser.add_argument(
        "-m", "--mode",
        choices=["mono8", "mono16", "rgb", "rggb", "lrgb", "yuyv422", "auto"],
        required=True,
        help=(
            "  mono8    (8-Bit-Graustufen)\n"
            "  mono16   (16-Bit-Graustufen)\n"
            "  rgb      (RGB)\n"
            "  rggb     (Bayer-Pattern)\n"
            "  lrgb     (Luminance + RGB)\n"
            "  yuyv422  (YUYV 4:2:2)\n"
            "  auto     (automatische Erkennung)"
        )
    )

    args = parser.parse_args()

    if (args.input or args.output) and args.directory:
        parser.print_help()
        print("\nFehler: Die Optionen -i und -o können nicht zusammen mit -d verwendet werden.")
        exit(1)
    if args.input and not args.output:
        parser.print_help()
        print("\nFehler: -o (Pfad zur TIFF-Ausgabedatei) ist erforderlich, wenn -i (Eingabedatei) angegeben wird.")
        exit(1)
    if not args.input and args.output:
        parser.print_help()
        print("\nFehler: -i (Pfad zur FITS-Eingabedatei) ist erforderlich, wenn -o (Ausgabedatei) angegeben wird.")
        exit(1)
    if not (args.input or args.directory):
        parser.print_help()
        print("\nFehler: Bitte entweder -i und -o für Einzeldateien oder -d für ein Verzeichnis angeben.")
        exit(1)

    summary = {"mono8": 0, "mono16": 0, "rgb": 0, "rggb": 0, "lrgb": 0, "yuyv422": 0}

    if args.input:
        start_time = time.time()
        fits_to_tiff(args.input, args.output, args.mode, summary)
        elapsed_time = time.time() - start_time
        print("\nKonvertierung abgeschlossen.")
        print(f"Benötigte Zeit: {elapsed_time:.2f} Sekunden")
        print("Zusammenfassung:")
        for mode, count in summary.items():
            if count > 0:
                print(f"  {mode}: {count} Dateien")
    elif args.directory:
        if not os.path.isdir(args.directory):
            print("Das angegebene Verzeichnis existiert nicht.")
            exit(1)

        fits_files = [
            os.path.join(args.directory, f) for f in os.listdir(args.directory) if f.lower().endswith(".fits")
        ]

        if not fits_files:
            print("Keine FITS-Dateien im angegebenen Verzeichnis gefunden.")
            exit(1)

        print(f"Gefundene FITS-Dateien: {len(fits_files)}")

        start_time = time.time()

        for i, fits_path in enumerate(fits_files, 1):
            output_path = fits_path.rsplit(".", 1)[0] + ".tif"
            fits_to_tiff(fits_path, output_path, args.mode, summary)
            print_progress(i, len(fits_files))

        elapsed_time = time.time() - start_time
        print("\nKonvertierung abgeschlossen.")
        print(f"Benötigte Zeit: {elapsed_time:.2f} Sekunden")
        print("Zusammenfassung:")
        for mode, count in summary.items():
            if count > 0:
                print(f"  {mode}: {count} Dateien")


if __name__ == "__main__":
    main()

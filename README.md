# dna_stadistics
This project is an application in Python to analyze and verify DNA and RNA statistics. Its part of the projects developed during the biomedical engineering bachelor course at the University of Alicante between 2021 and 2025.
Subject: Fundamentals of Programming


## What does the project do?

- **Generate random DNA sequences:**  
  Using `generador_cadenas`, the program creates a random DNA sequence made up of the bases A, T, C, and G.

- **Clean sequences:**  
  The `limpiarCadena` function removes invalid characters, leaving only valid bases.

- **Calculate base percentages:**  
  `contarBase` calculates the percentage of a specific base (or group of bases) within a sequence.

- **Transcribe DNA to RNA:**  
  The `transcripcion` and `transcripcion_ADN` functions generate messenger RNA (mRNA) from a DNA sequence following base pairing rules.

- **Translate RNA to protein:**  
  `traduccion` and `traduccion_ADN` translate an RNA sequence into its corresponding amino acid sequence using the standard genetic code.

- **Convert RNA to DNA:**  
  `ARN_ADN` converts an RNA strand back into a DNA sequence.

- **Generate the complementary DNA strand:**  
  `complementaria` returns the complementary strand of the input DNA sequence.

- **Display statistics and results:**  
  With `estadisticas_CADENA` and `estadisticas_PROTEINA`, users can view detailed information such as sequence type (DNA/RNA), base percentages, transcription, translation, and more.

## Technologies used

- **Python** – Main programming language  
- **Tkinter** – For the graphical user interface  
- **Basic bioinformatics algorithms** – Custom implemented for genetic analysis

## How to run it

1. Make sure Python is installed (Python 3.8 or higher recommended).
2. Clone this repository or download the source code.
3. Run the main script:

   ```bash
   python interfaz_grafica.py


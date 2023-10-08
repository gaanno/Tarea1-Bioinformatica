#archivo encargado de usar clustal omega

import subprocess
import os

class Clustal(object):
    _OUTPUT_FILENAME = "clustal_result"
    _CLUSTAL_DIR = os.path.join(os.getcwd(), 'clustal', 'clustalo.exe')

    def __init__(self, input_file, file_format):
        self._input_file = os.path.join(os.getcwd(), input_file)
        self._arguments = f"{self._CLUSTAL_DIR} --in {input_file} -t Protein --infmt fa --outfmt={file_format} -o {self._OUTPUT_FILENAME}.{file_format} --force --threads 12"

    def run(self):
        """
            Ejecuta el comando de clustal omega
        """
        if not os.path.isfile(self._input_file):
            print("Archivo no encontrado")
            return

        subprocess.run(self._arguments, shell=True, check=True)
    

    
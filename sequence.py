import requests


class Sequence(object):
    _sequences = {}
    _OUTPUT_FILENAME = "all_sequences"

    def __init__(self, FILE_FORMAT) -> None:
        """Constructor"""
        self._FILE_FORMAT = FILE_FORMAT

    def add_sequence(self, database, accession_number) -> None:
        """
            Agrega una nueva secuencia a la lista de secuencias

            @param database: base de datos de la secuencia
            @param accession_number: numero de acceso de la secuencia

            @return: None
        """
        self._sequences[accession_number] = self._get_data(
            database=database, accession_number=accession_number)

    def _get_data(self, database, accession_number) -> str:
        """
            Obtiene la secuencia desde GenBank

            @param database: base de datos de la secuencia
            @param accession_number: numero de acceso de la secuencia

            @return: str con la secuencia o None si no se encuentra
        """
        URL = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={database}&id={accession_number}&rettype={self._FILE_FORMAT}"
        response = requests.get(URL)

        return response.text if response.status_code == 200 else None

    def get_sequence(self, accession_number) -> str:
        """
            Obtiene una secuencia ya agregada

            @param accession_number: numero de acceso de la secuenci|a

            @return: str con la secuencia o None si no se encuentra
        """
        return self._sequences.get(accession_number, None)

    def print_sequence(self, accession_number) -> None:
        """Imprime una secuencia ya agregada"""
        print(self._sequences.get(accession_number, None))

    def get_all_sequences(self) -> str:
        """Obtiene todas las secuencias agregadas"""
        values = [str(v) for v in self._sequences.values()]
        return "".join(values)

    def export(self, accession_number) -> None:
        """Exporta una secuencia determinada al FILE_FORMAToe stablecido"""
        if self._sequences.get(accession_number, None) is None:
            print("Secuencia no encontrada")
            return None

        with open(f'{accession_number}.{self._FILE_FORMAT}', 'w') as f:
            f.write(str(self._sequences[accession_number]))

    def export_all(self) -> None:
        """Exporta todas las secuencias"""
        with open(f'{self._OUTPUT_FILENAME}.{self._FILE_FORMAT}', 'w') as f:
            f.write(str(self.get_all_sequences()))

    def get_output_filename(self) -> str:
        """Obtiene el nombre del archivo de salida"""
        return f'{self._OUTPUT_FILENAME}.{self._FILE_FORMAT}'
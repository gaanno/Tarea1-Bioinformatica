import requests

class Sequence(object):
    _sequences = {}
    _format = "fasta"

    def add_sequence(self, database, accesion_number) -> None:
        """
            Agrega una nueva secuencia a la lista de secuencias

            @param database: base de datos de la secuencia
            @param accesion_number: numero de acceso de la secuencia

            @return: None
        """
        self._sequences[accesion_number] = self._get_data(
            database=database, accesion_number=accesion_number)

    def _get_data(self, database, accesion_number) -> str:
        """
            Obtiene la secuencia desde GenBank

            @param database: base de datos de la secuencia
            @param accesion_number: numero de acceso de la secuencia

            @return: str con la secuencia o None si no se encuentra
        """
        URL = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db={database}&id={accesion_number}&rettype={self._format}"
        response = requests.get(URL)

        return response.text if response.status_code == 200 else None

    def get_sequence(self, accesion_number) -> str:
        """
            Obtiene una secuencia ya agregada

            @param accesion_number: numero de acceso de la secuenci|a

            @return: str con la secuencia o None si no se encuentra
        """
        return self._sequences.get(accesion_number, None)

    def print_sequence(self, accesion_number) -> None:
        """Imprime una secuencia ya agregada"""
        print(self._sequences.get(accesion_number, None))

    def get_all_sequences(self) -> dict:
        """Obtiene todas las secuencias agregadas"""
        values = [str(v) for v in self._sequences.values()]
        return "".join(values)

    def export(self, accesion_number) -> None:
        """Exporta una secuencia determinada al formatoe stablecido"""
        if self._sequences.get(accesion_number, None) is None:
            print("Secuencia no encontrada")
            return None

        with open(f'{accesion_number}.{self._format}', 'w') as f:
            f.write(str(self._sequences[accesion_number]))

    def export_all(self) -> None:
        """Exporta todas las secuencias"""
        for key in self._sequences.keys():
            self.export(key)
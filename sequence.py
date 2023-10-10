from Bio import Entrez
from Bio import SeqIO


class Sequence(object):
    _OUTPUT_FILENAME = "all_sequences"
    _sequence_pairs = []
    database: str
    file_format: str

    def __init__(self, email="") -> None:
        """Constructor"""
        Entrez.email = email

    def add_sequence(self, seqs: list) -> None:
        """
            Agrega una nueva secuencia a la lista de secuencias

            @param database: base de datos de la secuencia
            @param accession_number: numero de acceso de la secuencia

            @return: None
        """
        data = self._get_data(accession_number=seqs)
        self._sequence_pairs.append(data)

    def _get_data(self, accession_number) -> str:
        """
            Obtiene la secuencia desde GenBank

            @param database: base de datos de la secuencia
            @param accession_number: numero de acceso de la secuencia

            @return: list con la secuencia o None si no se encuentra
        """
        with Entrez.efetch(db=self.database, id=accession_number,
                           rettype=self.file_format) as handle:
            return [handle for handle in SeqIO.parse(handle, self.file_format)]

    def export(self) -> None:
        """Exporta todas las secuencias"""
        for pair in self._sequence_pairs:
            for seq in pair:
                with open(f"{seq.id}.{self.file_format}", "w") as output_handle:
                    SeqIO.write(seq, output_handle, self.file_format)

    def print_sequences(self) -> None:
        """Imprime todas las secuencias agregadas"""
        for pair in self._sequence_pairs:
            for sequence in pair:
                print(f"{sequence}")
                # print(f"{sequence.id} {sequence.seq}")

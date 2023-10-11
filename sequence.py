from Bio import Entrez
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

class Sequence(object):
    _OUTPUT_FILENAME = "all_sequences"
    _sequence_pairs = []
    _global_alignments = []
    _local_alignments = []

    def __init__(self, database, file_format="gb", email="") -> None:
        """Constructor"""
        Entrez.email = email
        self.database = database
        self.file_format = file_format

    def add_sequence(self, seqs: list) -> None:
        """
            Agrega una nueva secuencia a la lista de secuencias

            @param database: base de datos de la secuencia
            @param accession_number: numero de acceso de la secuencia

            @return: None
        """
        print(f"Getting sequences {seqs}")
        data = self._get_data(accession_number=seqs)
        print("Adding sequence")
        self._sequence_pairs.append(data)
        
    def export(self) -> None:
        """Exporta todas las secuencias"""
        for pair in self._sequence_pairs:
            for seq in pair:
                with open(f"{seq.id}.{self.file_format}", "w") as output_handle:
                    SeqIO.write(seq, output_handle, self.file_format)

    def do_alignment(self, type):
        aligner = PairwiseAligner()
        aligner.mode = type

        for pair in self._sequence_pairs:
            alignment = aligner.align(*pair)
            ids = [seq.id for seq in pair]
            if type == 'global':
                alignment = pairwise2.align.globalxx(pair[0].seq, pair[1].seq)
                self._global_alignments.append([ids, alignment])

            elif type == 'local':
                alignment = pairwise2.align.localxx(pair[0].seq, pair[1].seq)
                self._local_alignments.append([ids, alignment])

    def print_sequences(self) -> None:
        """Imprime todas las secuencias agregadas"""
        for pair in self._sequence_pairs:
            for sequence in pair:
                print(sequence)
                #print(f"{sequence.id} {sequence.seq}")

    def print_alignments(self, type) -> None:
        """Imprime todos los alineamientos globales"""
        print(f"--- {type} ---")
        if type == 'global':
            for alignment in self._global_alignments:
                print(f"Secuencias: {alignment[0]}\n{format_alignment(*alignment[1][0])}")

        if type == 'local':
            for alignment in self._local_alignments:
                print(f"Secuencias: {alignment[0]}\n{format_alignment(*alignment[1][0])}")

    def show_data(self):
        # muestra score, porcentaje de identidad y de similitud
        print("--- global ---")
        for gl in self._global_alignments:
            identity = (gl[1][0].score/gl[1][0].end)*100
            print(
                f"Elementos: {gl[0]}\nscore: {gl[1][0].score}\nIdentidad: {identity}%", end="\n\n")

        print("--- local ---")
        for lc in self._local_alignments:
            identity = (lc[1][0].score/lc[1][0].end)*100
            print(
                f"Elementos: {lc[0]}\nscore: {lc[1][0].score}\nIdentidad: {identity}%", end="\n\n")

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

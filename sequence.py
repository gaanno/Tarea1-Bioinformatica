from Bio import Entrez
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices


class Sequence(object):
    _OUTPUT_FILENAME = "all_sequences"
    _sequence_pairs = []

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
        self._sequence_pairs.append(data)

    def export(self) -> None:
        """Exporta todas las secuencias"""
        for pair in self._sequence_pairs:
            for seq in pair:
                with open(f"{seq.id}.{self.file_format}", "w") as output_handle:
                    SeqIO.write(seq, output_handle, self.file_format)

    def do_alignment(self):

        for pair in self._sequence_pairs:
            alignment = pairwise2.align.globalxx(pair[0].seq, pair[1].seq)
            score, identity = self._get_resume(alignment)
            print(
                f"Aligning {pair[0].id, pair[1].id} type:global {score=} {identity=}")

            alignment = pairwise2.align.localxx(pair[0].seq, pair[1].seq)
            score, identity = self._get_resume(alignment)
            print(
                f"Aligning {pair[0].id, pair[1].id} type:local {score=} {identity=}\n")

    def print_sequences(self) -> None:
        """Imprime todas las secuencias agregadas"""
        for pair in self._sequence_pairs:
            for sequence in pair:
                print(sequence)
                # print(f"{sequence.id} {sequence.seq}")

   
    def matrices(self):
        for pair in self._sequence_pairs:
            print(f"Aligning {pair[0].id, pair[1].id} for BLOSUM62 and PAM250")
            alignment = pairwise2.align.globaldx(
                pair[0].seq, pair[1].seq, substitution_matrices.load("BLOSUM62"))
            score, identity = self._get_resume(alignment)
            print(
                f"Blosum62 Aligning {pair[0].id, pair[1].id} score {score} identity {identity}")

            alignment = pairwise2.align.globaldx(
                pair[0].seq, pair[1].seq, substitution_matrices.load("PAM250"))
            score, identity = self._get_resume(alignment)
            print(
                f"PAM250 Aligning {pair[0].id, pair[1].id} score {score} identity {identity}\n")

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

    def _get_resume(self,alignment) -> tuple:
        score = str(format_alignment(*alignment[0])).count("|")
        identity = (score/alignment[0].end)*100
        return score, identity

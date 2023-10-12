from Bio import Entrez
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Align import substitution_matrices


class Sequence(object):
    _OUTPUT_FILENAME = "all_sequences"
    _sequence_pairs = []
    gap = -10
    match = 2
    mismatch = -4

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
            # alignment = pairwise2.align.globalxx(pair[0].seq, pair[1].seq)
            self._global_alignment(pair[0], pair[1])
            self._local_alignment(pair[0], pair[1])
            print()

    def print_sequences(self) -> None:
        """Imprime todas las secuencias agregadas"""
        for pair in self._sequence_pairs:
            for sequence in pair:
                print(sequence)

    def do_alignment_matrices(self):
        for pair in self._sequence_pairs:
            self._alignment_blosum62(pair[0], pair[1])
            self._alignment_pam250(pair[0], pair[1])
            self._global_alignment(pair[0], pair[1])
            self._local_alignment(pair[0], pair[1])
            print()

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

    def _alignment_pam250(self, seq1, seq2):
        alignment = pairwise2.align.globaldx(
            seq1.seq, seq2.seq, substitution_matrices.load("PAM250"))
        self._get_resume(alignment, seq1.id, seq2.id, "PAM250")

    def _alignment_blosum62(self, seq1, seq2):
        alignment = pairwise2.align.globaldx(
            seq1.seq, seq2.seq, substitution_matrices.load("BLOSUM62"))
        self._get_resume(alignment, seq1.id, seq2.id, "blosum62")

    def _global_alignment(self, seq1, seq2):
        alignment = pairwise2.align.globalms(
            seq1, seq2, self.match, self.mismatch, self.gap, self.gap)
        self._get_resume(alignment, seq1.id, seq2.id, "global")

    def _local_alignment(self, seq1, seq2):
        alignment = pairwise2.align.localms(
            seq1, seq2, self.match, self.mismatch, self.gap, self.gap)
        self._get_resume(alignment, seq1.id, seq2.id, "local")

    def _get_resume(self, alignment, seq1, seq2, type) -> tuple:
        score = str(format_alignment(*alignment[0])).count("|")
        identity = (score/alignment[0].end)*100
        print(f"Sequences: {seq1, seq2} {type=} {score=} {identity=}%")

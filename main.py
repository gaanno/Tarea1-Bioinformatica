# .\clustalo.exe --in prueba -t Protein --infmt fa --outfmt=fasta -o resultado.txt --force --threads 12
from sequence import Sequence
from exec_clustal import Clustal


if __name__ == '__main__':
    
    
    file_format = "fasta"

    seq = Sequence(file_format)
    clustal = Clustal(seq.get_output_filename(), file_format)

    seq.add_sequence(database="protein", accession_number="BAA20512.1")
    seq.add_sequence(database="protein", accession_number="CAA23748.1")
    seq.add_sequence(database="protein", accession_number="CAA24095.1")

    seq.export_all()
    clustal.run()
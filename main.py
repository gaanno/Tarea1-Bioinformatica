from sequence import Sequence
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-file', type=str, help='File containing accession numbers')

    if not (file := parser.parse_args().file) or file is None:
        parser.print_help()
        exit()
        
    return parser.parse_args()

def load_accession_numbers(file, seq: Sequence):
    with open(file) as f:
        for line in f.readlines():
            seq1, seq2 = line.split(' ')
            seq.add_sequence([seq1, seq2])

def main():
    parser = parse_args()
    seq = Sequence()
    seq.database = "protein"
    seq.file_format = "gb"

    load_accession_numbers(parser.file, seq)

    seq.print_sequences()
    seq.export()

if __name__ == '__main__':
    main()
from sequence import Sequence
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-file', type=str, help='File container of the accession numbers')

    if not (file := parser.parse_args().file) or file is None:
        parser.print_help()
        exit()
        
    return parser.parse_args()

def load_accession_numbers(file, seq):
    with open(file) as f:
        #return [line.strip() for line in f.readlines()]
        for line in f.readlines():
            acc, type = line.split(' ')
            seq.add_sequence(database=type, accession_number=acc)

def main():
    file_format = "fasta"
    parser = parse_args()
    seq = Sequence(file_format)

    load_accession_numbers(parser.file, seq)

    print(seq.get_all_sequences())


if __name__ == '__main__':
    main()
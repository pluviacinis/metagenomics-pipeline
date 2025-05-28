#!/usr/bin/env python3
import pyhmmer
from pyhmmer import easel, plan7
import argparse

def safe_decode(value, default="-"):
    return value.decode() if value is not None else default

def main():
    parser = argparse.ArgumentParser(description='PyHMMER scan для Phamb с доменами')
    parser.add_argument('hmm_db', help='Путь к HMM-базе')
    parser.add_argument('seq_file', help='FASTA-файл с последовательностями')
    parser.add_argument('-o', '--output', required=True, help='Выходной tbl-файл')
    parser.add_argument('--cpu', type=int, default=4, help='Число потоков')
    parser.add_argument('-E', type=float, default=10.0, help='Порог E-value')
    args = parser.parse_args()

    with plan7.HMMFile(args.hmm_db) as hmm_file:
        hmms = list(hmm_file)
        Z = len(hmms)
    with easel.SequenceFile(args.seq_file, digital=True) as seq_file:
        sequences = seq_file.read_block()

    results = pyhmmer.hmmscan(sequences, hmms, cpus=args.cpu, E=args.E, Z=Z)

    with open(args.output, 'w') as tblout:
        # Заголовок Phamb
        tblout.write("# target name\taccession\tquery name\taccession\tE-value\tScore\tBias\tDom_E-value\tDom_Score\tDom_Bias\n")
        
        for result in results:
            query_name = safe_decode(result.query.name)
            for hit in result:
                if hit.evalue > args.E:
                    continue
                target_acc = safe_decode(hit.accession)
                for domain in hit.domains:
                    line = (
                        f"{safe_decode(hit.name):<45}\t{target_acc}\t"
                        f"{query_name}\t{target_acc}\t"
                        f"{hit.evalue:.1g}\t{hit.score:.1f}\t{hit.bias:.1f}\t"
                        f"{domain.i_evalue:.1g}\t{domain.score:.1f}\t{domain.bias:.1f}\n"
                    )
                    tblout.write(line)

if __name__ == "__main__":
    main()

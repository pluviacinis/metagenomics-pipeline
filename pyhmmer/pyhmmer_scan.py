#!/usr/bin/env python3
import pyhmmer
from pyhmmer import easel, plan7
import argparse

def safe_decode(value, default="N/A") -> str:
    """Безопасное декодирование bytes в str с обработкой None"""
    return value.decode() if value is not None else default

def main():
    parser = argparse.ArgumentParser(
        description='PyHMMER scan: HMMER-compatible hmmscan implementation'
    )
    parser.add_argument('hmm_db', help='Path to HMM database (.hmm or .h3m)')
    parser.add_argument('seq_file', help='Input sequences (FASTA format)')
    parser.add_argument('-o', '--output', required=True, help='Output tbl file')
    parser.add_argument('--cpu', type=int, default=4, help='Number of threads')
    parser.add_argument('-E', type=float, default=10.0, 
                       help='E-value threshold (default: 10.0)')
    args = parser.parse_args()

    # 1. Загрузка HMM-профилей
    with plan7.HMMFile(args.hmm_db) as hmm_file:
        hmms = list(hmm_file)
        Z = len(hmms)

    # 2. Загрузка последовательностей
    with easel.SequenceFile(args.seq_file, digital=True) as seq_file:
        sequences = seq_file.read_block()

    # 3. Запуск hmmscan
    results = pyhmmer.hmmscan(
        sequences,
        hmms,
        cpus=args.cpu,
        E=args.E,
        Z=Z
    )

    # 4. Форматирование вывода с обработкой None
    with open(args.output, 'w') as tblout:
        tblout.write("# Target Sequence          Accession  E-value  Score  Bias\n")
        tblout.write("#------------------        ---------  -------  -----  ----\n")
        
        for result in results:
            for hit in result:
                if hit.evalue <= args.E:
                    name = safe_decode(hit.name, "<no_name>")
                    accession = safe_decode(hit.accession, "<no_accession>")
                    line = f"{name:<25} {accession:<10} " \
                           f"{hit.evalue:.1g} {hit.score:.1f} {hit.bias:.1f}\n"
                    tblout.write(line)

if __name__ == "__main__":
    main()

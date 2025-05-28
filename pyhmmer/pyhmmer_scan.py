#!/usr/bin/env python3
import pyhmmer
from pyhmmer import easel, plan7
import argparse

def safe_decode(value, default="-"):
    """Безопасное декодирование bytes в str с обработкой None"""
    return value.decode() if value is not None else default

def main():
    parser = argparse.ArgumentParser(
        description='PyHMMER scan: Phamb-compatible hmmscan'
    )
    parser.add_argument('hmm_db', help='Path to HMM database')
    parser.add_argument('seq_file', help='Input sequences (FASTA)')
    parser.add_argument('-o', '--output', required=True, help='Output tbl file')
    parser.add_argument('--cpu', type=int, default=4, help='Threads')
    parser.add_argument('-E', type=float, default=10.0, help='E-value threshold')
    args = parser.parse_args()

    # Загрузка данных
    with plan7.HMMFile(args.hmm_db) as hmm_file:
        hmms = list(hmm_file)
        Z = len(hmms)
    with easel.SequenceFile(args.seq_file, digital=True) as seq_file:
        sequences = seq_file.read_block()

    # Запуск hmmscan
    results = pyhmmer.hmmscan(sequences, hmms, cpus=args.cpu, E=args.E, Z=Z)

    # Форматирование вывода для Phamb
    with open(args.output, 'w') as tblout:
        # Заголовки Phamb
        header = """#                                                                                                   --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name                                            accession  query name           accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#                                    ------------------- ---------- -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------\n"""
        tblout.write(header)
        
        for result in results:
            # Получаем имя HMM-профиля (query)
            query_name = safe_decode(result.query.name)
            
            # Итерируемся по хитам в TopHits
            for hit in result:
                if hit.evalue > args.E:
                    continue
                
                # Основные поля
                target_name = safe_decode(hit.name)
                target_acc = safe_decode(hit.accession)
                full_eval = f"{hit.evalue:.1g}"
                full_score = f"{hit.score:.1f}"
                full_bias = f"{hit.bias:.1f}"
                
                # Поля доменов (заглушки)
                line = (
                    f"{target_name:<45} {target_acc:<10} "
                    f"{query_name:<20} {target_acc:<10} "
                    f"{full_eval:>9} {full_score:>6} {full_bias:>5}   "
                    "NaN      NaN    NaN   "  # Заглушка для доменных данных
                    "1.0   1   1   0   1   1   1   1 #\n"
                )
                tblout.write(line)

if __name__ == "__main__":
    main()

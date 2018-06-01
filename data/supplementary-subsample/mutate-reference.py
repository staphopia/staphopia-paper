#! /usr/bin/env python3
import argparse as ap
import os
import random
import sys


def split_seq(seq, length):
    return (seq[0 + i:length + i] for i in range(0, len(seq), length))


def create_insertion(reference_base, length):
    bases = [reference_base]
    for i in range(length):
        bases.append(random.choice('ATGC'))

    return ''.join(bases)


def create_snp(base):
    alternates = 'ATGC'.replace(base, '')
    return random.sample(list(alternates), 1)[0]


if __name__ == '__main__':
    parser = ap.ArgumentParser(prog='mutate-reference.py',
                               conflict_handler='resolve',
                               description="Mutate a given reference FASTA")

    parser.add_argument('reference', type=str, metavar="REFERENCE_FASTA",
                        help=('Input FASTA sequence to mutate.'))
    parser.add_argument('output', type=str, metavar="OUTPUT_DIRECTORY",
                        help='Directory to place output files.')
    parser.add_argument('--num_variants', default=100, type=int, metavar="INT",
                        help=('Number of variants to output. (Default 100)'))
    parser.add_argument('--snp', default=0.95, type=float, metavar="FLOAT",
                        help=('The percent of mutations that are SNPs. '
                              '(Default: 0.95, InDels=1-SNPs)'))
    parser.add_argument('--seed', default=0, type=int, metavar="INT",
                        help=('A seed to set for reproducible results.'))

    args = parser.parse_args()

    if args.seed:
        random.seed(args.seed)

    # Read FASTA file
    seqs = []
    with open(args.reference, "r") as f:
        for line in f:
            line = line.rstrip()
            if not line.startswith('>'):
                seqs.append(line.upper())

    seq = list(''.join(seqs))
    # Reduce the total number of possible variant sites incase there are excess
    # simualted deletions.
    seq_length = len(seq) - (args.num_variants * 5)

    # Create variant positions
    variants = sorted(random.sample(range(seq_length), args.num_variants))
    total = 0

    reference = os.path.basename(os.path.splitext(args.reference)[0])
    variant_file = '{0}/{1}-{2}-variants.txt'.format(
        args.output, reference, args.num_variants
    )
    fasta_file = '{0}/{1}-{2}.fasta'.format(
        args.output, reference, args.num_variants
    )

    args_file = '{0}/{1}-{2}-args.txt'.format(
        args.output, reference, args.num_variants
    )

    with open(variant_file, 'w') as f:
        for variant in variants:
            if seq[variant] is not 'X':
                # zero-based index, position = position + 1
                position = variant + 1

                reference_base = seq[variant]
                alternate_base = None

                if random.random() >= args.snp:
                    # InDel
                    length = random.randint(3, 15) - 1
                    if random.randint(0, 1):
                        # Insertion
                        alternate_base = create_insertion(reference_base,
                                                          length)
                        seq[variant] = alternate_base
                    else:
                        # Deletion
                        alternate_base = reference_base
                        reference_base = ''.join(seq[variant:variant + length])
                        total += len(reference_base) - len(alternate_base) - 1
                        for i in range(variant + 1, variant + length):
                            seq[i] = 'X'
                else:
                    # SNP
                    alternate_base = create_snp(reference_base)
                    seq[variant] = alternate_base

                f.write('{0}\t{1}\t{2}\n'.format(
                    position, reference_base, alternate_base
                ))

    with open(fasta_file, 'w') as f:
        f.write('>{0}:{1} variants\n'.format(reference, args.num_variants))
        seq = ''.join(seq).replace('X', '')
        for seq in split_seq(seq, 79):
            f.write('{0}\n'.format(seq))

    with open(args_file, 'w') as f:
        f.write('{0}\n'.format(' '.join(sys.argv)))

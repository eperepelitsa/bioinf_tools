import argparse
from pathlib import Path
from main import filter_fastq
from logger import setup_logger

def create_parser():
    parser = argparse.ArgumentParser(
        description='Filter FASTQs by GC content, length, and quality'
    )
    parser.add_argument('-i', '--input', required=True,
                        help='Input FASTQ file path')
    parser.add_argument('-o', '--output', required=True,
                        help='Output FASTQ file path')
    parser.add_argument('-g', '--gc-bounds', nargs=2, type=float,
                        default=[0, 100], metavar=('MIN', 'MAX'),
                        help='GC content bounds (default: 0-100)')
    parser.add_argument('-l', '--length-bounds', nargs=2, type=int,
                        default=[0, 2**32], metavar=('MIN', 'MAX'),
                        help='Length bounds (default: 0-2**32)')
    parser.add_argument('-q', '--quality-threshold', type=float,
                        default=0, help='Minimum average quality score')
    parser.add_argument('--log', default='filter_fastq.log',
                        help='Log file path (default: filter_fastq.log)')
    return parser

if __name__ == '__main__':
    parser = create_parser()
    args = parser.parse_args()
    logger = setup_logger(args.log)
    logger.info(f"Starting filtering: {args.input} -> {args.output}")
    try:
        stats = filter_fastq(
            input_file=args.input,
            output_file=args.output,
            gc_bounds=tuple(args.gc_bounds),
            length_bounds=tuple(args.length_bounds),
            quality_threshold=args.quality_threshold
        )
        logger.info(f"Filtering complete. {stats['passed']}/{stats['total']} sequences passed")
    except Exception as e:
        logger.error(f'Error: {e}')
        raise
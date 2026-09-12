"""Preparation-only CLI: no fit or training subcommand."""
import argparse


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    diagnostic=sub.add_parser('diagnose',help='deterministic correctness checks, no sampled data')
    diagnostic.add_argument('--out',required=True)
    a=p.parse_args()
    from .diagnostics import run
    run(a.out)


if __name__=='__main__':
    main()

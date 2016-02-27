import json
import logging
from logging.config import fileConfig

from genome import Genome

fileConfig('logging.conf')


def parse_input(filename):
    """Parses the input file and returns a parameter dictionary."""
    with open(filename) as f:
        data = json.load(f)
    return data


def main():
    filename = '../tests/input.json'
    logging.info("Reading input from %s" % filename)
    data = parse_input(filename)
    logging.info("Started: create genomes")
    genomes = [Genome(g['name'], g['accession_numbers'])
               for g in data['genomes']]
    logging.info("Finished: create genomes")

import sys
import os
import logging
import multiprocessing
import subprocess
import glob
import shutil
import shlex
import click

from snakemake import load_configfile
from ruamel.yaml import YAML
from mosaic import __version__
from mosaic.config import get_default_config,set_logger

set_logger()
logging.info(f'mosaic {__version__}')
logging.info(' '.join(map(shlex.quote, sys.argv)))

def log_exception(msg):
    logging.critical(msg)
    logging.info('Documentation is available at: '
                    'https://github.com/Mosaic/wiki')
    logging.info('Issues can be raised at: '
                    'https://github.com/Mosaic/issues')
    sys.exit(1)

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """
    mosaic - modular viral metagenomic analysis
    """

def get_snakefile(f="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), f)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

@cli.command(
    'run',
    context_settings=dict(ignore_unknown_options=True),
    short_help='run vcontact'
)
@click.option('-w',
    '--input-dir',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    help='output directory',
    default='.'
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='max # of jobs allowed in parallel.',
)
@click.option(
    '-n',
    '--dryrun',
    is_flag=True,
    default=False,
    show_default=True,
    help='Check rules to run and files to produce',
)
@click.option(
    '-r',
    '--reason',
    is_flag=True,
    default=False,
    show_default=True,
    help='Print the reason for each executed rule.',
)
def runWorkflow(input_dir, jobs, dryrun, reason):
    '''Training customized classifier model.
    '''

    DEFAULT_CONFIG = get_default_config()

    cmd = (
        'snakemake '
        '--config '
            'input_dir={input_dir} '
        '--jobs {jobs} --rerun-incomplete --latency-wait 600 {dryrun} {reason}'
    ).format(
        input_dir=input_dir,
        jobs=jobs,
        dryrun='--dryrun' if dryrun else '',
        reason='--reason' if reason else '',
    )
    logging.info('Executing: %s' % cmd)
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        exit(1)

@cli.command(
    'vcontact',
    context_settings=dict(ignore_unknown_options=True),
    short_help='run vcontact'
)

@click.option('-i',
    '--reference-contigs',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    help='output directory',
    default='None'
)
@click.option('-C',
    '--vcontact-reference-genomes',
    type=str,
    help='output directory',
    default="None"
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='max # of jobs allowed in parallel.',
)
@click.option(
    '-n',
    '--dryrun',
    is_flag=True,
    default=False,
    show_default=True,
    help='Check rules to run and files to produce',
)
@click.option(
    '-r',
    '--reason',
    is_flag=True,
    default=False,
    show_default=True,
    help='Print the reason for each executed rule.',
)
@click.argument(
    'snakemake_args',
    nargs=-1,
    type=click.UNPROCESSED,
)
def runVcontact2(reference_contigs, vcontact_reference_genomes, jobs, dryrun, reason,snakemake_args):
    '''Running vContact on user contigs.
    '''

    DEFAULT_CONFIG = get_default_config()

    cmd = (
        'cp {reference_contigs} {fasta_tot} \n'
        'cp {reference_contigs} {fasta_filtered} \n'
        'snakemake --use-conda -p runVcontact2 '
        '--config '
            'representative_contigs={fasta_tot} '
            'reference_genomes_vcontact={vcontact_reference_genomes} '
        '--jobs {jobs} --rerun-incomplete --latency-wait 600 {dryrun} {reason}'
        ' {args} '

    ).format(
        reference_contigs=reference_contigs,
        fasta_tot=os.path.dirname(os.path.abspath(reference_contigs)) + "/"+ reference_contigs.split("/")[-1].split(".")[0] + ".tot.fasta",
        fasta_filtered=os.path.dirname(os.path.abspath(reference_contigs)) + "/filtered_" + reference_contigs.split("/")[-1].split(".")[0] + ".tot.fasta",
        vcontact_reference_genomes=vcontact_reference_genomes,
        jobs=jobs,
        dryrun='--dryrun' if dryrun else '',
        reason='--reason' if reason else '',
        args=' '.join(snakemake_args),
    )
    logging.info('Executing: %s' % cmd)
    try:
        subprocess.run(cmd, check=True, shell=True)
    except subprocess.CalledProcessError as e:
        exit(1)

#!/usr/bin/perl -w

# Written by Alejandro Reyes
# El programa solo necesita llamarlo desde la carpeta donde están todos los cov_stats y el único parametro que lleva es el nombre del archivo de salida.
# Todos los archivos deben seguir el siguiente esquema: toss_SAMPLEID.covstats.txt y all_SAMPLEID.covstats.txt
# Donde SAMPLEID es el nombre de cada muestra.
# Por último asume que son PE 2x150, si eso cambia toca cambiar el código. Puedo poner una variable que el usuario dé, porque no hay forma de detectar eso desde los covstats.

use strict;
use List::Util qw(max min);

if (@ARGV != 3) {
    die "\nUsage: Make-vOTU-RPKM-Norm.pl <covTossFile> <covAllFile> <outfile>\n\n";
}

my ($CovTossFile, $CovAllFile, $outfile) = @ARGV;

open(OUT, ">$outfile") or die ("Couldn't open outfile $outfile\n");
print OUT "Contig\tRPKM\tCovRatio\tNeededRatio\tExpectedRatio\n";

&process_cov($CovTossFile, $CovAllFile);

close OUT;

sub process_cov {
    my ($CovTossFile, $CovAllFile) = @_;

    open (IN, "<$CovTossFile") or die ("Couldn't open file: $CovTossFile\n");
    open (COVALL, "<$CovAllFile") or die ("Couldn't open file: $CovAllFile\n");

    my $TotalReadsMap = 0;
    my %contig_len = ();
    my %extra_reads = ();
    my %needed_reads = ();
    my %toss_cov = ();
    my %deltaMapBase = ();
    my %contigs_all = ();
    my %contigs_toss = ();
    my %extra = ();
    my %CovRatio = ();
    my %NeededRatio = ();
    my %ExpRatio = ();

    while (my $line = <IN>) {
        chomp $line;
        next if ($line =~ /^#/);
        my @temp = split /\t+/, $line;
        my $name = shift(@temp);
        @{$contigs_toss{$name}} = @temp;
        $TotalReadsMap += $contigs_toss{$name}[6];
        $TotalReadsMap += $contigs_toss{$name}[5];
    }
    close IN;

    while (my $line = <COVALL>) {
        chomp $line;
        next if ($line =~ /^#/);
        my @temp = split /\t+/, $line;
        my $name = shift(@temp);
        die ("Name $name does not exist in Toss\n") unless $contigs_toss{$name};
        @{$contigs_all{$name}} = @temp;
        $contig_len{$name} = $contigs_all{$name}[1];
        $extra_reads{$name} = abs(($contigs_all{$name}[6] + $contigs_all{$name}[5]) - ($contigs_toss{$name}[6] + $contigs_toss{$name}[5]));
        $toss_cov{$name} = ($contigs_toss{$name}[4] > 0) ? ($contigs_toss{$name}[5] + $contigs_toss{$name}[6]) / $contigs_toss{$name}[4] : 0;
        $deltaMapBase{$name} = abs($contigs_all{$name}[4] - $contigs_toss{$name}[4]);
        $needed_reads{$name} = sprintf "%.0f", ($toss_cov{$name} * $deltaMapBase{$name} * 0.9);
        $extra{$name} = min(abs($extra_reads{$name}), abs($needed_reads{$name}));
        $TotalReadsMap += $extra{$name};
        $CovRatio{$name} = log(($contigs_all{$name}[3] + 1) / ($contigs_toss{$name}[3] + 1)) / log(10);
        $NeededRatio{$name} = log(($needed_reads{$name} + 1) / ($extra_reads{$name} + 1)) / log(10);
        $ExpRatio{$name} = ($contigs_all{$name}[3] > 0 && ($contigs_toss{$name}[6] + $contigs_toss{$name}[5]) > 0) ? log((100 * (1 - exp(-1 * (($contigs_toss{$name}[6] + $contigs_toss{$name}[5] + $needed_reads{$name}) * 150) / $contig_len{$name}))) / $contigs_all{$name}[3]) / log(10) : 0;
    }
    close COVALL;

    my $RPKM = 0;

    foreach my $k (keys %contigs_toss) {
        $RPKM = ($contigs_toss{$k}[5] + $contigs_toss{$k}[6] + $extra{$k}) * 1E9 / ($contig_len{$k} * $TotalReadsMap);
        print OUT "$k\t$RPKM\t$CovRatio{$k}\t$NeededRatio{$k}\t$ExpRatio{$k}\n";
    }
}
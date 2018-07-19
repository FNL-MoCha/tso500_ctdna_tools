#!/usr/bin/env perl
# -*- coding: utf-8 -*-
# Read a hotspots BED file (must be Ion formatted) and a MAF file (from the 
# TSO500 ctDNA pipeline) and annotate variants as being MOIs or not based on the
# Hotspots BED file and the NCI-MATCH MOI Rules.
#
# 2018/07/03 - D Sims
################################################################################
use warnings;
use strict;
use autodie;

use Getopt::Long qw( :config bundling auto_abbrev no_ignore_case );
use File::Basename;
use Data::Dump;
use List::Util qw(min max);
use Sort::Versions;
use Term::ANSIColor; # Optional colorized output to terminal.
use Log::Log4perl qw(get_logger);
use DateTime;

use constant DEBUG => 0;

my $tsg_file = "$ENV{'HOME'}/Dropbox/ngs_biofx_files/tso500_ctdna_v1.0/" .
    "TSG_cancergenecensus.txt";

#my $tsg_file = "$ENV{'HOME'}/Dropbox/ngs_biofx_files/tso500_ctdna_v1.0/" .
    #"match_tsgs.txt";

my $scriptname = basename($0);
my $version = "v0.9.071918";
my $description = <<"EOT";
Read in a hotspots BED file (Ion Torrent formatted), and annotate a MAF file 
with the matching variant ID.  Also, determine which variants are MOIs based on 
NCI-MATCH MOI rules.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] <maf_file(s)>
    -m, --mois_only   Only output variants that have passed our MOI rules.
    -b, --bed         Hotspots BED file to use for annotation.
    -v, --version     Version information
    -h, --help        Print this help information
EOT

my $help;
my $ver_info;
my $hs_bed;
my $outfile;
my $mois_only = 0;

GetOptions( 
    "mois_only|m"     => \$mois_only,
    "bed|b=s"       => \$hs_bed,
    "output|o=s"    => \$outfile,
    "version|v"     => \$ver_info,
    "help|h"        => \$help )
or die $usage;

sub help {
	printf "%s - %s\n\n%s\n\n%s\n", $scriptname, $version, $description, $usage;
	exit;
}

sub version {
	printf "%s - %s\n", $scriptname, $version;
	exit;
}

help if $help;
version if $ver_info;

# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 ) {
    print "ERROR: No MAF files passed to script!\n";
    print "$usage\n";
    exit 1;
}

my $outfh;
if ($outfile) {
    print "Writing output to $outfile.\n";
    open($outfh, ">", $outfile);
} else {
    $outfh = \*STDOUT;
}

# Set up a logger.
my $logfile = 'tso500_moi_annotator_' . now('short') . '.log';
my $logger_conf = qq(
    log4perl.logger    = DEBUG, Logfile
    log4perl.appender.Logfile    = Log::Log4perl::Appender::File
    log4perl.appender.Logfile.filename    = $logfile
    log4perl.appender.Logfile.mode    = append
    log4perl.appender.Logfile.layout = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.Logfile.layout.message_chomp_before_newlist = 0
    log4perl.appender.Logfile.layout.ConversionPattern = %d [ %p ]: %m%n
    log4perl.logger.main = DEBUG
);
Log::Log4perl->init(\$logger_conf);
my $logger = get_logger();
$logger->info( "Starting TSO500 MOI Annotation Script..." );

################------ END Arg Parsing and Script Setup ------#################
my @maf_files = @ARGV;
my $bed_data = read_hs_bed($hs_bed);

for my $maf_file (@maf_files) {
    print "Annotating '$maf_file'...\n" if DEBUG;
    $logger->info( "Annotating '$maf_file'..." );
    my $results = annotate_maf($maf_file, $bed_data);
    (my $new_file = $maf_file) =~ s/\.maf/.annotated.maf/;
    $logger->info( "Finished annotating. Printing results..." );
    print_results($results, $new_file, $mois_only);
    $logger->info("Done with $maf_file!\n");
}

sub print_results {
    my ($data, $filename, $filter) = @_;
    my $filter_status;
    ($filter) ? ($filter_status = 'on') : ($filter_status = 'off');

    $logger->info( "Writing results to $filename (filter is $filter_status)");
    print "Writing results to $filename (filter is $filter_status)\n";

    open(my $outfh, ">", $filename);
    select $outfh;

    # Print headers
    print shift @$data;
    print shift @$data;

    for my $var (@$data) {
        next if $var =~ /\.$/ and $filter;
        print $var;
    }
}

sub annotate_maf {
    my ($maf, $hotspots) = @_;
    my @results;

    open(my $fh, "<", $maf);
    # will have two header lines to deal with before we get to the data.
    my $header1 = <$fh>;
    push(@results, $header1);

    # Add moi_type to header
    chomp(my $header2 = <$fh>);
    push(@results, "$header2\tmoi_type\n");

    my $var_count = 0;
    my %moi_count = (
        'Hotspots' => 0,
        'Deleterious' => 0,
        'EGFR Inframe Indel' => 0,
        'ERBB2 Inframe Indel' => 0,
        'KIT Exons 9, 11, 13, or 14 Activating' => 0,
    );

    while (<$fh>) {
        chomp;
        $var_count++;
        my @elems = split(/\t/);
        my $moi_type = '.';
        
        if (DEBUG) {
            print "testing $elems[4]:$elems[5]-$elems[6];$elems[10]>$elems[12]\n";
        }

        # Try to get a hotspot ID
        my $hsid = map_variant($elems[4], $elems[5], $elems[6], $elems[10], 
            $elems[12], $hotspots);

        if (DEBUG) {
            print "> $elems[37]($elems[0]):$elems[34]:$elems[36] maps to  ", 
                "==> $hsid\n";
        }
        $elems[110] = $hsid;

        if ($hsid ne '.') {
            $moi_type = 'Hotspot';
            $moi_count{'Hotspots'}++;
        } else {
            $moi_type = run_nonhs_rules($elems[0], $elems[38], $elems[50], 
                \%moi_count);
        }

        if (DEBUG) {
            print "MOI category: $moi_type\n";
            print "-"x75;
            print "\n\n";
        }
        push(@results, join("\t", @elems, "$moi_type\n"));
    }
    $logger->info("Total variants in file: $var_count." );
    my $moi_count_string;
    for (sort keys %moi_count) {
        $moi_count_string .= sprintf("\t%-37s: %s\n", $_, $moi_count{$_});
    }
    $logger->info("MOI Counts:\n$moi_count_string" );
    return \@results;
}

sub read_tsgs {
    open(my $fh, "<", $tsg_file);
    return [map{ chomp; $_} <$fh>];
}

sub run_nonhs_rules {
    my ($gene, $location, $function, $moi_count) = @_;

    my $tsgs = read_tsgs();
    my $moi_type = '.';
    my $exon = (split(/\//, $location))[0];

    print "incoming => gene: $gene, exon: $exon, function: $function\n" if DEBUG;

    # TODO: add splice variant to list.
    if (grep $gene eq $_, @$tsgs and $function =~ /stop|frameshift/) {
        $moi_count->{'Deleterious'}++;
        return 'Deleterious in TSG';
    }
    elsif ($gene eq 'EGFR') {
        if ($exon eq '19' and $function eq 'inframe_deletion') {
            $moi_count->{'EGFR Inframe Indel'}++;
            return 'EGFR inframe deletion in Exon 19';
        }
        elsif ($exon eq '20' and $function eq 'inframe_insertion') {
            $moi_count->{'EGFR Inframe Indel'}++;
            return 'EGFR inframe insertion in Exon 20';
        }
    }
    elsif ($gene eq 'ERBB2' 
        and $exon eq '20' 
        and $function eq 'inframe_insertion') {
            $moi_count->{'ERBB2 Inframe Indel'}++;
            return 'ERBB2 inframe insertion in Exon 20';
    }
    elsif ($gene eq 'KIT') {
        if ((grep { $exon eq $_  } ('9', '11', '13', '14'))
            && $function =~ /inframe.*/ || $function eq 'missense_variant') {
            $moi_count->{'KIT Exons 9, 11, 13, or 14 Activating'}++;
            return 'KIT Mutation in Exons 9, 11, 13, or 14';
        }
    }
    return $moi_type;
}

sub map_variant {
    my ($chr, $maf_start, $maf_end, $ref, $alt, $hotspots) = @_;
    my $varid = '.';

    if ($hotspots->{$chr}) {
        for my $range (keys %{$hotspots->{$chr}}) {
            my ($r1, $r2) = split(/-/, $range);
            if ($maf_start >= $r1 and $maf_end <= $r2) {
                $varid = match_variant($hotspots->{$chr}{$range}, $maf_start, 
                    $maf_end, $ref, $alt);
               last;
            }
        }
    }
    return $varid;
}

sub match_variant {
    my ($vars, $maf_start, $maf_end, $ref, $alt) = @_;
    for my $var (@$vars) {
        my ($chr, $hs_start, $hs_end, $hs_ref, $hs_alt, $hsid) = split(/;/, $var);
        # MAF file is 1-based, while the HS BED file is 0-based. Also will have 
        # to use different position for mapping if indel compared to snv.
        my $pos;
        ($ref eq '-') ? ($pos = $maf_start) : ($pos = $maf_start-1);
        if ($pos == $hs_start and $ref eq $hs_ref and $alt eq $hs_alt) {
            return $hsid;
        }
    }
    return '.';
}

sub read_hs_bed {
    my $bedfile = shift;
    my @hotspots;
    my %positions;
    open(my $fh, '<', $bedfile);
    my $header = <$fh>;
    while (<$fh>) {
        my ($chr, $start, $end, $id, $allele, $amp) = split(/\t/);
        my ($ref, $alt) = $allele =~ /REF=([ACTG]+)?;OBS=([ACTG]+)?(?:;.*)?/;
        $ref //= '-';
        $alt //= '-';
        push(@{$positions{$chr}}, $start, $end);
        push(@hotspots, [$chr, $start, $end, $ref, $alt, $id]);
    }
    close $fh;
    return build_hs_table(\%positions, \@hotspots);

}

sub build_hs_table {
    my ($positions, $hotspots) = @_;
    my %ranged_hs_table;
    my $bin_width = 10000000;

    for my $chr (sort { versioncmp($a, $b) } keys %$positions) {
        my $min = min(@{$$positions{$chr}});
        my $max = max(@{$$positions{$chr}});
        my $floor = round($min, $bin_width, 'floor');

        while ($floor < $max) {
            my $range = sprintf( "%s-%s", $floor+1, $floor += $bin_width);
            $ranged_hs_table{$chr}->{$range} = [];
        }
    }

    # Insert hotspots where they belong in our hash.
    for my $var (@$hotspots) {
        if ($ranged_hs_table{$var->[0]}) {
            insert_var($ranged_hs_table{$var->[0]}, $var);
        }
    }
    return \%ranged_hs_table;
}

sub insert_var {
    my ($regions, $var) = @_;
    for my $range (keys %$regions) {
        my ($start, $end) = split(/-/, $range);
        if ($var->[1] > $start and $var->[2] < $end) {
            push(@{$$regions{$range}}, join(';', @$var)) and return;
        }
    }
}

sub round {
    # Currently just a floor routine, but can add a ceil too if needed.
    my ($num, $bin, $op) = @_;
    if ($op eq 'floor') {
        return int($num/$bin) * $bin;
    }
}
    
sub now {
    my $format = shift;
    my $now = DateTime->now;
    ($format eq 'short') 
        ? return $now->strftime('%Y%m%d')
        : return $now->strftime('%Y%m%d %H:%m:%s');
}

sub __commify {
    # From Perl Cookbook.  Just for dev and testing to make visualization easier.
    my $number = reverse shift;
    $number =~ s/(\d{3})(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $number;
}

sub __exit__ {
    # Better exit routine for dev and debugging purposes.
    my ($line, $msg) = @_;
    $msg //= '';
    print "\n\n";
    print colored("Got exit message at line $line with message: $msg", 
        'bold white on_green');
    print "\n";
    exit;
}

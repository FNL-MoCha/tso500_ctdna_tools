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
use Text::CSV;

use constant DEBUG => 1;
use constant TRUE => 1;
use constant FALSE => 0;

my $POPULATION_THRESHOLD = 0.05;

my $version = "v0.24.061819";

my $scriptdir = dirname($0);

# Default lookup files.
my $tsg_file = "$scriptdir/resource/tsg_gene_list.txt";
my $hs_bed = "$scriptdir/resource/mocha_tso500_ctdna_hotspots_v1.072018.bed";
my $oncokb_file = "$scriptdir/resource/oncokb_lookup.txt";

for my $resource_file ($tsg_file, $hs_bed, $oncokb_file) {
    die "ERROR: Can not locate necessary resource file '$resource_file'! ",
        "Check that this file is in your package.\n" unless -e $resource_file;
}

my $scriptname = basename($0);
my $description = <<"EOT";
Read in a hotspots BED file (Ion Torrent formatted) or an OncoKB variants file 
(preferred!), and annotate a MAF file with the matching variant ID.  Also,
determine which variants are MOIs based on NCI-MATCH MOI rules.
EOT

my $usage = <<"EOT";
USAGE: $scriptname [options] -a <annotation_method> <maf_file(s)>
    Filtering Options:
    -a, --annot        Method to use for annotation. Select from "hs_bed" or 
                       "oncokb".
    -b, --bed          Hotspots BED file to use for annotation. DEFAULT: $hs_bed
    -o, --oncokb_file  OncoKB file to use for annotation. DEFAULT: $oncokb_file

    Output Options
    -m, --mois_only    Only output variants that have passed our MOI rules.
    -t, --trim_file    Output a trimmed data file in addition to the MAF file.
                       Will only contain some of the most basic output data.

    Other Options
    -v, --version      Version information
    -h, --help         Print this help information
EOT

my $help;
my $ver_info;
my $outfile;
my $annot_method;
my $mois_only = 0;
my $trim_file = 1;
my $verbose;

GetOptions( 
    "annot|a=s"     => \$annot_method,
    "mois_only|m"   => \$mois_only,
    "bed|b=s"       => \$hs_bed,
    "oncokb|o=s"    => \$oncokb_file,
    "version|v"     => \$ver_info,
    "trim_file|t"   => \$trim_file,
    "Verbose|V"     => \$verbose,
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
 
# Add some color output info
my $err   = colored("ERROR:", 'bold red on_black');
my $warn  = colored("WARN:", 'bold yellow on_black');
my $info  = colored("INFO:", 'bold cyan on_black');
my $debug = colored("DEBUG:", 'bold magenta on_black');

# Make sure enough args passed to script
if ( scalar( @ARGV ) < 1 ) {
    print "$err No MAF files passed to script!\n";
    print "$usage\n";
    exit 1;
}

unless ($annot_method) {
    print "$err You must choose a method to use for annotation of these ",
        "data!\n\n";
    print "$usage\n";
    exit 1;
}

$verbose = 1 if DEBUG;

# Set up a logger.
my $logfile = 'tso500_moi_annotator_' . now('short') . '.log';
my $logger_conf = qq(
    log4perl.logger = DEBUG, Logfile, Screen
    log4perl.logger.main = DEBUG

    log4perl.appender.Logfile    = Log::Log4perl::Appender::File
    log4perl.appender.Logfile.filename    = $logfile
    log4perl.appender.Logfile.mode    = append
    log4perl.appender.Logfile.layout = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.Logfile.layout.message_chomp_before_newlist = 0
    log4perl.appender.Logfile.layout.ConversionPattern = %d [ %p ]: %m{indent=4}%n
);

my $extra_conf = qq(
    log4perl.appender.Screen = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr = 0
    log4perl.appender.Screen.layout = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.Screen.layout.message_chomp_before_newlist = 0
    log4perl.appender.Screen.layout.ConversionPattern = %d [ %p ]: %m{indent=4}%n
);

$logger_conf .= $extra_conf if $verbose;

Log::Log4perl->init(\$logger_conf);
my $logger = get_logger();
my $intro_str = sprintf("\n%s\n\t\t    Starting TSO500 MOI Annotation Script\n%s\n",
    "="x80, "="x80);
$logger->info($intro_str);


################------ END Arg Parsing and Script Setup ------#################

# Load up annotation data from either hotspots BED or oncokb file.
my $annotation_data;
if ($annot_method eq 'hs_bed') {
    print "got here!\n";
    $logger->info("Using Hotspot BED file " . basename($hs_bed));
    $annotation_data = read_hs_bed($hs_bed);
} else {
    $logger->info("Using OncoKB file " . basename($oncokb_file));
    $annotation_data = read_oncokb_file($oncokb_file);
}

# dd $annotation_data;
# exit;

# Load up TSGs for the non-hs rules
$logger->info("Using TSG file " . basename($tsg_file));
my $tsgs = read_tsgs($tsg_file);

# Process each MAF file
for my $maf_file (@ARGV) {
    print "Annotating '$maf_file'...\n" if DEBUG;
    $logger->info( "Annotating '$maf_file'..." );
    my $results;
    # XXX
    # $results = annotate_maf($maf_file, $annotation_data, $tsgs);
    $results = annotate_maf2($maf_file, $annotation_data, $tsgs);

    # XXX
    dd $results;
    __exit__(__LINE__, '');

    # Print results.
    (my $new_file = $maf_file) =~ s/\.maf/.annotated.maf/;
    $logger->info( "Finished annotating. Printing results..." );
    print_results($results, $new_file, $mois_only, $trim_file);
    $logger->info("Done with $maf_file!\n\n");
}

sub read_oncokb_file {
    my $oncokb_file = shift;
    my %data;
    open(my $fh, "<", $oncokb_file);
    my $okb_version = (split(' ', readline($fh)))[1];
    $logger->info("OncoKB lookup file version: v$okb_version.");
    my $header = <$fh>;
    while (<$fh>) {
        chomp(my @fields = split(/\t/));
        $data{$fields[0]}->{$fields[2]} = [@fields[3..5]];
    }
    return \%data;
}
    
sub print_results {
    my ($data, $filename, $filter, $trim_file) = @_;
    my $filter_status;
    ($filter) ? ($filter_status = 'on') : ($filter_status = 'off');


    $logger->info( "Writing results to $filename (filter is $filter_status)");
    print "Writing results to $filename (filter is $filter_status)\n";

    open(my $outfh, ">", $filename);

    my $trimfh;
    if ($trim_file) {
        $logger->info( "Also generating a trimmed output file.");
        (my $trim_outfile = $filename) =~ s/\.maf/_trimmed.tsv/;
        open($trimfh, ">", $trim_outfile);
    };

    # Print headers
    print {$outfh} shift @$data;

    my @header = split(/\t/, shift @$data);
    my %header_elems = map{ $_ => $header[$_] } 0..$#header;
    print {$outfh} join("\t", @header);

    my @trim_indices = qw(0 4 5 6 8 9 10 12 13 15 34 35 37 85 71 72 76 99 123 
        133 137 139 140);

    # Get header for trim file if we want one.
    print {$trimfh} join("\t", @header_elems{@trim_indices});

    for my $var (@$data) {
        next if $var =~ /\.$/ and $filter;
        print {$outfh} "$var\n";
        if ($trim_file) {
            my @trim_data = (split(/\t/, $var))[@trim_indices];
            print {$trimfh} join("\t", @trim_data), "\n";
        }
    }
}

sub annotate_maf2 {
    my ($maf, $hotspots, $tsgs) = @_;
    my @results;

    my $csv = Text::CSV->new({ sep_char => "\t" });

    open(my $fh, "<", $maf);
    my $maf_ver = $csv->getline($fh);
    my $header = $csv->getline($fh);

    if ($annot_method eq 'hs_bed') {
        push(@$header, qw(MOI_Type Count));
    }
    elsif ($annot_method eq 'oncokb') {
        push(@$header, qw(MOI_Type Oncogenicity Effect));
    }
    my ($var_count, $filter_count) = 0;
    # TODO: Redo this. Have a set list of categories or something and just have
    # the hash generated as categories are filed rather than this kludgy hard
    # coded way.
    my %moi_count = (
        'Hotspots' => 0,
        'Deleterious' => 0,
        'EGFR Inframe Indel' => 0,
        'ERBB2 Inframe Indel' => 0,
        'KIT Exons 9, 11, 13, 14, or 17 Mutations' => 0,
        'TP53 DBD Mutations' => 0,
        'PIK3CA Exon 20 Mutations' => 0,
        'MED12 Exons 1 or 2 Mutations' => 0,
        'CALR C-terminal truncating' => 0,
        'CALR Exon 9 indels' => 0,
        'NOTCH1 Truncating Mutations' => 0,
        'NOTCH2 Truncating Mutations' => 0,
        'CCND1 Truncating Mutations' => 0,
        'CCND3 Truncating Mutations' => 0,
        'PPM1D Truncating Mutations' => 0,
    );

    while (my $elems = $csv->getline($fh)) {
        $var_count++;
        print("\n", "-"x75, "\n") if DEBUG;

        my %var_data;
        @var_data{@$header} = @$elems;

        next unless $var_data{'Hugo_Symbol'} eq 'MTOR';

        # Filter out SNPs, Intronic Variants, etc.
        $filter_count++ and next unless filter_var2(\%var_data, 'gnomad');
        next;
    }
    __exit__(__LINE__, "");

=cut
        my $moi_type = '.';
        my ($gene, $chr, $start, $end, $ref, $alt, $hgvs_c,$hgvs_p, $tscript_id, 
            $exon, $function) = @elems[0,4,5,6,10,12,34,36,37,38,50];

        # Try to get a hotspot ID
        my ($hsid, $oncogenicity, $effect); 
        # Annotate with a Hotspots BED file
        if ($annot_method eq 'hs_bed') {

            print "testing $chr:$start-$end:$ref>$alt\n" if DEBUG;
            $hsid = map_variant_hs($chr, $start, $end, $ref, $alt, 
                $hotspots);
        } 
        # Annotate with an OncoKB Lookup file. 
        else {
            print "testing $gene:$hgvs_p\n" if DEBUG;
            ($hsid, $oncogenicity, $effect) = map_variant_oncokb($gene, $hgvs_p,
                $hotspots);
        }
=cut
}

sub annotate_maf {
    my ($maf, $hotspots, $tsgs) = @_;
    my @results;

    open(my $fh, "<", $maf);
    # will have two header lines to deal with before we get to the data.
    my $header1 = <$fh>;
    push(@results, $header1);

    # Add moi_type to header
    chomp(my $header2 = <$fh>);
    if ($annot_method eq 'hs_bed') {
        push(@results, "$header2\tMOI_Type\tCount\n");
    } 
    elsif ($annot_method eq 'oncokb') {
        push(@results, "$header2\tMOI_Type\tOncogenicity\tEffect\n");
    }

    my ($var_count, $filter_count) = 0;
    my %moi_count = (
        'Hotspots' => 0,
        'Deleterious' => 0,
        'EGFR Inframe Indel' => 0,
        'ERBB2 Inframe Indel' => 0,
        'KIT Exons 9, 11, 13, 14, or 17 Mutations' => 0,
        'TP53 DBD Mutations' => 0,
        'PIK3CA Exon 20 Mutations' => 0,
        'MED12 Exons 1 or 2 Mutations' => 0,
        'CALR C-terminal truncating' => 0,
        'CALR Exon 9 indels' => 0,
        'NOTCH1 Truncating Mutations' => 0,
        'NOTCH2 Truncating Mutations' => 0,
        'CCND1 Truncating Mutations' => 0,
        'CCND3 Truncating Mutations' => 0,
        'PPM1D Truncating Mutations' => 0,
    );

    while (<$fh>) {
        $var_count++;
        print("\n", "-"x75, "\n") if DEBUG;

        chomp(my @elems = split(/\t/));

        # Filter out SNPs, Intronic Variants, etc.
        $filter_count++ and next unless filter_var(\@elems);

        my $moi_type = '.';
        my ($gene, $chr, $start, $end, $ref, $alt, $hgvs_c,$hgvs_p, $tscript_id, 
            $exon, $function) = @elems[0,4,5,6,10,12,34,36,37,38,50];

        # Try to get a hotspot ID
        my ($hsid, $oncogenicity, $effect); 
        # Annotate with a Hotspots BED file
        if ($annot_method eq 'hs_bed') {

            print "testing $chr:$start-$end:$ref>$alt\n" if DEBUG;
            $hsid = map_variant_hs($chr, $start, $end, $ref, $alt, 
                $hotspots);
        } 
        # Annotate with an OncoKB Lookup file. 
        else {
            print "testing $gene:$hgvs_p\n" if DEBUG;
            ($hsid, $oncogenicity, $effect) = map_variant_oncokb($gene, $hgvs_p,
                $hotspots);
        }

        if (DEBUG) {
            print "> $tscript_id($gene):$hgvs_c:$hgvs_p: maps to ==> $hsid\n\n";
        }

        $elems[110] = $hsid;
        if ($hsid ne '.') {
            $moi_type = 'Hotspot';
            $moi_count{'Hotspots'}++;
        } 
        else {
            ($moi_type, $oncogenicity, $effect) = run_nonhs_rules($gene, $exon, 
                $hgvs_p, $function, \%moi_count, $tsgs);
        }

        if (DEBUG) {
            print "MOI category: $moi_type\n";
            print "-"x75;
            print "\n\n";
        }
        ($annot_method eq 'hs_bed') 
            ? push(@results, join("\t", @elems, $moi_type))
            : push(@results, join("\t", @elems, $moi_type, $oncogenicity, 
                $effect));
    }

    $logger->info(sprintf("Variant Summary\nTotal variants:\t%3s\nFiltered out:\t%3s\n" .
            "Retained:\t\t%3s", $var_count, $filter_count, $var_count-$filter_count));
    my $moi_count_string;
    for (sort keys %moi_count) {
        $moi_count_string .= sprintf("%-42s: %s\n", $_, $moi_count{$_});
    }
    $logger->info("MOI Counts:\n$moi_count_string" );

    return \@results;
}

sub filter_var2 {
    # Remove any variants that are SNPs, Intronic, etc.  May have already been
    # filtered out prior to getting here, but nevertheless the buck stops here!
    my ($variant, $pop_filter) = @_;

    # Functional annotation filter. 
    my $var_class = $variant->{'Variant_Classification'};
    if (grep { /$var_class/} qw(Intron UTR Silent Flank IGR)) {
        $logger->debug(sprintf("Filtering out variant %s because it's '%s'", 
            __gen_hgvs($variant, 'hgvs_c')->[0], $var_class));
        return FALSE;
    }

    # Give some options to the population data you can use to filter. 
    my %pop_fields = (
        'gnomad'  => [qw(gnomAD_AF gnomAD_AFR_AF gnomAD_AMR_AF gnomAD_ASJ_AF
            gnomAD_EAS_AF gnomAD_FIN_AF gnomAD_NFE_AF gnomAD_OTH_AF 
            gnomAD_SAS_AF)],
        'exac'    => [qw(ExAC_AF ExAC_AF_Adj ExAC_AF_AFR ExAC_AF_AMR ExAC_AF_EAS
            ExAC_AF_FIN ExAC_AF_NFE ExAC_AF_OTH ExAC_AF_SAS)],
        '1000G'   => [qw(1000G_ALL 1000G_AFR 1000G_AMR 1000G_EAS 1000G_EUR
            1000G_SAS)],
    );

    # Get the population frequency data for the set requested, and sort them
    # largest to smallest for comparison.
    my @pop_freqs = sort {$b <=> $a } 
        map{ ($variant->{$_}) ? $variant->{$_} : 0 } @{$pop_fields{$pop_filter}};
    if ($pop_freqs[0] > $POPULATION_THRESHOLD) {
        $logger->debug(sprintf("Filtering out variant %s because %s frequency ",
                "too high.", __gen_hgvs($variant, 'hgvs_c')->[0], $pop_filter));
        return FALSE;
    }
    return TRUE;
}

# TODO: REmove me
sub filter_var {
    # Remove any variants that are SNPs, Intronic, etc.  May have already been
    # filtered out prior to getting here, but nevertheless the buck stops here!
    my $variant = shift;

    # Functional annotation filter. 
    if (grep { /$variant->[8]/ } qw(Intron UTR Silent Flank)) {
        $logger->debug(sprintf("Filtering out variant %s because it's '%s'", 
            __gen_hgvs($variant, 'hgvs_c')->[0], $variant->[8]));
        return FALSE;
    }

    # GnomAD Filter
    my @afs = sort { $b <=> $a } 
        map {($variant->[$_]) ? $variant->[$_] : 0} (124..131); # gnomAD index
    if ($afs[0] > 0.05) {
        $logger->debug(sprintf("Filtering out variant %s because GnomAD Freq too". 
            " high.", __gen_hgvs($variant, 'hgvs_c')->[0]));
        return FALSE;
    }
    return TRUE;
}

sub map_variant_oncokb {
    my ($gene, $hgvsp, $lookup_data) = @_;
    my $oncogenicity = '.';
    my $effect = '.';
    my $varid = '.';
    if (exists $lookup_data->{$gene}) {
        if (exists $lookup_data->{$gene}{$hgvsp}) {
            ($oncogenicity, $effect, $varid) = @{$lookup_data->{$gene}{$hgvsp}};
        }
    }
    return ($varid, $oncogenicity, $effect);
}

sub read_tsgs {
    open(my $fh, "<", $tsg_file);
    readline($fh); # Throw out header.
    my $tsg_version = (split(' ', readline($fh)))[1];
    $logger->info("TSG lookup version: $tsg_version");
    return [map{ chomp; $_} <$fh>];
}

sub run_nonhs_rules {
    # Look for non-hotspot MOIs in the data. Return after the first hit, even
    # though some variants may fit into more than one category.
    # TODO:
    #     It's not very flexible or maintainable to have these rules hardcoded.
    #     Move this ruleset to a flat JSON file that can be read in and more
    #     easily maintained.
    my ($gene, $location, $hgvs_p, $function, $moi_count, $tsgs) = @_;

    my $moi_type = '.';
    my $exon = (split(/\//, $location))[0];
    my ($aa_start, $aa_end) = $hgvs_p =~ /^p\.[\*A-Z]+(\d+)(?:_[A-Z]+(\d+))?.*/;
    $aa_end //= $aa_start; # only get end if there is a range from indel.

    # DEBUG
    print "incoming => gene: $gene, exon: $exon, function: $function\n" if DEBUG;
    
    # Deleterious / Truncating in TSG
    if (grep $gene eq $_, @$tsgs and $function =~ /stop|frameshift/) {
        if ($gene eq 'NOTCH1' and $aa_end <= 2250) {
            $moi_count->{'NOTCH1 Truncating Mutations'}++;
            return ('NOTCH1 Truncating Mutations', 'Oncogenic', 
                'Likely Loss-of-function');
        }
        elsif ($gene eq 'NOTCH2') {
            if ($aa_end <= 2009) {
                $moi_count->{'NOTCH2 Truncating Mutations'}++;
                return ('NOTCH2 Truncating Mutations', 'Likely Oncogenic', 
                    'Likely Loss-of-function');
            }
            elsif ($aa_start > 2009 and $aa_end <= 2471) {
                $moi_count->{'NOTCH2 Truncating Mutations'}++;
                return ('NOTCH2 Truncating Mutations', 'Likely Oncogenic', 
                    'Likely Gain-of-function');
            }
        } else {
            $moi_count->{'Deleterious'}++;
            return ('Deleterious in TSG', 'Likely Onogenic', 
                'Likely Loss-of-function');
        }
    }
    # EGFR Exon 19 inframe del and Exon 20 inframe ins.
    elsif ($gene eq 'EGFR') {
        if ($exon eq '19' and $function eq 'inframe_deletion') {
            $moi_count->{'EGFR Inframe Indel'}++;
            return ('EGFR inframe deletion in Exon 19', 'Oncogenic', 
                'Gain-of-function');
        }
        elsif ($exon eq '20' and $function eq 'inframe_insertion') {
            $moi_count->{'EGFR Inframe Indel'}++;
            return ('EGFR inframe insertion in Exon 20', 'Oncogenic',
                'Gain-of-function');
        }
    }
    # ERBB2 Exon 20 inframe indel
    elsif ($gene eq 'ERBB2' 
        and $exon eq '20' 
        and $function eq 'inframe_insertion') {
            $moi_count->{'ERBB2 Inframe Indel'}++;
            return ('ERBB2 inframe insertion in Exon 20', 'Likely Oncogenic',
                'Likely Gain-of-function');
    }
    # Kit Exons, 9, 11, 13, 14, or 17 mutations.
    elsif ($gene eq 'KIT') {
        if ((grep { $exon eq $_  } ('9', '11', '13', '14', '17'))
            && $function =~ /inframe.*/ || $function eq 'missense_variant') {
            $moi_count->{'KIT Exons 9, 11, 13, 14, or 17 Mutations'}++;
            return ('KIT Mutation in Exons 9, 11, 13, 14, or 17',
                'Likely Oncogenic', 'Likely Gain-of-function');
        }
    }
    # TP53 DBD mutations (AA 102-292).
    elsif ($gene eq 'TP53' and ($aa_start > 102 and $aa_end < 292)) {
        $moi_count->{'TP53 DBD Mutations'}++;
        return ('TP53 DBD Mutations', 'Oncogenic', 'Loss-of-function');
    }
    # CCND1 Truncating
    elsif ($gene eq 'CCND1' && ($aa_start >= 256 and $aa_end <= 286)) {
        $moi_count->{'CCND1 Truncating Mutations'}++;
        return ('CCND1 Truncating Mutations', 'Likely Oncogenic', 
            'Likely Gain-of-function');
    }
    # CCND3 Truncating 
    elsif ($gene eq 'CCND3' && ($aa_start >= 286 and $aa_end <= 292)) {
        $moi_count->{'CCND3 Truncating Mutations'}++;
        return ('CCND3 Truncating Mutations', 'Likely Oncogenic', 
            'Likely Gain-of-function');
    }
    # CXCR4 Truncating.
    elsif ($gene eq 'CXCR4' && ($aa_start >= 332 and $aa_end <= 352)) {
        $moi_count->{'CXCR4 Truncating Mutations'}++;
        return ('CXCR4 Truncating Mutations', 'Likely Oncogenic', 
            'Likely Gain-of-function');
    }
    # PPM1D Truncating.
    elsif ($gene eq 'PPM1D' && ($aa_end >= 422 and $aa_end <= 605)) {
        $moi_count->{'PPM1D Truncating Mutations'}++;
        return ('PPM1D Truncating Mutations', 'Likely Oncogenic', 
            'Likely Gain-of-function');
    } 
    # PIK3CA Exon 20 mutations.
    elsif ($gene eq 'PIK3CA' and $exon eq 20) {
        $moi_count->{'PIK3CA Exon 20 Mutations'}++;
        return ('PIK3CA Exon 20 Mutations', 'Likely Oncogenic', 
            'Likely Gain-of-function');
    }
    # MED12 Exon1 or 2 mutation
    elsif ($gene eq 'MED12' and ( grep{ $exon eq $_} qw(1 2))) {
        $moi_count->{'MED12 Exons 1 or 2 Mutations'}++;
        my @ret_val = ('MED12 Exons 1 or 2 Mutations', 'Oncogenic');
        ($exon eq '1')
            ? (push(@ret_val, 'Gain-of-function') and return @ret_val)
            : (push(@ret_val, 'Likely Gain-of-function') and return @ret_val);
    }
    # CALR Exon 9 indels and truncating.
    elsif ($gene eq 'CALR') {
        if ($exon eq '9') {
            if ($function =~ /frameshift|stop/) {
                $moi_count->{'CALR C-terminal truncating'}++;
                return ('CALR C-terminal truncating', 'Likely Oncogenic',
                    'Likely Loss-of-function');
            }
            elsif ($function =~ /insert|delet/) {
                $moi_count->{'CALR Exon 9 indels'}++;
                return ("CALR Exon 9 indels", 'Likely Oncogenic', 
                    'Likely Gain-of-function');
            }
        }
    }
    # JAK2 Exon 12 Alterations
    elsif ($gene eq 'JAK2' and $exon eq '12') {
        $moi_count->{'JAK2 Exon 12 Alterations'}++;
        return ("JAK2 Exon 12 alterations", "Likely Oncogenic", "Gain-of-function");
    }
    # MPL Exon 10 mutations.
    elsif ($gene eq 'MPL' and $exon eq '10') {
        $moi_count->{'MPL Exon 10 Mutations'}++;
        return ("MPL Exon 10 mutations", "Likely Oncogenic", 
            "Likely Gain-of-function");
    }
    # FLT3 Kinase Domain mutations.
    elsif ($gene eq 'FLT3' && ($aa_start >= 604 and $aa_end <= 958)) {
        $moi_count->{'FLT3 TyrK mutations'}++;
        return ("FLT3 TyrK Domain mutions", "Likely Oncogenic", 
            "Likely Gain-of-function");
    }
    return ($moi_type, '.', '.');
}

sub map_variant_hs {
    my ($chr, $maf_start, $maf_end, $ref, $alt, $hotspots) = @_;
    my $varid = '.';

    if ($hotspots->{$chr}) {
        for my $range (keys %{$hotspots->{$chr}}) {
            my ($r1, $r2) = split(/-/, $range);
            if ($maf_start >= $r1 and $maf_end <= $r2) {
                $varid = match_variant($hotspots->{$chr}{$range}, 
                    $maf_start, $maf_end, $ref, $alt);
                last;
            }
        }
    }
    return $varid;
}

sub match_variant {
    my ($vars, $maf_start, $maf_end, $ref, $alt) = @_;
    for my $var (@$vars) {
        my ($chr, $hs_start, $hs_end, $hs_ref, $hs_alt, $hsid, 
            $count) = split(/;/, $var);
        # MAF file is 1-based, while the HS BED file is 0-based. Also will have 
        # to use different position for mapping if indel compared to snv.
        my $pos;
        ($ref eq '-') ? ($pos = $maf_start) : ($pos = $maf_start-1);
        # TODO: May need to add an AA mapping algorithm here to match hotspots 
        # with different base changes.
        if ($pos == $hs_start and $ref eq $hs_ref and $alt eq $hs_alt) {
            return $hsid;
        }
    }
    return '.'
}

sub read_hs_bed {
    my $bedfile = shift;
    my @hotspots;
    my %positions;
    open(my $fh, '<', $bedfile);
    my $header = <$fh>;
    while (<$fh>) {
        my ($chr, $start, $end, $id, $allele, $amp) = split(/\t/);

        # If we have a count metric from Rajesh's PublicData file, then include it.
        my $count = '.';
        if ( $amp =~ /COUNT=(\d+)$/ ) {
            $count = $1;
        }

        my ($ref, $alt) = $allele =~ /REF=([ACTG]+)?;OBS=([ACTG]+)?(?:;.*)?/;
        $ref //= '-';
        $alt //= '-';
        
        push(@{$positions{$chr}}, $start, $end);
        push(@hotspots, [$chr, $start, $end, $ref, $alt, $id, $count]);
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

sub __gen_hgvs {
    # Input the MAF variant string, and an optional hgvs output type, and 
    # retrieve a fully formatted HGVS annotation string. Not including the
    # output type will result in a list containg gHGVS, cHGVS, and pHGVS vals.
    # Valid return types are 'hgvs_g', 'hgvs_c', 'hgvs_p'.
    my ($var_elems, $ret_type) = @_;

    my %g_refseq = (
        chr1 => 'NC_000001.10', chr2 => 'NC_000002.11', chr3 => 'NC_000003.11',
        chr4 => 'NC_000004.11', chr5 => 'NC_000005.9', chr6 => 'NC_000006.11',
        chr7 => 'NC_000007.13', chr8 => 'NC_000008.10', chr9  => 'NC_000009.11',
        chr10 => 'NC_000010.10', chr11 => 'NC_000011.9', chr12 => 'NC_000012.11', 
        chr13 => 'NC_000013.10', chr14 => 'NC_000014.8', chr15 => 'NC_000015.9',
        chr16 => 'NC_000016.9', chr17 => 'NC_000017.10', chr18 => 'NC_000018.9',
        chr19 => 'NC_000019.9', chr20 => 'NC_000020.10', chr21 => 'NC_000021.8',
        chr22 => 'NC_000022.10', chrX => 'NC_000023.10', chrY => 'NC_000024.9',
    );
    my @wanted_keys = qw(Hugo_Symbol Chromosome Start_Position Reference_Allele
        Tumor_Seq_Allele2 HGVSp HGVSp_Short RefSeq);
    my %data;

    @data{qw(gene chr start ref alt cds aa refseq)} = @$var_elems{@wanted_keys};

    # TODO: What if we have more than one transcript ID?  How to handle that case?

    # We have a splice variant or something that doesn't have a change.
    $data{'aa'} = 'p.?' unless ($data{'aa'});

    my %hgvs_annots = (
        'hgvs_g' => "$g_refseq{$data{'chr'}}:g.$data{'start'}$data{'ref'}>" .
           "$data{'alt'}",
        'hgvs_c' => "$data{'refseq'}($data{'gene'}):$data{'cds'}",
        'hgvs_p' => "$data{'refseq'}($data{'gene'}):$data{'aa'}",
    );

    ($ret_type) 
        ? return [$hgvs_annots{$ret_type}] 
        : return [@hgvs_annots{qw( hgvs_g hgvs_c hgvs_p )}];
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

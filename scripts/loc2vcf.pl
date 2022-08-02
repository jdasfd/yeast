#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Path::Tiny;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

loc2vcf.pl - convert gene pos to genome loc and output vcf

=head1 SYNOPSIS

    perl loc2vcf.pl -r <region file> -l <list file> -a <aln format file> [options]
      Options:
        --help          -?          brief help message
        --region        -r  STR     region file in bed4 format
        --list          -l  STR     list file in gene-pos-REF-ALT format (1 col)
        --aln           -a  STR     aln format file contains alignment info

    perl loc2vcf.pl -r region.tsv -l gene.lst -a gene.aln.fa

=cut

GetOptions(
    'help|?'        => sub { Getopt::Long::HelpMessage(0) },
    'region|r=s'    => \( my $region_file ),
    'list|l=s'      => \( my $list_file ),
    'aln|a=s'       => \( my $aln_file ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $region_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($region_file)->is_file ) {
    die "Error: can't find file [$region_file]";
}

if ( !defined $list_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($list_file)->is_file ) {
    die "Error: can't find file [$list_file]";
}

if ( !defined $aln_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($aln_file)->is_file ) {
    die "Error: can't find file [$aln_file]";
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

my $loc_num;
my ($chr, $reg_start, $reg_end);
my $genom_loc;

$region_file = path($region_file)->absolute->stringify;
$list_file = path($list_file)->absolute->stringify;
$aln_file = path($aln_file)->absolute->stringify;

open my $r_in, '<', $region_file;
while(<$r_in>){
    chomp;
    ($chr, $reg_start, $reg_end) = (split/\t/, $_)[1,2,3];
}
close $r_in;

open my $a_in, '<', $aln_file;
while(<$a_in>){
    chomp;
    next if /^>/;
    next if /^[ATCG]/;
    my $length = /^(-+)[ATCG]/;
    $loc_num = length($length);
}
close $a_in;

$genom_loc = $reg_start + $loc_num;

open my $l_in, '<', $list_file;
while(<$l_in>){
    chomp;
    my ($ch_loc, $ref, $alt) = (split /-/, $_)[1,2,3];
    $ch_loc = $ch_loc + $genom_loc - 1;
    print "$chr\t$ch_loc\t$ref\t$alt\n";
}
close $l_in;

__END__

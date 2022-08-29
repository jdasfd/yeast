#!/usr/bin/perl
use strict;
use warnings;

use Getopt::Long;
use Path::Tiny;
use Data::Dumper;

#----------------------------------------------------------#
# GetOpt section
#----------------------------------------------------------#

=head1 NAME

loc2vcf.pl - convert gene pos to genome loc and output vcf

=head1 SYNOPSIS

    perl loc2vcf.pl [options]
      Options:
        --help          -?          brief help message
        --blast         -b  STR     blast result
        --tsv           -t  STR     tsv file in gene-pos-REF-ALT-fit-mut format

    perl loc2vcf.pl -r region.tsv -t gene.tsv

    All output is using stdout. Please using pipe or redirect in linux.

=cut

GetOptions(
    'help|?'        => sub { Getopt::Long::HelpMessage(0) },
    'blast|b=s'     => \( my $blast_file ),
    'tsv|t=s'       => \( my $tsv_file ),
    'output|o=s'    => \( my $output_file ),
) or Getopt::Long::HelpMessage(1);

if ( !defined $blast_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($blast_file)->is_file ) {
    die "Error: can't find file [$blast_file]";
}

if ( !defined $tsv_file ) {
    die Getopt::Long::HelpMessage(1);
}
elsif ( !path($tsv_file)->is_file ) {
    die "Error: can't find file [$tsv_file]";
}

#----------------------------------------------------------#
# init
#----------------------------------------------------------#

my %gene_info;
my ($gene, $chr, $start, $end, $strand);

# blast results input and save them to a hash
open my $b_in, '<', $blast_file;

while(<$b_in>){
    chomp;
    ($gene, $chr, $start, $end, $strand) = (split/\t/, $_)[0,2,3,4,5];
    $gene_info{$gene} = [$chr, $start, $end, $strand];
}

close $b_in;

# print Dumper(\%gene_info);

# tsv2vcf according to genome location
open my $t_in, '<', $tsv_file;

while(<$t_in>){
    chomp;
    my ($in_gene, $loc, $ref, $alt, $fit, $type) = split/\t/, $_;
    if ($gene_info{$in_gene} -> [3] eq "+"){
        my $start = $gene_info{$in_gene} -> [1];
        my $genom_loc = $start + $loc - 1;
        my $chr = $gene_info{$in_gene} -> [0];
        print "$chr\t$genom_loc\t$ref\t$alt\t$fit\t$type\n";
    }
    elsif ($gene_info{$in_gene} -> [3] eq "-"){
        my $start = $gene_info{$in_gene} -> [2];
        my $genom_loc = $start - $loc + 1;
        my $ref_re = &rebase($ref);
        my $alt_re = &rebase($alt);
        my $chr = $gene_info{$in_gene} -> [0];
        print "$chr\t$genom_loc\t$ref_re\t$alt_re\t$fit\t$type\n";
    }
}

close $t_in;

#----------------------------------------------------------#
# Subroutines
#----------------------------------------------------------#

sub rebase {
    my $base = shift;
    my $re_base;

    if ($base eq "A"){
        $re_base = "T";
    }
    elsif($base eq "T"){
        $re_base = "A";
    }
    elsif($base eq "G"){
        $re_base = "C";
    }
    elsif($base eq "C"){
        $re_base = "G";
    }
    else{
        $re_base = "N";
    }

    return $re_base;
}

__END__

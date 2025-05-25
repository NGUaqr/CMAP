#!/usr/bin/perl
use strict;
use warnings;
use Data::Dump qw(dump);
use List::Util qw(max min);

my $drugfile = shift;     # eg druglist.txt
my $diseasefile = shift;  # eg diseaselist.txt
my $wdsize = shift;

# Read drug response data
open my $in, "<", $drugfile or die "Cannot open $drugfile: $!";
my $pt_dresponse = {};
while (<$in>) {
    chomp(my ($gene, $logFC) = split /\t/);
    $pt_dresponse->{uc($gene)}->{'logFC'} = $logFC;
}
close $in;

# Rank from big to small for drug response fold change
my $tmprank = 0;
foreach my $gene (sort { $pt_dresponse->{$b}->{'logFC'} <=> $pt_dresponse->{$a}->{'logFC'} } keys %$pt_dresponse) {
    $pt_dresponse->{$gene}->{'rank'} = $tmprank + 1;
    $tmprank++;
}
my $dgene = $tmprank;
# Read disease response data
open $in, "<", $diseasefile or die "Cannot open $diseasefile: $!";
my $head = <$in>;
my $pt_disease = {};
while (<$in>) {
    chomp(my ($gene, $logFC, $other) = split /\t/);
    $pt_disease->{uc($gene)} = $logFC;
}
close $in;

# Rank from big to small for disease response fold change
my @upregulated;
my @downregulated;
foreach my $gene (keys %$pt_disease) {
    if ($pt_disease->{$gene} > 0) {
        push @upregulated, [$gene, $pt_disease->{$gene}];
    } elsif ($pt_disease->{$gene} < 0) {
        push @downregulated, [$gene, $pt_disease->{$gene}];
    }
}
# Sort upregulated genes first by logFC, then by gene name
@upregulated = sort { $b->[1] <=> $a->[1] || $a->[0] cmp $b->[0] } @upregulated;
# Sort downregulated genes first by logFC, then by gene name
@downregulated = sort { $a->[1] <=> $b->[1] || $a->[0] cmp $b->[0] } @downregulated;

# Output sorted genes for debugging
print "Upregulated genes:\n";
foreach my $gene (@upregulated) {
    print "$gene->[0]\t$gene->[1]\n";
}
print "Downregulated genes:\n";
foreach my $gene (@downregulated) {
    print "$gene->[0]\t$gene->[1]\n";
}

# Determine top and bottom genes
my $up_size = int($wdsize * scalar @upregulated);
my $down_size = int($wdsize * scalar @downregulated);

# Ensure at least 1 gene is selected for top and bottom
$up_size = 1 if $up_size < 1;
$down_size = 1 if $down_size < 1;

print "Up size: $up_size, Down size: $down_size\n"; # Debugging output

my @top = @upregulated[0..($up_size - 1)];
my @bottom = @downregulated[0..($down_size - 1)];

# Output top and bottom genes for debugging
print "Top genes:\n";
dump(\@top);
print "Bottom genes:\n";
dump(\@bottom);

# Calculate KS test scores
my $kstop = KStest(\@top, $pt_dresponse);
my $ksbottom = KStest(\@bottom, $pt_dresponse);

my $connscore;
if ($kstop != $ksbottom) {
    $connscore = $kstop - $ksbottom;
} else {
    $connscore = 0;
}
print "$connscore\n";

open my $out, ">", "out_cmap.txt" or die "Cannot open out_cmap.txt: $!";
print $out "$connscore\n";
close $out;

sub KStest {
    my ($ref, $pt_drug) = @_;
    my @topin = @$ref;
    my @topcom;
    for my $i (0..$#topin) {
        if (exists $pt_drug->{$topin[$i][0]}) {
            push @topcom, $topin[$i];
        }
    }
    my @topa;
    my @topb;
    for my $i (0..$#topcom) {
        my $drank = $pt_drug->{$topcom[$i][0]}->{'rank'};
        if (defined $drank) {
            print "Gene: $topcom[$i][0], Drank: $drank\n";
        } else {
            print "Gene: $topcom[$i][0] has no rank in drug data\n";
            next;
        }
		my $a = (($i + 1) / ($#topcom + 1)) - ($drank / $dgene);
        my $b = ($drank / $dgene) - ($i / ($#topcom + 1));
		print "a: $a, b: $b\n";
        push @topa, $a;
        push @topb, $b;
    }

    my $ks;
    if (max(@topa) > max(@topb)) {
        $ks = max(@topa);
    } else {
        $ks = -max(@topb);
    }
    return $ks;
}

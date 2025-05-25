use strict;
use warnings;
use List::Util qw(shuffle);
use POSIX qw(floor);
use File::Slurp;
use File::Path qw(make_path);

# Function to generate random float between min and max
sub rand_float {
    my ($min, $max) = @_;
    return $min + rand($max - $min);
}

# Function to generate random gene names
sub generate_gene_names {
    my ($num_genes) = @_;
    my @gene_names;
    for my $i (1..$num_genes) {
        push @gene_names, "Gene_" . $i;
    }
    return @gene_names;
}

# Generate drug data
sub generate_drug_data {
    my ($num_genes) = @_;
    my @gene_names = generate_gene_names($num_genes);
    my @drug_data;
    for my $gene (@gene_names) {
        my $activity = rand_float(-1, 1);
        push @drug_data, "$gene\t$activity\n";
    }
    return @drug_data;
}

# Generate disease data
sub generate_disease_data {
    my ($num_genes) = @_;
    my @gene_names = generate_gene_names($num_genes);
    my @disease_data;
    push @disease_data, "Gene_Symbol\tlogFC\tpv\n";
    for my $gene (@gene_names) {
        my $logFC = rand_float(-1, 1);
        my $pv = rand_float(0,0.1);
        push @disease_data, "$gene\t$logFC\t$pv\n";
    }
    return @disease_data;
}

# Parameters
my $num_genes_drug = 100; # Number of genes per drug file
my $num_drugs = 50; # Number of drug files
my $num_diseases = 20; # Number of disease files

# Create directories if they don't exist
make_path('druglist-pos');
make_path('diseaselist-deg2funciton');

# Generate drug files
for my $drug (1..$num_drugs) {
    my @drug_data = generate_drug_data($num_genes_drug);
    write_file("druglist-pos/drug_$drug.txt", @drug_data);
}

# Generate disease files
for my $disease (1..$num_diseases) {
    my $num_genes_disease = int(rand(51)) + 50; # Number of genes for disease (50 to 100)
    my @disease_data = generate_disease_data($num_genes_disease);
    write_file("diseaselist-deg2funciton/disease_$disease.txt", @disease_data);
}

print "Drug files written to druglist-pos/\n";
print "Disease files written to diseaselist-deg2funciton/\n";

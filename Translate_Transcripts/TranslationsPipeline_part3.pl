#!/usr/bin/perl
use warnings;
use strict;
my $species = $ARGV[0];
system "cat TRANSLATE_FILES/*.new.pep > $species.pep";
system "cat [0-9]*.cdna > $species.cdna";
system "rm translate* genewise*";
my %seqs;
open my $cdnafile, "<", $species . ".cdna";
my $sid;
while(<$cdnafile>){
        chomp;
        if(/>/){
                $sid = $_;
                next;
        }
        $seqs{$sid} .= $_;
}
close $cdnafile;
my(%g)=('TCA'=>'S','TCC'=>'S','TCG'=>'S','TCT'=>'S','TTC'=>'F','TTT'=>'F','TTA'=>'L','TTG'=>'L','TAC'=>'Y','TAT'=>'Y','TAA'=>'_','TAG'=>'_','TGC'=>'C','TGT'=>'C','TGA'=>'_','TGG'=>'W','CTA'=>'L','CTC'=>'L','CTG'=>'L','CTT'=>'L','CCA'=>'P','CCC'=>'P','CCG'=>'P','CCT'=>'P','CAC'=>'H','CAT'=>'H','CAA'=>'Q','CAG'=>'Q','CGA'=>'R','CGC'=>'R','CGG'=>'R','CGT'=>'R','ATA'=>'I','ATC'=>'I','ATT'=>'I','ATG'=>'M','ACA'=>'T','ACC'=>'T','ACG'=>'T','ACT'=>'T','AAC'=>'N','AAT'=>'N','AAA'=>'K','AAG'=>'K','AGC'=>'S','AGT'=>'S','AGA'=>'R','AGG'=>'R','GTA'=>'V','GTC'=>'V','GTG'=>'V','GTT'=>'V','GCA'=>'A','GCC'=>'A','GCG'=>'A','GCT'=>'A','GAC'=>'D','GAT'=>'D','GAA'=>'E','GAG'=>'E','GGA'=>'G','GGC'=>'G','GGG'=>'G','GGT'=>'G', 'GCN'=>'A', 'CGN'=>'R', 'GGN'=>'G', 'CCN'=>'P', 'TCN'=>'S', 'ACN'=>'T', 'GTN'=>'V');

open my $OUT, ">",  $species . ".pep";
open my $OUT2, ">", $species . ".fixed";
open my $stopfile, ">", $species . ".seqswithstops.txt";
foreach my $seqid (sort keys %seqs){
        my $protein;    
        my $cdna;
        my $codon;
        for(my $i=0;$i<(length($seqs{$seqid})-2);$i+=3){
                $codon=substr($seqs{$seqid},$i,3);
                $codon= uc $codon;
                if (exists $g{$codon}){
                        if($g{$codon} eq "\_"){
                               # $protein .= "\*";
				print $stopfile "$seqid\n";
                        }
                        else{   
                                $protein .= $g{$codon};
                        }
                }
                else{
                        print "Bad codon: $codon\n";
                }
                unless(exists $g{$codon}){
                        if ($codon =~ /N/){
                                $protein .= "X";
                        }
                }
                if (exists $g{$codon} && $g{$codon} eq "\_"){
                        print "STOP is codon: $codon\n";
                }
                else{
                        $cdna .= $codon;
                }
        }
        print $OUT "$seqid\n$protein\n";
        print $OUT2 "$seqid\n$cdna\n";
}
close $OUT;
close $OUT2;
close $stopfile;
system "rm -r TRANSLATE_FILES/ *.set.txt";
system "rm RUNGENEWISE.* [0-9]*.cdna";
system "mv $species.fixed $species.cdna";


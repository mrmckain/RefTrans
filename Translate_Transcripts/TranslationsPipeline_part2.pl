#!/usr/bin/perl -w
use strict;
my %blastpairs;
open my $setfile, "<", $ARGV[0];
while(<$setfile>){
        chomp;
        my ($seqid, $hitid) = split /\s+/;
	$blastpairs{$seqid}=1;
}

my $sets = substr($ARGV[0], 0, index($ARGV[0], "."));

open my $seqfile, "<", $ARGV[1];
my %orthos;
my $geneseqid;
while(<$seqfile>){
	chomp;
    if(/>/){
	if(/\s/){
		/>(.*?)\s.+/;
    		$geneseqid = $1;
	}
	else{
		$_ =~ />(.+)/;
		$geneseqid = $1;
	}
	unless(exists $blastpairs{$geneseqid}){
		$geneseqid=();
	}
    }
    elsif($geneseqid){
    	$orthos{$geneseqid} .= $_;
    }
}
close $seqfile;
%blastpairs=();

open $setfile, "<", $ARGV[0];
while(<$setfile>){
	chomp;
	my ($seqid, $hitid) = split /\s+/;
    
    open my $TEMPCDNA, ">", "TRANSLATE_FILES/$seqid.cdna";
	print $TEMPCDNA ">$seqid\n$orthos{$seqid}\n";
	close $TEMPCDNA;
	
	
    my $tempid = substr($hitid, 4);
    my $temppos = index($tempid, "\|");
    #my $fileid = substr($tempid, 0, $temppos);
    my $geneid = substr($tempid, ($temppos+1));

    my $starter=0;
    my $seq_temp;
    my $stopper=0;
   
    my $temppep;
    open $temppep, "<", $ARGV[3] or die "Can't open protein file\n";
    while (<$temppep>){
    	chomp;
        if ($_ =~ /$geneid\s+/ || $_ =~ /$geneid$/){
        	$starter=1;
            next;
        }
        if ($stopper == 0 && $starter == 1){
        	if ($_ =~ />/){
            	$stopper=1;
            }
            else{
            	$seq_temp .= $_;
            }
        }
    }
    close $temppep;
    
    open my $TEMPOUT, ">", "TRANSLATE_FILES/$seqid.pep";
 	print $TEMPOUT ">$hitid\n$seq_temp\n";
    close $TEMPOUT;
    
    my $cdnafile = "TRANSLATE_FILES/" . $seqid . ".cdna";
	my $pepfile = "TRANSLATE_FILES/" . $seqid . ".pep";

	system "basedir/wise2.2.0/src/bin/genewise $pepfile $cdnafile -both -splice flat -pep -cdna > TRANSLATE_FILES/$seqid.genewise.txt";
	
	open GENE, "<", "TRANSLATE_FILES/$seqid.genewise.txt";
	my (%pepids, $pepid, $pepseq, %cdnaids, $cdnaid, $cdnaseq);
	my $pick = 0;
	my $seqer = 0;
	while(<GENE>){
        chomp;
        if($_ =~ />/ && $_ =~ /\.pep/){
        	$pepid = $_;
        	$pepid =~ s/.pep//g;
        	$seqer = 1;
        	$pick = 2;
        	next;
    	}
    	if($pick == 2 && $_ !~ /\.\[/ && $seqer == 1){
        	if($_ =~ /\/\//){
            	if(exists $pepids{$pepid}){
                	if(length($pepseq)>length($pepids{$pepid})){
                    	$pepids{$pepid} = $pepseq;
                    	$pepseq = ();
                	}
            	}
            	else{
                	$pepids{$pepid} = $pepseq;
                	$pepseq = ();
            	}
        	}
        	else {
            	$pepseq .= $_;
            	next;
        	}
    	}
  		if($_ =~ />/ && $_ =~ /\.\[/ && $_ !~ /\.pep/){
        	$cdnaid = substr($_, 0, index($_, "\.\["));
        	$pick =1;
        	next;
    	}
    	if($pick ==1){
        	if($_ =~ /\/\//){
            	if(exists $cdnaids{$cdnaid}){
                	if(length($cdnaseq)>length($cdnaids{$cdnaid})){
                    	$cdnaids{$cdnaid} = $cdnaseq;
                        $cdnaseq = ();
                    }
                }
                else{
                	$cdnaids{$cdnaid} = $cdnaseq;
                	$cdnaseq = ();
            	}
        	}
        	else {
            	$cdnaseq .= $_;
            	next;
        	}
		}
	}
	
	open my $fcdnaout, ">>", "$sets.cdna";
	for my $cdnasids (sort keys %cdnaids){
    	while($cdnaids{$cdnasids} =~ /intron/){
    		my $ipos = index($cdnaids{$cdnasids},"intron");
    		$cdnaids{$cdnasids} = substr($cdnaids{$cdnasids}, ($ipos+6));
    	}
        print $fcdnaout "$cdnasids\n$cdnaids{$cdnasids}\n";
	}

system "rm TRANSLATE_FILES/$seqid.cdna TRANSLATE_FILES/$seqid.genewise.txt TRANSLATE_FILES/$seqid.pep";
}

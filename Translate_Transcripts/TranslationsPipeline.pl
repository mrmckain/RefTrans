#!/usr/bin/perl -w
use strict;
use Parallel::ForkManager;
my $pm = Parallel::ForkManager->new($ARGV[4]);
system "mkdir TRANSLATE_FILES";

my $blastfile = $ARGV[0];
my $file = $ARGV[1];
my $ortho = $ARGV[2];

my %queryseqs;
my $tempquid;
open my $fileseq, "<", $ARGV[1];
while(<$fileseq> ){
	chomp;
	if(/^>/){
		if(/\s/){
			/>(.*?)\s.+/;	
			$tempquid = $1;
		}
		else{
			$_ =~ />(.+)/;
			$tempquid = $1;
		}
	}
	else{
		$queryseqs{$tempquid}.= $_;
	}
}
for my $quid (keys %queryseqs){
	$queryseqs{$quid}=length($queryseqs{$quid});
}	
my $sets=0;

my %besthit;
my %blastscore;
open my $blast, "<", $blastfile;
open my $bestblast, ">", $blastfile . ".besthits";
my $m = 1;
my $bestcurline;
my %fullseqs;
my $tempquery;
my $tempsubj;
while(<$blast>){
        chomp;
	my @tarray=split/\t/;
        #$tarray[1] =~ s/.p$/.t/g;
        #$tarray[1] =~ s/_P(\d+)$/_T$1/g;
        if($tarray[10] > 1e-10){
		next;
	}
	if($queryseqs{$tarray[0]} > 50){
		my $qstart;
		my $qstop;
                my $qhit = $tarray[7]-$tarray[6];
		if($qhit < 0){
			$qstart=$tarray[7];
			$qstop=$tarray[6];
			$qhit=abs($qhit);
		}
		else{
			$qstart=$tarray[6];
			$qstop=$tarray[7];
		}
		
                my $shit = $tarray[9]-$tarray[8];
	if(!$tempquery){
              #  my $qover = length($query{$tarray[0]});
                $tempquery=$tarray[0];
		my $sover = $queryseqs{$tarray[0]};
                if(exists $fullseqs{$tarray[0]}{$tarray[1]}){
                        for (my $i=$qstart-1; $i <= ($qstop-1); $i++){
                                @{$fullseqs{$tarray[0]}{$tarray[1]}}[$i]++;
                        }
                }
                else{
                        for (my $j=0; $j<$sover; $j++){
                                @{$fullseqs{$tarray[0]}{$tarray[1]}}[$j]=0;
                        }
                        for (my $i=$qstart-1; $i <= ($qstop-1); $i++){
                                @{$fullseqs{$tarray[0]}{$tarray[1]}}[$i]++;
                        }
                }
	}
	else{
		if($tarray[0] eq $tempquery){
		my $sover = $queryseqs{$tarray[0]};
                if(exists $fullseqs{$tarray[0]}{$tarray[1]}){
                        for (my $i=$qstart-1; $i <= ($qstop-1); $i++){
                                @{$fullseqs{$tarray[0]}{$tarray[1]}}[$i]++;
                        }
                }
                else{
                        for (my $j=0; $j<$sover; $j++){
                                @{$fullseqs{$tarray[0]}{$tarray[1]}}[$j]=0;
                        }
                        for (my $i=$qstart-1; $i <= ($qstop-1); $i++){
                                @{$fullseqs{$tarray[0]}{$tarray[1]}}[$i]++;
                        }
                }
		}
		else{
			&besthittest($tempquery);
			delete $fullseqs{$tempquery};
			my $sover = $queryseqs{$tarray[0]};
			for (my $j=0; $j<$sover; $j++){
                                @{$fullseqs{$tarray[0]}{$tarray[1]}}[$j]=0;
                        }
                        for (my $i=$qstart-1; $i <= ($qstop-1); $i++){
                                @{$fullseqs{$tarray[0]}{$tarray[1]}}[$i]++;
                        }
			$tempquery = $tarray[0];
		}
			
        }
}

=item        my ($query, $hit, $pid, $allen, $mis, $gap, $qstart, $qend, $hstart, $hend, $eval, $bit) = split /\s+/;
        if (exists $besthit{$query}){
		if($eval < $blastscore{$query}){

			$besthit{$query} = $hit;
			$blastscore{$query} = $eval;
			$bestcurline=$_;
		}
		else{
			next;
		}
	}
	else{
		$besthit{$query} = $hit;
                $blastscore{$query} = $eval;
		if($bestcurline){
			print $bestblast "$bestcurline\n";
		}
		$bestcurline=$_;
	}
=cut
	
}
close $blast;
sub besthittest{
	my $querid = $_[0];
	my $best_sub;
	my $best_score;
	for my $subid (keys %{$fullseqs{$querid}}){
        	print "$querid\t$subid\n";
		my $sover = $queryseqs{$querid};
        	my $hit_count;
        	for my $count (@{$fullseqs{$querid}{$subid}}){
                	if($count > 0){
                        	$hit_count++;
                	}
        }

        if($hit_count/$sover >= 0.5){
                if($best_sub){
			if($hit_count/$sover > $best_score){
				$best_sub = $subid;
				$best_score = $hit_count/$sover;
			}
		}
		else{
			$best_sub = $subid;
                        $best_score = $hit_count/$sover;
		}
	}
	}
	if($best_sub){
		print $bestblast "$querid\t$best_sub\n";
		$besthit{$querid}=$best_sub;
	}
        
}

#print $bestblast "$bestcurline\n";
my %newhit;
for my $queryid (sort keys %besthit){
	$newhit{$m}{$queryid}=$besthit{$queryid};
	$m++;
}

my $pairs = scalar keys %newhit;
if($pairs%1000 == 0){
        $sets= $pairs/1000;
}
else{
        $sets= (int($pairs/1000)+1);
}
my $j=1;
for (my $i=1; $i <= $sets; $i++){
	open my $setfile, ">", "$i.set.txt";
    	while($j<=$pairs && $j <= $i*1000){
    		for my $qid (sort keys %{$newhit{$j}}){
        		print $setfile "$qid\t$newhit{$j}{$qid}\n";
            		$j++;
        	}
    	}
}
my @tasks = (1..$sets);
TASKS:
for my $task (@tasks){

        $pm->start and next TASKS;
        `perl basedir/TranslationsPipeline_part2.pl $task.set.txt $file $ortho $ARGV[3]`;

        $pm->finish;
}
$pm->wait_all_children;
`basedir/TranslationsPipeline_part3.pl $ortho`;

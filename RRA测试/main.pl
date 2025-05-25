#Programmed by guofeifei
use strict;
use Data::Dump qw(dump);
#读入疾病列表
opendir(INF,'./diseaselist-deg2function') or die;
my @difiles=readdir(INF);
closedir INF;
#读入药物列表
opendir(INF,'./druglist-pos') or die;
my @files=readdir(INF);
closedir INF;

#建立log文件
open(LOG,">log.txt") or die;

foreach my $difile(@difiles){
	if($difile=~/(.*)\.txt$/){
		
		my $di=$1;
		if(-e "result_$di.txt"){
			print "result_$di.txt has been existed\n";
			next;
		}
		my $pt_score;
		foreach my $file(@files){
			if($file=~/txt$/){
				print "$di $file started\n";
				#print $file;
				system("python parallel.py ./druglist-pos/$file ./diseaselist-deg2function/$difile");
				#计算fisher score
				#print "python SS_fisher.py ./druglist-pos/$file ./diseaselist-deg2function/$difile";
				open(IN,"out_fisher.txt") or die;
				$pt_score->{$file}->{'SSfisher'}=<IN>;
				close IN;
				unlink("out_fisher.txt");
				
				#计算corr score
				#print "python SS_corr.py ./druglist-pos/$file ./diseaselist-deg2function/$difile";
				open(IN,"out_corr.txt") or die;
				my $score1=<IN>;
				$score1 =~ s/\n//g;
				$pt_score->{$file}->{'SSpearson_corr'}=$score1;
				$pt_score->{$file}->{'SSspearman_corr'}=<IN>;
				close IN;
				unlink("out_corr.txt");
				
				#计算lincs score
				#print "python SS_LINCS.py ./druglist-pos/$file ./diseaselist-deg2function/$difile";
				open(IN,"out_lincs.txt") or die;
				my $score1=<IN>;
				$score1 =~ s/\n//g;
				$pt_score->{$file}->{'SSlincs'}=$score1;
				$pt_score->{$file}->{'SSlincs_Norm'}=<IN>;
				close IN;
				unlink("out_lincs.txt");
				
				#计算Zhang score
				#print "python SS_ZhangScore.py ./druglist-pos/$file ./diseaselist-deg2function/$difile";
				open(IN,"out_zhang.txt") or die;
				my $score1=<IN>;
				$score1 =~ s/\n//g;
				$pt_score->{$file}->{'SSzhang'}=$score1;
				$pt_score->{$file}->{'SSzhang_Norm'}=<IN>;
				close IN;
				unlink("out_zhang.txt");
				
				#计算CMAP score
				#print "perl SS_CMAP.pl ./druglist-pos/$file ./diseaselist-deg2function/$difile 0.2";
				open(IN,"out_cmap.txt") or die;
				my $score1=<IN>;
				$score1 =~ s/\n//g;
				$pt_score->{$file}->{'SScmap'}=$score1;
				close IN;
				unlink("out_cmap.txt");
			}
		}
		closedir INF;
		print dump $pt_score;
		open(OUT,">./output/result_$di.txt") or die;
		foreach my $drug(sort keys %$pt_score){
			foreach my $method(sort keys %{$pt_score->{$drug}}){
				print OUT "\t$method";
			}
			print OUT "\n";
			last;
		}

		foreach my $drug(sort keys %$pt_score){
			print OUT "$drug";
			foreach my $method(sort keys %{$pt_score->{$drug}}){
				print OUT "\t$pt_score->{$drug}->{$method}";
			}
			print OUT "\n";
		}
		close OUT;	
	}	
}

close LOG;
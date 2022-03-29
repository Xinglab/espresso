use strict;
use threads;
use Getopt::Long;

my $version = 'S_alpha1.2.2';
my ($help, $sam, $in, $list_samples, $fa, $anno, $out, $SJ_bed, $read_num_cutoff, $read_ratio_cutoff, $group_extra_range, $num_thread, $mapq_cutoff, $keep_tmp, 
	$cont_del_max, $chrM, $inserted_cont_cutoff, $SJFS_dist, $SJFS_dist_add);


Getopt::Long::GetOptions (

	'list_samples|L=s'			=>	\$list_samples,
	'fa|F=s'					=>	\$fa,
	'anno|A=s'					=>	\$anno,
	'SJ_bed|B=s'				=>	\$SJ_bed,
	'out|O=s'					=>	\$out,

	'help|H!'					=>	\$help,
	'keep_tmp|K!'				=>	\$keep_tmp,

	'read_ratio_cutoff|R=f'		=>	\$read_ratio_cutoff,
	'read_num_cutoff|N=i'		=>	\$read_num_cutoff,

	'num_thread|T=i'			=>	\$num_thread,
	'group_extra_range|E=i'		=>	\$group_extra_range,
	'mapq_cutoff|Q=i'			=>	\$mapq_cutoff,
	'cont_del_max|C=i'			=>	\$cont_del_max,
	'chrM|M=s'					=>	\$chrM,

	'inserted_cont_cutoff=i'	=>	\$inserted_cont_cutoff,
	'SJFS_dist=i'				=>	\$SJFS_dist,
	'SJFS_dist_add=i'			=>	\$SJFS_dist_add

);

if (defined($help)) {
	print "
Program:  ESPRESSO (Error Statistics PRomoted Evaluator of Splice Site Options)
Version:  $version
Contact:  Yuan Gao <gaoy\@email.chop.edu, gy.james\@163.com>

Usage:    perl ESPRESSO_S.pl -L samples.tsv -F ref.fa -A anno.gtf -O work_dir

Arguments:

    -L, --list_samples
          tsv list of sample(s) (each file in a line with 1st column as sorted 
          BAM/SAM file and 2nd column as sample name; required)
    -F, --fa
          FASTA file of all reference sequences. Please make sure this file is
          the same one provided to mapper. (required)
    -A, --anno
          input annotation file in GTF format (optional)
    -B, --SJ_bed
          input custom reliable splice junctions in BED format (optional; each
          reliable SJ in one line, with the 1st column as chromosome, the 2nd
          column as upstream splice site 0-base coordinate, the 3rd column as 
          downstream splice site and 6th column as strand)
    -O, --out
          work directory (existing files in this directory may be OVERWRITTEN; 
          default: ./)

    -H, --help
          show this help information

    -N, --read_num_cutoff
          min perfect read count for denovo detected candidate splice junctions
          (default: 2)
    -R, --read_ratio_cutoff
          min perfect read ratio for denovo detected candidate splice junctions: 
          Set this as 1 for completely GTF-dependent processing (default: 0)
    
    -C, --cont_del_max
          max continuous deletion allowed; intron will be identified if longer
          (default: 50)
    -M, --chrM
          tell ESPRESSO the ID of mitochondrion in reference file (default:
          chrM)

    -T, --num_thread
          thread number (default: minimum of 5 and sam file number)  
    -Q, --mapq_cutoff
          min mapping quality for processing (default: 1)

";
} elsif ( !defined($list_samples) or !defined($fa) ) {
	print "The following parameter(s) are required:\n";
	if (!defined($list_samples)) {
		print "\t--list_samples/-L";
	} if (!defined($fa)) {
		print "\t--fa/-F";
	}
	print "\nPlease use the --help or -H option to get usage information.\n";
} else {

	$out = '.' if !defined $out;
	$read_num_cutoff = 2 if !defined $read_num_cutoff;
	$read_ratio_cutoff = 0 if !defined $read_ratio_cutoff;


	$group_extra_range = 0 if !defined $group_extra_range;

	$mapq_cutoff = 1 if !defined $mapq_cutoff;
	$inserted_cont_cutoff = 20 if !defined $inserted_cont_cutoff;
	$SJFS_dist = 10 if !defined $SJFS_dist;
	$SJFS_dist_add = 15 if !defined $SJFS_dist_add;
	$num_thread = 5 if !defined $num_thread;

	$chrM = 'chrM' if !defined $chrM;
	$cont_del_max = 50 if !defined $cont_del_max;


	my @warn_reason;
	my @die_reason;

	if ( !-d $out ) {
		push @warn_reason, " Work directory $out does not exist. Make it by myself\n";
		mkdir $out or die "cannot mkdir $out: $!";
	}

	my @samtools_versions = qx(samtools --version | head -1);
	my $tag_samtools = 0;
	if (@samtools_versions > 0) {
		my @line = split /\s+/, $samtools_versions[0];
		my @version_tiny = split /[.]/, $line[1];
		$tag_samtools = 1 if $version_tiny[0] >= 1 and $version_tiny[1] >= 6;
	}
	push @die_reason, "Please make sure samtools 1.6 or higher version is installed and included in \$PATH!\n" if $tag_samtools == 0;

	my %sam_size;
	my (%file_ID, @files_ID_sort, %sam_sample);
	my ($split_sam_first_read, $sam_distr_ref, $sam_or_bam);
		
	open LIST, "<", $list_samples or die "cannot open $list_samples: $!";

	while(<LIST>) {
		chomp;
		my @line = split /\t/;

		if (exists $sam_sample{$line[0]} and $sam_sample{$line[0]} ne $line[1]) {
			push @die_reason, "$line[0] in $list_samples corresponds to multiple samples: $sam_sample{$line[0]}, $line[1]!\n";
		} elsif (-f $line[0] and -r $line[0]) {
			$sam_size{$line[0]} = -s $line[0];
			$sam_sample{$line[0]} = $line[1];

			push @files_ID_sort, $line[0];
			$file_ID{$line[0]} = $#files_ID_sort;
		} else {
			push @die_reason, "$line[0] in $list_samples does not exist or is not readable!\n";
		}
		if ($line[0] !~ /[sb]am$/i) {
			push @die_reason, "Don't know what format $line[0] is!\n";
		}
	}
	close LIST;

	printf  '['.scalar(localtime)."] Calculating how to assign %g files into %g threads\n", scalar(keys %sam_size), $num_thread if scalar(keys %sam_size) > 1;
	$num_thread = &min(scalar(keys %sam_size), $num_thread);
	$sam_distr_ref = &thread_distribution( $num_thread, \%sam_size );

	if ( @warn_reason >= 1 ) {
		print @warn_reason;
	}
	if ( @die_reason >= 1 ) {
		die @die_reason;
	}

	my $list_samples2 = substr($list_samples, rindex($list_samples,'/')+1);
	open LIST2, ">", "$out/${list_samples2}.updated" or die "cannot write $out/${list_samples2}.updated: $!";
	while (my ($file, $ID) = each %file_ID) {
		print LIST2 "$file\t$sam_sample{$file}\t$ID\n";
	}
	close LIST2;

	my %strand_digit = ('+' => 0, '-' => 1);

	my %splicing_signals = ('0' => {'GT' => 1, 'GC' => 1, 'AG' => -1, 'AT' => 2, 'AC' => -2}, #+
		'1' => {'CT' => 1, 'AC' => -1, 'GC' => -1, 'GT' => 2, 'AT' => -2});	#-

	my %splicing_signals_strand = ('GTAG' => '0', 'GCAG' => '0', 'ATAC' => '0', 'CTAC' => '1', 'CTGC' => '1', 'GTAT' => '1');

	my $sep_end_blast = '/';
	my $sep_SJ = ':';

	
	opendir DIR, $out or die "cannot opendir $out: $!";

	for my $file (readdir DIR) {
		if ( -d "$out/$file" and $file =~ /^\d+$/ ) {
			opendir DIR2, "$out/$file" or die "cannot opendir $out/$file: $!";
			for my $file2 (readdir DIR2) {
				if (-f "$out/$file/$file2" and (index($file, 'sam.list')==0 or $file eq 'sj.list' or $file eq 'sam.list' or $file eq 'group.list')){
					unlink "$out/$file/$file2";
				}
			}
			closedir DIR2;
		} elsif (-f "$out/$file" and ($file eq 'SJ_group_all.fa' or $file =~ /SJ_simplified\.list$/)) {
			unlink "$out/$file";
		}
	}
	closedir DIR;

	my %chr_seq;
	my %chr_seq_len;
	my $current_chr;

	print  '[', scalar(localtime), "] Loading reference\n";
	open FA, "<", $fa or die "cannot open $fa: $!";
	while (<FA>) {
		s/\r\n//;
		chomp;
		if (/^>/) {
			my $info = substr($_, 1, length($_)-1);
			($current_chr, undef) = split /\s+/, $info;
		} else {
			$chr_seq{$current_chr} .= uc($_);
			$chr_seq_len{$current_chr} += length($_);
		}
	}
	close FA;

	#open READ_SUMM, ">", $out."/read_summary.txt" or die "cannot write $out/read_summary.txt: $!";

	if (defined $SJ_bed) {
		print  '[', scalar(localtime), "] Checking custom splice junction\n";
		open BED, "<", $SJ_bed or die "cannot open $SJ_bed: $!";
		while(<BED>) {
			chomp;
			my @line = split /\t/;
			if (exists $strand_digit{$line[5]}) {
				if (!exists $chr_seq_len{$line[0]}) {
					die "Chromosome is not found for $_";
				} elsif ( $chr_seq_len{$line[0]} <= $line[2] ) {
					die "Downstream splice site coordinate should be less than chromosome length($chr_seq_len{$line[0]} <= $line[2]): $_";
				} if ( $line[1]<=0 ) {
					die "Upstream splice site coordinate should be positive($line[1] <= 0): $_";
				}
			} else {
				die "Strand is not recognized for $_";
			}
		}
		close BED;
	}


	my @output_titles = ('group_ID','line_num','readID','read_length','flag','chr','start','mapq','end','notSameStrand','mappedGenome','clip_ends','exonIntronRef','IDS_SJ_ref','NM_num','insNumBg','delNumBg','substNumBg','totalNumBg','SJcorSeqRef','readSeq'); #'SAM', 
	my $output_title = join "\t", @output_titles;
	my %output_titles_ID;
	$output_titles_ID{$output_titles[$_]} = $_ for 0 .. $#output_titles;
	my @ths;

	for my $th_ref (@{$sam_distr_ref}) {
		# $th_ref is the reference to the array of distributed sam files in this thread
		my $th = threads -> new({'context' => 'void'}, \&parallel_scan_prep, [$th_ref,1]);
		my $th_id = $th->tid();
		print " Worker $th_id begins to scan: \n @{$th_ref}\n";
		unshift @ths, $th;
	}


	for (@ths) {
		my $th_id = $_->tid();
		$_ -> join();
		print " Worker $th_id finished reporting.\n";
	}
	@ths = ();

	print  '['.scalar(localtime)."] Re-cluster all reads\n";


	my (%SJ_all_info, %SJ_group);
	my (@group_sort_update, @group_all_info, %chr_group, %group_sep_in_all);

	my (%group_sep);
	my @sam_sort = sort {$a <=> $b or $a cmp $b} keys %sam_size;

	for my $sam_file (@sam_sort) {
		my $num = $file_ID{$sam_file};
		push @files_ID_sort, $sam_file;
		my $file_base = substr($sam_file, rindex($sam_file,'/')+1);

		#open OUT, ">", "$out/$file_base.list" or die "cannot write tmp $out/$file_base.list: $!";
		open GROUP, "<", "$out/$num/group.list" or die "cannot open tmp $out/$num/group.list: $!";
		#open SJ, "<", "$out/$file_base.sj.list" or die "cannot write tmp $out/$file_base.sj.list: $!";

		while (<GROUP>) {
			chomp;
			my @line = split /\t/;
			$group_sep{$line[1]}{$num.'_'.$line[0]} = [$line[2], $line[3]];
		}
		close GROUP;
		if (!defined $keep_tmp) {
			unlink "$out/$num/group.list";
		}
	}
	
	while (my ($chr, $group_ref) = each %group_sep) {
		if (keys %{$group_ref} > 0) {
			my @group_sort = sort {${$group_ref}{$a}[0] <=> ${$group_ref}{$b}[0]} keys %{$group_ref};
			push @group_sort_update, [$group_sort[0]];
			push @group_all_info, [${$group_ref}{$group_sort[0]}[0], ${$group_ref}{$group_sort[0]}[1]];
			$group_sep_in_all{$group_sort[0]} = $#group_sort_update;
			push @{$chr_group{$chr}}, $#group_sort_update;
			for my $i (1 .. $#group_sort) {
				if ( ${$group_ref}{$group_sort[$i]}[0] > $group_all_info[-1][1] ) {
					push @group_sort_update, [$group_sort[$i]];
					push @group_all_info, [${$group_ref}{$group_sort[$i]}[0], ${$group_ref}{$group_sort[$i]}[1]];
					push @{$chr_group{$chr}}, $#group_sort_update;
				} elsif (${$group_ref}{$group_sort[$i]}[1] > $group_all_info[-1][1]) {
					push @{$group_sort_update[-1]}, $group_sort[$i];
					$group_all_info[-1][1] = ${$group_ref}{$group_sort[$i]}[1];
				} else {
					push @{$group_sort_update[-1]}, $group_sort[$i];
				}
				$group_sep_in_all{$group_sort[$i]} = $#group_sort_update;
			}
		}
		
	}

	
	for my $sam_file (@sam_sort) {
		my $num = $file_ID{$sam_file};

		my $file_base = substr($sam_file, rindex($sam_file,'/')+1);

		open SJ, "<", "$out/$num/sj.list" or die "cannot open tmp $out/$num/sj.list: $!";

		while (<SJ>) {
			chomp;
			my @line = split /\t/;
			
			if (exists $SJ_all_info{$line[1]}) {
				$SJ_all_info{$line[1]}[3] += $line[5];
				$SJ_all_info{$line[1]}[4] += $line[6];
			} else {
				my $n = $group_sep_in_all{$num.'_'.$line[0]};
				$SJ_all_info{$line[1]} = [$line[2], $line[3], $line[4], $line[5], $line[6]];
				push @{$SJ_group{$n}}, $line[1];
			}
		}
		close SJ;
	}

	my (%anno_SJ);

	if (defined $anno) {
		print  '[', scalar(localtime), "] Loading annotation\n";
		my (%anno_exons, %isoform_info);
		open ANNO, "<", $anno or die "cannot open $anno: $!";
		while(<ANNO>) {
			chomp;
			my @line = split /\t/;
			next if $line[2] ne 'exon';
			my ($current_isoform, $current_gene);

			if($line[8] =~ /transcript_id \"(\S+)\"/) {
				$current_isoform = $1;
			} else {
				die "no transcript_id found in $_";
			}
			if($line[8] =~ /gene_id \"(\S+)\"/) {
				$current_gene = $1;
			} else {
				die "no transcript_id found in $_";
			}
			$isoform_info{$current_isoform} = [$current_gene, $line[6], $line[0]];
			$anno_exons{$current_isoform}{$line[3]-1} = $line[4];
		}
		close ANNO;

		while (my ($isoform, $exon_start_ref) = each %anno_exons) {
			my @exon_start_sort = sort {$a <=> $b} keys %{$exon_start_ref};
			my $strand = $isoform_info{$isoform}[1];
			my $chr = $isoform_info{$isoform}[2];
			for my $i (1 .. $#exon_start_sort){
				my $SJ2_no_strand = $chr.$sep_SJ.${$exon_start_ref}{$exon_start_sort[$i-1]}.$sep_SJ.$exon_start_sort[$i];
				my $SJ2 = $SJ2_no_strand.$sep_SJ.$strand_digit{$strand};
				if (!exists $anno_SJ{$chr}{$SJ2}) {
					$anno_SJ{$chr}{$SJ2} = [$strand_digit{$strand},${$exon_start_ref}{$exon_start_sort[$i-1]},$exon_start_sort[$i]]; #,[$isoform]
				} 
			}
		}
	}
	if (defined $SJ_bed) {
		print  '[', scalar(localtime), "] Loading custom splice junction\n";
		open BED, "<", $SJ_bed or die "cannot open $SJ_bed: $!";
		while(<BED>) {
			chomp;
			my @line = split /\t/;
			my $SJ2 = $line[0].$sep_SJ.$line[1].$sep_SJ.$line[2].$sep_SJ.$strand_digit{$line[5]};
			$anno_SJ{$line[0]}{$SJ2} = [$strand_digit{$line[5]},$line[1],$line[2]] if !exists $anno_SJ{$line[0]}{$SJ2};
		}
		close BED;
	}

	print  '[', scalar(localtime), "] Summarizing annotated splice junctions for each read group\n";

	
	#my %SJ_updated_all;
 
 	open SJ_FA, ">", $out."/SJ_group_all.fa" or die "cannot write $out/SJ_group_all.fa: $!";
	while ( my ($chr, $group_ID_ref) = each %chr_group ) {
		if (length($chr)<=5) {
			open SJ2, ">", "$out/${chr}_SJ_simplified.list" or die "cannot open tmp $out/${chr}_SJ_simplified.list: $!";
		} else {
			open SJ2, ">>", "$out/other_SJ_simplified.list" or die "cannot open tmp $out/other_simplified.list: $!";
		}
		
		my %annotated_SJ_group;
		my @annotated_SJ_chr = sort {$anno_SJ{$chr}{$a}[1] <=> $anno_SJ{$chr}{$b}[1] or $anno_SJ{$chr}{$a}[2] <=> $anno_SJ{$chr}{$b}[2]} keys %{$anno_SJ{$chr}};

		my $start_SJ_index = 0;

		for my $n (@{$group_ID_ref}) {
			#my $n = $group_all_info[]

			my ($group_start, $group_end) = ($group_all_info[$n][0], $group_all_info[$n][1]);
			for my $i ($start_SJ_index .. $#annotated_SJ_chr) {
				if ($anno_SJ{$chr}{$annotated_SJ_chr[$i]}[1] > $group_end) {
					last;
				} elsif ($anno_SJ{$chr}{$annotated_SJ_chr[$i]}[2] >= $group_start and $anno_SJ{$chr}{$annotated_SJ_chr[$i]}[1] <= $group_end) {
					push @{$annotated_SJ_group{$n}}, $annotated_SJ_chr[$i];
				} elsif ($anno_SJ{$chr}{$annotated_SJ_chr[$i]}[2] < $group_start){
					$start_SJ_index = $i;
				}
			}

			my (%recorded_SJ, %SJ_updated);
			for my $SJ (@{$SJ_group{$n}}) {
				my ($chr2, $SJ_start, $SJ_end, $perfect_count, $all_count) = @{$SJ_all_info{$SJ}};
				die "$chr ne $chr2: $SJ" if $chr ne $chr2;
				my $upstream_2nt = substr($chr_seq{$chr}, $SJ_start, 2);
				my $downstream_2nt = substr($chr_seq{$chr}, $SJ_end-2, 2);
				my $tag = 1;
				my %notSameStrand_SJ_freq;
				if (exists $splicing_signals_strand{$upstream_2nt.$downstream_2nt}){
					$notSameStrand_SJ_freq{$splicing_signals_strand{$upstream_2nt.$downstream_2nt}} ++;
				} 
				if (exists $anno_SJ{$chr}{$SJ.$sep_SJ.'0'}){
					$tag = 2;
					$notSameStrand_SJ_freq{'0'} ++;
					$recorded_SJ{$SJ.$sep_SJ.'0'} = 1;
				}
				if (exists $anno_SJ{$chr}{$SJ.$sep_SJ.'1'}){
					$tag = 2;
					$notSameStrand_SJ_freq{'1'} ++;
					$recorded_SJ{$SJ.$sep_SJ.'1'} = 1;
				}
				my ($notSameStrand_SJ);
				
				if (exists $notSameStrand_SJ_freq{'1'} and exists $notSameStrand_SJ_freq{'0'}) {
					$notSameStrand_SJ = 'x';
				} elsif (exists $notSameStrand_SJ_freq{'1'}) {
					$notSameStrand_SJ = '1';
				} elsif (exists $notSameStrand_SJ_freq{'0'}) {
					$notSameStrand_SJ = '0';
				} else {
					$notSameStrand_SJ = 'x';
					$tag = 0;
				}
				my $SJ2 = $SJ.$sep_SJ.$notSameStrand_SJ;
				$SJ_updated{$SJ2} = [$SJ_start, $SJ_end, $notSameStrand_SJ, $perfect_count, $all_count, $upstream_2nt, $downstream_2nt, $tag, "yes", "no", 0];
				#print SJ2 "$n\t$SJ2\t$chr2\t$SJ_start\t$SJ_end\t$notSameStrand_SJ\t$perfect_count\t$all_count\t$upstream_2nt\t$downstream_2nt\tyes\t";
				if (exists $recorded_SJ{$SJ.$sep_SJ.'0'} or $recorded_SJ{$SJ.$sep_SJ.'1'}) {
					#print SJ2 "yes\n";
					$SJ_updated{$SJ2}[-2] = "yes";
				}
				
			}

			for my $SJ2 (@{$annotated_SJ_group{$n}}) {
				if (!exists $recorded_SJ{$SJ2}) {
					my ($notSameStrand_SJ, $SJ_start, $SJ_end) = @{$anno_SJ{$chr}{$SJ2}};
					$SJ_updated{$SJ2} = [$SJ_start, $SJ_end, $notSameStrand_SJ, 0, 0, "TBD", "TBD", 2, "no", "yes", 0];
					#print SJ2 "$n\t$SJ2\t$chr\t$SJ_start\t$SJ_end\t$notSameStrand_SJ\t0\t0\tTBD\tTBD\tno\tyes\n";
				}
			}

			next if (keys %SJ_updated == 0);
			my @sorted_SJ = sort { $SJ_updated{$a}[0] <=> $SJ_updated{$b}[0] } keys %SJ_updated;
			my @sorted_SJ_group_f = ([$sorted_SJ[0]]);
			for my $i ( 1 .. $#sorted_SJ ){
				if ( $SJ_updated{$sorted_SJ[$i]}[0] > $SJ_updated{$sorted_SJ[$i-1]}[0]+$SJFS_dist_add+2*$SJFS_dist ) {
					push @sorted_SJ_group_f, [$sorted_SJ[$i]];
				} else {
					push @{$sorted_SJ_group_f[-1]}, $sorted_SJ[$i];
				}
			}

			my (@sorted_SJ_cluster, @sorted_SJ_cluster_ends);
			for my $SJ_group_ref ( @sorted_SJ_group_f ){
				my @sorted_SJ_group_r = sort { $SJ_updated{$a}[1] <=> $SJ_updated{$b}[1]} @{$SJ_group_ref};
				push @sorted_SJ_cluster, [$sorted_SJ_group_r[0]];
				push @sorted_SJ_cluster_ends, [ $SJ_updated{$sorted_SJ_group_r[0]}[0], $SJ_updated{$sorted_SJ_group_r[0]}[1] ];
				for my $i ( 1 .. $#sorted_SJ_group_r ){
					my $current_SJ_info_ref = $SJ_updated{$sorted_SJ_group_r[$i]};
					if ( ${$current_SJ_info_ref}[1] > $SJ_updated{$sorted_SJ_group_r[$i-1]}[1]+$SJFS_dist_add+2*$SJFS_dist ) {
						push @sorted_SJ_cluster, [$sorted_SJ_group_r[$i]];
						push @sorted_SJ_cluster_ends, [${$current_SJ_info_ref}[0], ${$current_SJ_info_ref}[1]];
					} else {
						push @{$sorted_SJ_cluster[-1]}, $sorted_SJ_group_r[$i];
						$sorted_SJ_cluster_ends[-1][1] = ${$current_SJ_info_ref}[1];
						if ($sorted_SJ_cluster_ends[-1][0] > ${$current_SJ_info_ref}[0]){
							$sorted_SJ_cluster_ends[-1][0] = ${$current_SJ_info_ref}[0];
						}
					}
				}
			}

			my @SJ_cluster_r_index_sort = sort { $sorted_SJ_cluster_ends[$a][1] <=> $sorted_SJ_cluster_ends[$b][1] } (0 .. $#sorted_SJ_cluster);
			#my @SJ_cluster_f_index_sort = sort { $sorted_SJ_cluster_ends[$b][0] <=> $sorted_SJ_cluster_ends[$a][0] } (0 .. $#sorted_SJ_cluster);

			for my $l (0 .. $#SJ_cluster_r_index_sort) {
				my $index_ori = $SJ_cluster_r_index_sort[$l];
				my @SJs_current_cluster = @{$sorted_SJ_cluster[$index_ori]};
				print SJ2 "SJ_cluster\t$n\t$index_ori\t$l\t$chr\t$sorted_SJ_cluster_ends[$index_ori][0]\t$sorted_SJ_cluster_ends[$index_ori][1]\n";
				for my $SJ2 (@SJs_current_cluster) {
					push @{$SJ_updated{$SJ2}}, $index_ori;
					my ($SJ_start, $SJ_end, $notSameStrand_SJ, $perfect_count, $all_count, $upstream_2nt, $downstream_2nt, $tag, $isPutative, $isAnno, $isHighConfidence, $cluster_index_ori) = @{$SJ_updated{$SJ2}};
					if ( $tag == 2 or ($tag == 1 and $perfect_count >= $read_num_cutoff and $perfect_count >= $all_count*$read_ratio_cutoff) ) {
						#$isHighConfidence = 1;
						$SJ_updated{$SJ2}[-2] = 1;
						print SJ_FA ">$SJ2 SJclst:$index_ori: group:$n:\n";
						print SJ_FA substr($chr_seq{$chr}, $SJ_start-$SJFS_dist-$SJFS_dist_add, $SJFS_dist+$SJFS_dist_add).substr($chr_seq{$chr}, $SJ_end, $SJFS_dist+$SJFS_dist_add)."\n";
					}
					print SJ2 "$n\t$SJ2\t";
					print SJ2 "\t$_" for @{$SJ_updated{$SJ2}};
					print SJ2 "\n";
					#$SJ_updated_all{$SJ2} = $SJ_updated{$SJ2};
				}

			}

		}
		close SJ2;
	}

	close SJ_FA;


	%chr_seq = ();
	%anno_SJ = ();
	my @ths2;


	print "$_($file_ID{$_})\n" for keys %file_ID;
	my $x = 0;
	for (keys %group_sep_in_all) {
		$x ++;
		last if $x>5;
		print "$_($group_sep_in_all{$_})\n";
	}
	for my $th_ref (@{$sam_distr_ref}) {
		# $th_ref is the reference to the array of distributed sam files in this thread
		my $th = threads -> new({'context' => 'void'}, \&parallel_scan_prep, [$th_ref,2]);
		my $th_id = $th->tid();
		print " Worker $th_id begins to scan: \n @{$th_ref}\n";
		push @ths2, $th;
	}

	
	for (@ths2) {
		my $th_id = $_->tid();
		$_ -> join();
		print " Worker $th_id finished reporting.\n";
	}
	@ths2 = ();
	print  '[', scalar(localtime), "] ESPRESSO_S finished its work.\n";

	sub parallel_scan_prep {
		my ($files_ref, $scan_num, $ref_file, $ref_group_sep_in_all);
		($files_ref, $scan_num) = @{$_[0]};

		for my $file (@{$files_ref}) {
			if ($file =~ /sam$/i) {
				if ($scan_num == 1) {
					&parallel_scan1([$file, 'sam']);
				} else {
					&parallel_scan2([$file, 'sam']);
				}
				
			} else {
				if ($scan_num == 1) {
					&parallel_scan1([$file, 'bam']);
				} else {
					&parallel_scan2([$file, 'bam']);
				}
				
			}
		}
	}

	sub parallel_scan2 {
		my ($file, $format) = @{$_[0]};

		my $num = $file_ID{$file};
		print "$file\t$num\n";

		my $exit_sort = system("sort -k 3,3 $out/$num/sam.list > $out/$num/sam.list2");
		if ($exit_sort != 0) {
			die "Failed to sort -k 3,3 $out/$num/sam.list. Exit code is $exit_sort";
		}

		my %read_info;
		my ($last_read, @last_read_alignments);
		my $right_length_index = -1;
		open READ2, "<", "$out/$num/sam.list2" or die "cannot write tmp $out/$num/sam.list2: $!";
		while (<READ2>) {
			chomp;
			my @line = split /\t/;
			next if $line[$output_titles_ID{'readID'}] eq 'readID';
			if (defined $last_read and $last_read ne $line[$output_titles_ID{'readID'}]) {
				my @sort_alignment = sort {${$b}[1] <=> ${$a}[1] or ${$b}[2] <=> ${$a}[2]} @last_read_alignments;
				#my @sort_alignment_length = sort {length(${$b}[4]) <=> length(${$a}[4])} @last_read_alignments;
				if (@sort_alignment > 1) {
					$read_info{$last_read} = [$sort_alignment[0][0], $sort_alignment[0][3], $sort_alignment[0][1]]; #read_length, line_num, mappedGenome
				}
				if (length($sort_alignment[0][4]) != $sort_alignment[0][0]) {
					if ($right_length_index >= 0) {
						$read_info{$last_read} = [$sort_alignment[0][0], $sort_alignment[0][3], $sort_alignment[0][1], $last_read_alignments[$right_length_index][4]];
					} else {
						$read_info{$last_read} = [$sort_alignment[0][0], $sort_alignment[0][3], $sort_alignment[0][1], 'short'];
					}
				}
				@last_read_alignments = ();
				$right_length_index = -1;
			}

			push @last_read_alignments, [$line[$output_titles_ID{'read_length'}], $line[$output_titles_ID{'mappedGenome'}], $line[$output_titles_ID{'mapq'}], $line[$output_titles_ID{'line_num'}], $line[$output_titles_ID{'readSeq'}]];
			$right_length_index = $#last_read_alignments if length($last_read_alignments[-1][-1]) == $last_read_alignments[-1][0] and $last_read_alignments[-1][-1] ne 'NA';
			$last_read = $line[$output_titles_ID{'readID'}];
		}
		close READ2;
			if (defined $last_read) {
				my @sort_alignment = sort {${$b}[1] <=> ${$a}[1] or ${$b}[2] <=> ${$a}[2]} @last_read_alignments;
				#my @sort_alignment_length = sort {length(${$b}[4]) <=> length(${$a}[4])} @last_read_alignments;
				if (@sort_alignment > 1) {
					$read_info{$last_read} = [$sort_alignment[0][0], $sort_alignment[0][3], $sort_alignment[0][1]]; #read_length, line_num, mappedGenome
				}
				if (length($sort_alignment[0][4]) != $sort_alignment[0][0]) {
					if ($right_length_index >= 0) {
						$read_info{$last_read} = [$sort_alignment[0][0], $sort_alignment[0][3], $sort_alignment[0][1], $last_read_alignments[$right_length_index][4]];
					} else {
						$read_info{$last_read} = [$sort_alignment[0][0], $sort_alignment[0][3], $sort_alignment[0][1], 'short'];
					}
				}
			}

		#unlink("$out/$num/sam.list2");

		if ($format eq 'sam') {
			open IN, "<", "$file" or die "cannot open $file: $!";
		} else {
			open IN, "samtools view -h $file |" or die "cannot open $file using samtools: $!";
		}	
		while (<IN>) {
			chomp;
			my @line = split /\t/;
			if ( exists $read_info{$line[0]} and defined $read_info{$line[0]}[3] and $read_info{$line[0]}[3] eq 'short' and length($line[9]) == $read_info{$line[0]}[0] ) {
				my $notSameStrand = &ten2b($line[1], 5);
				if ($notSameStrand == 1) {
					$read_info{$line[0]}[3] = &comp_rev($line[9]);
				} else {
					$read_info{$line[0]}[3] = $line[9];
				}
			}
		}
		close IN;

		open READ, "<", "$out/$num/sam.list" or die "cannot open tmp $out/$num/sam.list: $!";
		open READ3, ">", "$out/$num/sam.list3" or die "cannot write tmp $out/$num/sam.list3: $!";
		my ($pre_chr, %SJ_updated_all);
		#@output_titles
		while (<READ>) {
			chomp;
			my @line = split /\t/;
			next if $line[$output_titles_ID{'readID'}] eq 'readID';
			next if exists $read_info{$line[$output_titles_ID{'readID'}]} and $read_info{$line[$output_titles_ID{'readID'}]}[1] != $line[$output_titles_ID{'line_num'}];
			my $chr = $line[$output_titles_ID{'chr'}];
			if ($chr ne $pre_chr) {
				%SJ_updated_all=();
				if (length($line[$output_titles_ID{'chr'}]) <= 5) {
					open SJ2, "<", "$out/${chr}_SJ_simplified.list" or die "cannot open tmp $out/${chr}_SJ_simplified.list: $!";
				} else {
					open SJ2, "<", "$out/other_SJ_simplified.list" or die "cannot open tmp $out/other_SJ_simplified.list: $!";
				}
				while (<SJ2>) {
					chomp;
					my @line_SJ = split /\t/;
					next if $line_SJ[0] eq 'SJ_cluster';
					$SJ_updated_all{$line_SJ[1]} = [$line_SJ[-2], $line_SJ[-1]];
				}

				close SJ2;
			}

			$pre_chr = $chr;
			$line[$output_titles_ID{'group_ID'}] = $group_sep_in_all{$num.'_'.$line[$output_titles_ID{'group_ID'}]};
			if ($line[$output_titles_ID{'SJcorSeqRef'}] ne 'NA') {
				my @SJcorSeq_strings = split ',', $line[$output_titles_ID{'SJcorSeqRef'}];
				for my $i (0 .. $#SJcorSeq_strings) {
					my $SJcorSeq_string = $SJcorSeq_strings[$i];
					my @SJcor_info = split ';', $SJcorSeq_string;
					if (exists $SJ_updated_all{$SJcor_info[1].$sep_SJ.'0'}) {
						$SJcor_info[1] = $SJcor_info[1].$sep_SJ.'0';
					} elsif (exists $SJ_updated_all{$SJcor_info[1].$sep_SJ.'1'}) {
						$SJcor_info[1] = $SJcor_info[1].$sep_SJ.'1';
					} elsif (exists $SJ_updated_all{$SJcor_info[1].$sep_SJ.'x'}) {
						$SJcor_info[1] = $SJcor_info[1].$sep_SJ.'x';
					} else {
						die "did not find $SJcor_info[1] in $_";
					}
					$SJcor_info[2] = $SJ_updated_all{$SJcor_info[1]}[-2];
					push @SJcor_info, $SJ_updated_all{$SJcor_info[1]}[-1];
					$SJcorSeq_strings[$i] = join ';', @SJcor_info;
				}
				
				$line[$output_titles_ID{'SJcorSeqRef'}] = join ',', @SJcorSeq_strings;
			}
			for my $i (0 .. $#line-1) {
				print READ3 "$line[$i]\t";
			}

			if ( $line[$output_titles_ID{'readSeq'}] ne 'NA' ) {
				print READ3 "$line[-1]\n";
			} elsif ( $line[$output_titles_ID{'readSeq'}] eq 'NA' and exists $read_info{$line[$output_titles_ID{'readID'}]} and defined $read_info{$line[$output_titles_ID{'readID'}]}[3] ) {
				print READ3 "$read_info{$line[$output_titles_ID{'readID'}]}[3]\n";
			} elsif ( $line[$output_titles_ID{'readSeq'}] eq 'NA' ) {
				print READ3 "read_not_recorded\n";
			}

		}
		close READ;
		close READ3;

		if (!defined $keep_tmp){
			unlink "$out/$num/sam.list";
			unlink "$out/$num/sam.list2";
		}

	}

	sub parallel_scan1 {
		my ($file, $format, $key_read);
		if (@{$_[0]} == 3) {
			($file, $format, $key_read) = @{$_[0]};
		} else {
			($file, $format) = @{$_[0]};
		}
		my $num = $file_ID{$file};
		my $file_base = substr($file, rindex($file,'/')+1);
		my $pinhead;

		if ($format eq 'sam') {
			open IN, "<", "$file" or die "cannot open $file: $!";
		} else {
			open IN, "samtools view -h $file |" or die "cannot open $file using samtools: $!";
		}
		if (!-d "$out/$num") {
			mkdir "$out/$num" or die "cannot mkdir $out/$num: $!";
		}
		open OUT, ">", "$out/$num/sam.list" or die "cannot write tmp $out/$num/sam.list: $!";
		open GROUP, ">", "$out/$num/group.list" or die "cannot write tmp $out/$num/group.list: $!";
		open SJ, ">", "$out/$num/sj.list" or die "cannot write tmp $out/$num/sj.list: $!";
		#open DEFAULT, ">", "$out/$file.out" or die "cannot write tmp $out/$file.out: $!";
		#print OUT "$output_title\n";

		my (%SJ_read, %group_info, $pre_chr_end, @pre_group);
		my $group_ID = 0;
		my $line_num = 0;

		while (<IN>) {
			chomp;
			$line_num++;
			my @line = split /\t/;
			if (!defined $pinhead) {
				if (defined $key_read and $line[0] eq $key_read) {
					$pinhead = 1;
				} elsif (!defined $key_read and substr($line[0],0,1) ne '@'){
					$pinhead = 1;
				}
			}
			if (defined $pinhead) {
				my $notSameStrand = &ten2b($line[1], 5);

				if ($line[4] >= $mapq_cutoff and $line[2] ne $chrM) {

					my $msidn = &MSIDN_border($line[5], $cont_del_max);
					
					$line[11] =~ /NM:i:(\d+)/;
					my $NM_num = $1;
					die if !exists $chr_seq_len{$line[2]};
					if (${$msidn}{'max_I'} >= $inserted_cont_cutoff){
						next;
					} elsif ($line[3] - 1 + ${$msidn}{'mappedGenomeN'} > $chr_seq_len{$line[2]}) {
						next;
					}

					my $end = $line[3]-1;
					$end += $_ for @{${$msidn}{'exonIntronRef'}};
					my $out;
					$out .= "$line[0]\t${$msidn}{'read_length'}\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$end\t$notSameStrand\t${$msidn}{'mappedGenome'}\t";
					$out .= "$_," for @{${$msidn}{'clip_ends'}};
					$out .= "\t";

					if ( !defined($pre_chr_end) or $line[2] ne ${$pre_chr_end}[0] or $line[3] > ${$pre_chr_end}[2]+$group_extra_range ){

						if (defined($pre_chr_end)){
							$group_ID ++;
							$group_info{$group_ID} = [${$pre_chr_end}[0], ${$pre_chr_end}[1], ${$pre_chr_end}[2]];
							print GROUP "$group_ID\t${$pre_chr_end}[0]\t${$pre_chr_end}[1]\t${$pre_chr_end}[2]\t";
							for my $read(@pre_group) {
								print GROUP "$read,";
							}
							print GROUP "\n";
							if (scalar(keys %SJ_read)>0) {
								while (my ($SJ, $SJ_info_ref) = each %SJ_read) {
									print SJ "$group_ID\t$SJ\t${$SJ_info_ref}[0]\t${$SJ_info_ref}[1]\t${$SJ_info_ref}[2]\t"; #${$SJ_info_ref}[3]\t
									my @perfect_reads = grep {${$SJ_info_ref}[-1]{$_}==1} keys %{${$SJ_info_ref}[-1]};
									printf SJ "%g\t%g\t", scalar(@perfect_reads), scalar(keys %{${$SJ_info_ref}[-1]});
									if (@perfect_reads>0){
										print SJ "$_," for @perfect_reads;
									} else {
										print SJ "NA";
									}
									
									print SJ "\t";
									print SJ "$_," for keys %{${$SJ_info_ref}[-1]};
									print SJ "\n";
								}
							}
						}
						%SJ_read = ();
						@pre_group = ();
						$pre_chr_end = [ $line[2], &max($line[3]-$group_extra_range, 0), &min($end+$group_extra_range, $chr_seq_len{$line[2]}) ];
					} elsif (${$pre_chr_end}[2] < &min($end+$group_extra_range, $chr_seq_len{$line[2]})) {
						${$pre_chr_end}[2] = &min($end+$group_extra_range, $chr_seq_len{$line[2]});
					}
					push @pre_group, $line[0];

					print DEFAULT "$line[0]\n";

					my ($msidn_all);

					my (@IDS_SJ, @SJ_cor_seq);
					my $exon_end = $line[3]-1;
					my $exon_end2 = $line[3]-1;

					my ($insert_num_total, $delete_num_total) = (${$msidn}{'insertion_nt'}, ${$msidn}{'deletion_nt'});
					my ($subst_num_total) = $NM_num-$insert_num_total- $delete_num_total;
					my ($insert_num_SJ, $delete_num_SJ, $substi_num_SJ) = (0,0,0);
					my ($insert_num_bg, $delete_num_bg, $substi_num_bg, $total_bg) = (0,0,0,0);
					for my $j (0 .. $#{${$msidn}{'ID_from_current_SJ'}}){
						
						my @next_exon_ID = @{${$msidn}{'ID_from_current_SJ'}[$j]};
						my @next_exon_M = @{${$msidn}{'M_from_current_SJ'}[$j]};

						for my $ID_info(@next_exon_ID) {
							if (${$ID_info}[1] > $SJFS_dist and ${$ID_info}[2] > $SJFS_dist) {
								if (${$ID_info}[0] > 0) {
									$insert_num_bg += ${$ID_info}[0];
								} else {
									$delete_num_bg -= ${$ID_info}[0];
								}
								$total_bg += abs(${$ID_info}[0]);
							}
						}
						if ($j >= 1){
							$exon_end2 += ${$msidn}{'exonIntronRef'}[($j-1)*2];
							$exon_end2 += ${$msidn}{'exonIntronRef'}[($j-1)*2+1];
						}
						print DEFAULT "(out)SJ:$j\t$exon_end2:";	
						for my $M_info(@next_exon_M) {
							if (${$M_info}[1]+${$M_info}[0] > $SJFS_dist and ${$M_info}[2]+${$M_info}[0] > $SJFS_dist) {
								print DEFAULT "\t";
								print DEFAULT "$_;" for @{$M_info};
								my ($seq_ref, $seq_read);
								my $substr_current_num = 0;
								my $substr_length = ${$M_info}[0];
								if (${$M_info}[1] < $SJFS_dist){
									$substr_length -= $SJFS_dist-${$M_info}[1];
								}

								if (${$M_info}[2] < $SJFS_dist) {
									$substr_length -= $SJFS_dist-${$M_info}[2];
									if ($substr_length > 0) {
										$seq_ref = substr($chr_seq{$line[2]}, $exon_end2+$SJFS_dist, $substr_length);
										$seq_read = substr($line[9], ${$M_info}[3]+($SJFS_dist-${$M_info}[2]), $substr_length);
									} else {
										print DEFAULT "\tNA";
										next;
									}
									
									
								} else {
									if ($substr_length > 0) {
										$seq_ref = substr($chr_seq{$line[2]}, $exon_end2 + ${$M_info}[2], $substr_length);
										$seq_read = substr($line[9], ${$M_info}[3], $substr_length);
									} else {
										print DEFAULT "\tNA";
										next;
									}
									
								}

								$substr_current_num += &hamming_distance("\U$seq_read", "\U$seq_ref");
								$substi_num_bg += $substr_current_num;
								$total_bg += $substr_length;
								print DEFAULT "\t$seq_ref,$seq_read,$substr_current_num";
							}
						}
						print DEFAULT "\n";
						next if $j == 0;

						$exon_end += ${$msidn}{'exonIntronRef'}[($j-1)*2];
						print DEFAULT "SJ:$j\t$exon_end:";			
						my @pre_exon_ID = @{${$msidn}{'ID_from_current_SJ'}[$j-1]};
						my @pre_exon_M = @{${$msidn}{'M_from_current_SJ'}[$j-1]};
						my ($insert_num, $deletion_num, $subst_num) = (0,0,0);
						my $substr_current_num = 0;
						for my $ID_info(@pre_exon_ID) {
							if (${$ID_info}[1] < $SJFS_dist) {
								if (${$ID_info}[0] > 0) {
									$insert_num += ${$ID_info}[0];
								} elsif (${$ID_info}[1] - ${$ID_info}[0] > $SJFS_dist) {
									$deletion_num += $SJFS_dist-${$ID_info}[1];
								} else {
									$deletion_num -= ${$ID_info}[0];
								}
							}
						}
						for my $M_info(@pre_exon_M) {
							if (${$M_info}[1] < $SJFS_dist) {
								print DEFAULT "\t";
								print DEFAULT "$_;" for @{$M_info};
								my ($seq_ref, $seq_read);
								my $substr_current_num = 0;
								if (${$M_info}[1] + ${$M_info}[0] > $SJFS_dist) {
									$seq_ref = substr($chr_seq{$line[2]}, $exon_end-$SJFS_dist, $SJFS_dist-${$M_info}[1]);
									$seq_read = substr($line[9], ${$M_info}[3]+${$M_info}[0]-($SJFS_dist-${$M_info}[1]), $SJFS_dist-${$M_info}[1]);
									$substr_current_num += &hamming_distance("\U$seq_read", "\U$seq_ref");
								} else {
									$seq_ref = substr($chr_seq{$line[2]}, $exon_end - ${$M_info}[1]-${$M_info}[0], ${$M_info}[0]);
									$seq_read = substr($line[9], ${$M_info}[3], ${$M_info}[0]);
									$substr_current_num += &hamming_distance("\U$seq_read", "\U$seq_ref");
								}
								$subst_num += $substr_current_num;
								print DEFAULT "\t$seq_ref,$seq_read,$substr_current_num";
							}
						}
						for my $ID_info(@next_exon_ID) {
							if (${$ID_info}[2] < $SJFS_dist) {
								if (${$ID_info}[0] > 0) {
									$insert_num += ${$ID_info}[0];
								} elsif (${$ID_info}[2] - ${$ID_info}[0] > $SJFS_dist) {
									$deletion_num += $SJFS_dist-${$ID_info}[2];
								} else {
									$deletion_num -= ${$ID_info}[0];
								}
							}
						}

						my $SJ = $line[2].$sep_SJ.$exon_end.$sep_SJ;
						$exon_end += ${$msidn}{'exonIntronRef'}[($j-1)*2+1];
						$SJ .= $exon_end;

						for my $M_info(@next_exon_M) {
							if (${$M_info}[2] < $SJFS_dist) {
								print DEFAULT "\t";
								print DEFAULT "$_;" for @{$M_info};
								my ($seq_ref, $seq_read);
								my $substr_current_num = 0;
								if (${$M_info}[2] + ${$M_info}[0] > $SJFS_dist) {
									$seq_ref = substr($chr_seq{$line[2]}, $exon_end + ${$M_info}[2], $SJFS_dist-${$M_info}[2]);
									$seq_read = substr($line[9], ${$M_info}[3], $SJFS_dist-${$M_info}[2]);
									$substr_current_num += &hamming_distance("\U$seq_read", "\U$seq_ref");
								} else {
									$seq_ref = substr($chr_seq{$line[2]}, $exon_end + ${$M_info}[2], ${$M_info}[0]);
									$seq_read = substr($line[9], ${$M_info}[3], ${$M_info}[0]);
									$substr_current_num += &hamming_distance("\U$seq_read", "\U$seq_ref");
								}
								$subst_num += $substr_current_num;
								print DEFAULT "\t$seq_ref,$seq_read,$substr_current_num";
							}
						}
						print DEFAULT "\n";

						push @SJ_cor_seq, [${$msidn}{'SJ_dist_read'}[$j-1], $SJ, 0];

						push @IDS_SJ, [$insert_num, $deletion_num, $subst_num];

						my $isPerfect = 0;
						if ($insert_num+$deletion_num+$subst_num == 0) {
							$isPerfect = 1;
							#$SJ_read{$SJ}{$line[0]} = 1;
						}
						if (!exists $SJ_read{$SJ}) {
							my @SJ_info = split $sep_SJ, $SJ;
							$SJ_read{$SJ} = [$SJ_info[0], $SJ_info[1], $SJ_info[2], $SJ_info[2]-$SJ_info[1], {$line[0] => $isPerfect}]; #$SJ_info[3], 
						} elsif ($isPerfect == 1 or !exists $SJ_read{$SJ}[-1]{$line[0]}) {
							$SJ_read{$SJ}[-1]{$line[0]} = $isPerfect;
						}
						$insert_num_SJ += $insert_num;
						$delete_num_SJ += $deletion_num;
						$substi_num_SJ += $subst_num;
					}

					$out .= "$_," for @{${$msidn}{'exonIntronRef'}};
					$out .= "\t";

					if (@IDS_SJ > 0) {
						$out .= "${$_}[0];${$_}[1];${$_}[2]," for @IDS_SJ;
					} else {
						$out .= "NA";
					}
					
					$out .= join '', ("\t", $NM_num, "\t", $insert_num_bg, "\t", $delete_num_bg, "\t", $substi_num_bg, "\t", $total_bg, "\t");

					if (@SJ_cor_seq > 0) {
						$out .= "${$_}[0];${$_}[1];${$_}[2]," for @SJ_cor_seq;
					} else {
						$out .= "NA";
					}
					my $tag_read_seq = 'OK';
					if ( ${$msidn}{'read_length'} == length($line[9]) ){
						if ($notSameStrand == 1) {
							print OUT $group_ID+1,"\t$line_num\t$out\t",&comp_rev($line[9]),"\n";
						} else {
							print OUT $group_ID+1,"\t$line_num\t$out\t$line[9]\n";
						}
					} elsif ( ${$msidn}{'read_length'} > length($line[9]) ) {
						print OUT $group_ID+1,"\t$line_num\t$out\tNA\n";
						$tag_read_seq = 'short';
					} else {
						die "${$msidn}{'read_length'}:\t$_\n";
					}
					#push @{$read_info{$line[0]}}, [${$msidn}{'read_length'}, ${$msidn}{'mappedGenome'}, $line[4], $line_num, $tag_read_seq];
				}
			}
		}

		close IN;
		close OUT;


						if (defined($pre_chr_end)){
							$group_ID ++;
							$group_info{$group_ID} = [${$pre_chr_end}[0], ${$pre_chr_end}[1], ${$pre_chr_end}[2]];
							print GROUP "$group_ID\t${$pre_chr_end}[0]\t${$pre_chr_end}[1]\t${$pre_chr_end}[2]\t";
							for my $read(@pre_group) {
								print GROUP "$read,";
							}
							print GROUP "\n";
							if (scalar(keys %SJ_read)>0) {
								while (my ($SJ, $SJ_info_ref) = each %SJ_read) {
									print SJ "$group_ID\t$SJ\t${$SJ_info_ref}[0]\t${$SJ_info_ref}[1]\t${$SJ_info_ref}[2]\t"; #${$SJ_info_ref}[3]\t
									my @perfect_reads = grep {${$SJ_info_ref}[-1]{$_}==1} keys %{${$SJ_info_ref}[-1]};
									printf SJ "%g\t%g\t", scalar(@perfect_reads), scalar(keys %{${$SJ_info_ref}[-1]});
									if (@perfect_reads>0){
										print SJ "$_," for @perfect_reads;
									} else {
										print SJ "NA" for @perfect_reads;
									}
									
									print SJ "\t";
									print SJ "$_," for keys %{${$SJ_info_ref}[-1]};
									print SJ "\n";
								}
							}
						}

		close GROUP;
		close SJ;
		#unlink("$out/$file_base");
		#close DEFAULT;
		
	}

	sub thread_distribution {
		my $t = $_[0];
		my %file_size = %{$_[1]};	# e.g. 'SRR10000' => 1999999
		my @ordered_file = sort { $file_size{$b} <=> $file_size{$a} } (keys %file_size);
		my (@thread_sum, @file_distributed);

		for my $file (@ordered_file) {
			my ($min_thread, undef) = sort { $thread_sum[$a] <=> $thread_sum[$b] } (0 .. $t-1);
			unshift @{$file_distributed[ $min_thread ]}, $file;
			$thread_sum[ $min_thread ] += $file_size{$file};
		}
		my @sort_thread = map {$file_distributed[$_]} (sort {$thread_sum[$b] <=> $thread_sum[$a]} (0 .. $t-1));
		\@sort_thread;
	}

	sub MSIDN_border {
		my @counts = split /[MSIDHN]/, $_[0];
		my $read_length_val = 0;
		my ($not_mapped_num, $mapped_nt_genome, $mapped_nt_genomeN, $insertion_nt, $deletion_nt) = (0,0,0,0);
		my $continuous_deletion_max = $_[1];

		$_[0] =~ s/H/S/g;
		my @styles = split /\d+/, $_[0];
		shift @styles;

		my @exonIntron = (0);
		my @SJ_dist_read = (0);
		my @ID_from_current_SJ;
		my @M_from_current_SJ;
		push @ID_from_current_SJ, [];
		push @M_from_current_SJ, [];
		my $dist_from_prev_SJ = 0;
		my $dist_from_read_start = 0;
		my $max_I = 0;
		my @clip_ends = (0, 0);
		for my $i (0 .. $#styles) {
			if ($styles[$i] eq 'N' or ($counts[$i] > $continuous_deletion_max and $styles[$i] eq 'D')) { # or ($styles[$i] eq 'D' and $counts[$i] >= $continuous_deletion_max)
				push @exonIntron, $counts[$i];
				push @exonIntron, 0;
				push @SJ_dist_read, $SJ_dist_read[-1];
				$dist_from_prev_SJ = 0;
				push @ID_from_current_SJ, [];
				push @M_from_current_SJ, [];
				$mapped_nt_genomeN += $counts[$i];
			} elsif ($styles[$i] eq 'M') {
				$exonIntron[-1] += $counts[$i];
				$SJ_dist_read[-1] += $counts[$i];
				$read_length_val += $counts[$i];
				$mapped_nt_genome += $counts[$i];
				$mapped_nt_genomeN += $counts[$i];
				@{$M_from_current_SJ[-1]} = map {[${$_}[0],${$_}[1]+$counts[$i],${$_}[2],${$_}[3]]} @{$M_from_current_SJ[-1]} if @{$M_from_current_SJ[-1]}>=1;
				push @{$M_from_current_SJ[-1]}, [$counts[$i], 0, $dist_from_prev_SJ, $dist_from_read_start];
				$dist_from_prev_SJ += $counts[$i];
				$dist_from_read_start += $counts[$i];
				@{$ID_from_current_SJ[-1]} = map {[${$_}[0],${$_}[1]+$counts[$i],${$_}[2]]} @{$ID_from_current_SJ[-1]} if @{$ID_from_current_SJ[-1]}>=1;
			} elsif ($styles[$i] eq 'D') {
				$exonIntron[-1] += $counts[$i];
				$mapped_nt_genome += $counts[$i];
				$mapped_nt_genomeN += $counts[$i];
				@{$ID_from_current_SJ[-1]} = map {[${$_}[0],${$_}[1]+$counts[$i],${$_}[2]]} @{$ID_from_current_SJ[-1]} if @{$ID_from_current_SJ[-1]}>=1;
				@{$M_from_current_SJ[-1]} = map {[${$_}[0],${$_}[1]+$counts[$i],${$_}[2],${$_}[3]]} @{$M_from_current_SJ[-1]} if @{$M_from_current_SJ[-1]}>=1;
				$deletion_nt += $counts[$i];
				push @{$ID_from_current_SJ[-1]}, [0-$counts[$i], 0, $dist_from_prev_SJ];
				$dist_from_prev_SJ += $counts[$i];
			} elsif ($styles[$i] eq 'S') {
				$SJ_dist_read[-1] += $counts[$i];
				$read_length_val += $counts[$i];
				$not_mapped_num += $counts[$i];
				$dist_from_read_start += $counts[$i];
				if ($mapped_nt_genome > 0) {
					$clip_ends[1] += $counts[$i];
				} else {
					$clip_ends[0] += $counts[$i];
				}
			} elsif ($styles[$i] eq 'I') {
				$SJ_dist_read[-1] += $counts[$i];
				$read_length_val += $counts[$i];
				$insertion_nt += $counts[$i];
				$max_I = $counts[$i] if $counts[$i] > $max_I;
				push @{$ID_from_current_SJ[-1]}, [$counts[$i], 0, $dist_from_prev_SJ];
				$dist_from_read_start += $counts[$i];
			}
		}
		pop @SJ_dist_read;
		{'exonIntronRef' => \@exonIntron,
		'read_length' => $read_length_val,
		'not_mapped_num' => $not_mapped_num, 
		'mappedGenome' => $mapped_nt_genome,
		'mappedGenomeN' => $mapped_nt_genomeN, 
		'insertion_nt' => $insertion_nt, 
		'deletion_nt' => $deletion_nt, 
		'ID_from_current_SJ' => \@ID_from_current_SJ, 
		'M_from_current_SJ' => \@M_from_current_SJ, 
		'SJ_dist_read' => \@SJ_dist_read,
		'clip_ends' => \@clip_ends,
		'max_I' => $max_I
		};
	}

	sub comp_rev {
		my $seq = reverse($_[0]);
		$seq =~ tr/ATCG/TAGC/;
		$seq;
	}

	sub ten2b {
		my $b_string = sprintf("%b", $_[0]);
		if ($_[1] <= length($b_string)) {
			substr(reverse($b_string), $_[1]-1, 1);
		} else {
			0;
		}
	}

	sub min {
		if ($_[0] < $_[1]) {
			$_[0];
		} else {
			$_[1];
		}
	}

	sub max {
		if ($_[0] > $_[1]) {
			$_[0];
		} else {
			$_[1];
		}
	}

	sub hamming_distance{ length( $_[0] ) - ( ($_[0] ^ $_[1]) =~ tr[\0][\0] ) };

}

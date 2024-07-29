use strict;
use threads;
use threads::shared;
use Thread::Queue;
use Getopt::Long;
use Storable qw(freeze thaw);

use File::Basename qw(dirname);
use lib dirname(__FILE__);
use ESPRESSO_Version;

my $version_number = ESPRESSO_Version::get_version_number();
my $version_string = "S_$version_number";
my ($help, $sam, $in, $list_samples, $fa, $anno, $out, $SJ_bed, $read_num_cutoff, $read_ratio_cutoff, $group_extra_range, $num_thread, $mapq_cutoff, $keep_tmp,
	$cont_del_max, $chrM, $inserted_cont_cutoff, $SJFS_dist, $SJFS_dist_add, $sort_buffer_size, $alignment_read_groups);

my $arguments_before_parsing = "@ARGV";
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
	'SJFS_dist_add=i'			=>	\$SJFS_dist_add,
	'sort_buffer_size=s'	=>	\$sort_buffer_size,
	'alignment_read_groups'	=>	\$alignment_read_groups

);

if (defined($help)) {
	print "
Program:  ESPRESSO (Error Statistics PRomoted Evaluator of Splice Site Options)
Version:  $version_string
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
          thread number. At most 1 thread can be used per input alignment file.
          (default: minimum of 5 and sam file number)
    -Q, --mapq_cutoff
          min mapping quality for processing (default: 1)
    --sort_buffer_size
          memory buffer size for running 'sort' commands (default: 2G)
    --alignment_read_groups
          use overlapping alignment coordinates to determine read groups

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

	if (defined $sort_buffer_size) {
		$sort_buffer_size = "--buffer-size=$sort_buffer_size";
	} else {
		$sort_buffer_size = "--buffer-size=2G";
	}

	my @worker_threads;
	my @thread_input_queues;
	my @thread_output_queues;
	for (1 .. $num_thread) {
		push @thread_input_queues, Thread::Queue->new();
		push @thread_output_queues, Thread::Queue->new();
		my $thread = threads->new({'context' => 'void'}, \&thread_work_loop, [$thread_input_queues[-1], $thread_output_queues[-1]]);
		push @worker_threads, $thread;
	}

	my @warn_reason;
	my @die_reason;

	if ( !-d $out ) {
		push @warn_reason, " Work directory $out does not exist. Make it by myself\n";
		mkdir $out or die "cannot mkdir $out: $!";
	}

	# The summary file will be written as a final step
	my $summary_path = "$out/espresso_s_summary.txt";
	my %summary_data = ();
	&initialize_summary_data(\%summary_data);
	$summary_data{'perl_command'} = "$^X " . __FILE__ . " $arguments_before_parsing";

	my $samtools_version_error = ESPRESSO_Version::check_samtools_version();
	if ($samtools_version_error ne '') {
		push @die_reason, "$samtools_version_error\n";
	}

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

	if (!$alignment_read_groups and !$anno) {
		push @die_reason, "--anno is required unless --alignment_read_groups\n";
	}

	# --gtf_read_groups was added to the code as an optional feature.
	# Now it's the default and --alignment_read_groups enables the old behavior.
	# $alignment_read_groups is only used for argument parsing.
	# $gtf_read_groups is used in the remainder of the code.
	my $gtf_read_groups = !$alignment_read_groups;
	undef $alignment_read_groups;

	print_with_timestamp(sprintf("Calculating how to assign %g files into %g threads", scalar(keys %sam_size), $num_thread)) if scalar(keys %sam_size) > 1;
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

	my %splicing_signals_strand = ('GTAG' => '0', 'GCAG' => '0', 'ATAC' => '0', 'CTAC' => '1', 'CTGC' => '1', 'GTAT' => '1');

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

	print_with_timestamp("Loading reference");
	open FA, "<", $fa or die "cannot open $fa: $!";
	while (<FA>) {
		s/\r\n//;
		chomp;
		if (/^>/) {
			my $info = substr($_, 1, length($_)-1);
			($current_chr, undef) = split /\s+/, $info;
		} elsif (defined $current_chr) {
			$chr_seq{$current_chr} .= uc($_);
			$chr_seq_len{$current_chr} += length($_);
		}
	}
	close FA;
	my $num_chr_seqs = scalar(keys %chr_seq);
	if ($num_chr_seqs == 0) {
		die "No sequence information found in $fa";
	}
	# sort chr names to get consistent output
	my @chr_sort = sort {$a cmp $b} keys %chr_seq;

	# $anno_SJ{$chr}{$SJ} = [$strand, $SJ_start, $SJ_end]
	# where $SJ = chr:start:end:strand
	#       $SJ_start is the 0-based 1st position of the junction
	#       $SJ_end is the 1-based last position of the junction
	my %anno_SJ;
	my %anno_coords_by_chr_by_gene;
	if (defined $anno) {
		print_with_timestamp("Loading annotation");
		&load_annotation_file($anno, $fa, \%anno_coords_by_chr_by_gene, \%anno_SJ, \%summary_data, \%chr_seq, \%strand_digit, $sep_SJ, $gtf_read_groups);
	} else {
		$summary_data{'num_chrs_only_in_fa'} = $num_chr_seqs;
	}

	if (defined $SJ_bed) {
		print_with_timestamp("Loading custom splice junction");
		&load_sj_bed_file($SJ_bed, \%strand_digit, \%chr_seq_len, \%anno_SJ, $sep_SJ);
	}
	$summary_data{'num_annotated_splice_junctions'} += (scalar keys %{$anno_SJ{$_}}) for keys %anno_SJ;

	my @output_titles = ('group_ID','line_num','readID','sample','read_length','flag','chr','start','mapq','end','notSameStrand','mappedGenome','clip_ends','exonIntronRef','IDS_SJ_ref','NM_num','insNumBg','delNumBg','substNumBg','totalNumBg','SJcorSeqRef','readSeq');
	my %output_titles_ID;
	$output_titles_ID{$output_titles[$_]} = $_ for 0 .. $#output_titles;

	{ # limit scope of variables used to pass arguments to threads
		my %inputs_for_all_threads;
		$inputs_for_all_threads{'scan_number'} = 1;
		$inputs_for_all_threads{'file_ID'} = \%file_ID;
		$inputs_for_all_threads{'out'} = $out;
		$inputs_for_all_threads{'mapq_cutoff'} = $mapq_cutoff;
		$inputs_for_all_threads{'chrM'} = $chrM;
		$inputs_for_all_threads{'cont_del_max'} = $cont_del_max;
		$inputs_for_all_threads{'chr_seq_len'} = \%chr_seq_len;
		$inputs_for_all_threads{'inserted_cont_cutoff'} = $inserted_cont_cutoff;
		$inputs_for_all_threads{'group_extra_range'} = $group_extra_range;
		$inputs_for_all_threads{'SJFS_dist'} = $SJFS_dist;
		$inputs_for_all_threads{'chr_seq'} = \%chr_seq;
		$inputs_for_all_threads{'sep_SJ'} = $sep_SJ;
		$inputs_for_all_threads{'sam_sample'} = \%sam_sample;
		my $serialized_inputs_for_all_threads :shared;
		$serialized_inputs_for_all_threads = freeze(\%inputs_for_all_threads);
		%inputs_for_all_threads = ();

		for my $thread_i (0 .. $#{$sam_distr_ref}) {
			my $files_for_thread = $sam_distr_ref->[$thread_i];
			my $serialized_files_for_thread :shared;
			$serialized_files_for_thread = freeze($files_for_thread);
			my %thread_input :shared;
			$thread_input{'files_for_thread'} = $serialized_files_for_thread;
			$thread_input{'inputs_for_all_threads'} = $serialized_inputs_for_all_threads;
			$thread_input_queues[$thread_i]->enqueue(\%thread_input);
			print "Worker $thread_i begins to scan: \n @{$files_for_thread}\n";
		}
	}

	&wait_for_threads_to_complete_work(scalar(@{$sam_distr_ref}), \@thread_input_queues, \@thread_output_queues, \@worker_threads, \%summary_data);

	print_with_timestamp("Re-cluster all reads");

	my (%group_sep);
	my @sam_sort = sort {$a <=> $b or $a cmp $b} keys %sam_size;

	for my $sam_file (@sam_sort) {
		my $num = $file_ID{$sam_file};
		push @files_ID_sort, $sam_file;

		my $group_path = "$out/$num/group.list";
		open GROUP, "<", $group_path or die "cannot open tmp $group_path: $!";
		while (<GROUP>) {
			chomp;
			my @line = split /\t/;
			my $orig_group_id = $line[0];
			my $chr = $line[1];
			my $start = $line[2];
			my $end = $line[3];
			$group_sep{$chr}{$num.'_'.$orig_group_id} = [$start, $end];
		}
		close GROUP;
		if (!defined $keep_tmp) {
			unlink $group_path;
		}
	}

	my (@group_all_info, %chr_group, %group_sep_in_all);
	my $num_read_groups = &finalize_read_groups($gtf_read_groups, \%group_sep, \@group_all_info, \%chr_group, \%group_sep_in_all, \@chr_sort, \%anno_coords_by_chr_by_gene);
	$summary_data{'number_of_read_groups'} = $num_read_groups;

	my (%SJ_all_info, %SJ_group);
	for my $sam_file (@sam_sort) {
		my $num = $file_ID{$sam_file};
		open SJ, "<", "$out/$num/sj.list" or die "cannot open tmp $out/$num/sj.list: $!";
		while (<SJ>) {
			chomp;
			my @line = split /\t/;
			my $group_id = $line[0];
			my $sj = $line[1];
			my $chr = $line[2];
			my $start = $line[3];
			my $end = $line[4];
			my $num_perfect = $line[5];
			my $num_all = $line[6];

			if (exists $SJ_all_info{$sj}) {
				$SJ_all_info{$sj}[3] += $num_perfect;
				$SJ_all_info{$sj}[4] += $num_all;
			} else {
				my @possible_groups = &lookup_final_group($num, $group_id, $chr, $start, $end, \%group_sep_in_all, \@group_all_info, \%chr_group, $gtf_read_groups);
				if (scalar(@possible_groups) > 0) {
					$SJ_all_info{$sj} = [$chr, $start, $end, $num_perfect, $num_all];
				}
				for my $n (@possible_groups) {
					push @{$SJ_group{$n}}, $sj;
				}
			}
		}
		close SJ;
	}

	print_with_timestamp("Summarizing annotated splice junctions for each read group");

 	open SJ_FA, ">", $out."/SJ_group_all.fa" or die "cannot write $out/SJ_group_all.fa: $!";
	while ( my ($chr, $group_ID_ref) = each %chr_group ) {
		if (length($chr)<=5) {
			open SJ2, ">", "$out/${chr}_SJ_simplified.list" or die "cannot open tmp $out/${chr}_SJ_simplified.list: $!";
		} else {
			open SJ2, ">>", "$out/other_SJ_simplified.list" or die "cannot open tmp $out/other_simplified.list: $!";
		}

		my @annotated_SJ_chr = sort {$anno_SJ{$chr}{$a}[1] <=> $anno_SJ{$chr}{$b}[1] or $anno_SJ{$chr}{$a}[2] <=> $anno_SJ{$chr}{$b}[2]} keys %{$anno_SJ{$chr}};

		my %sj_has_been_counted_in_a_group = ();
		my $start_SJ_index = 0;
		for my $n (@{$group_ID_ref}) {
			my @annotated_SJ_group;
			my ($group_start, $group_end) = ($group_all_info[$n][0], $group_all_info[$n][1]);
			for my $i ($start_SJ_index .. $#annotated_SJ_chr) {
				if ($anno_SJ{$chr}{$annotated_SJ_chr[$i]}[1] > $group_end) {
					last;
				} elsif ($anno_SJ{$chr}{$annotated_SJ_chr[$i]}[2] >= $group_start and $anno_SJ{$chr}{$annotated_SJ_chr[$i]}[1] <= $group_end) {
					push @annotated_SJ_group, $annotated_SJ_chr[$i];
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

			for my $SJ2 (@annotated_SJ_group) {
				if (!exists $recorded_SJ{$SJ2}) {
					my ($notSameStrand_SJ, $SJ_start, $SJ_end) = @{$anno_SJ{$chr}{$SJ2}};
					$SJ_updated{$SJ2} = [$SJ_start, $SJ_end, $notSameStrand_SJ, 0, 0, "TBD", "TBD", 2, "no", "yes", 0];
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

			for my $l (0 .. $#SJ_cluster_r_index_sort) {
				my $index_ori = $SJ_cluster_r_index_sort[$l];
				my @SJs_current_cluster = @{$sorted_SJ_cluster[$index_ori]};
				print SJ2 "SJ_cluster\t$n\t$index_ori\t$l\t$chr\t$sorted_SJ_cluster_ends[$index_ori][0]\t$sorted_SJ_cluster_ends[$index_ori][1]\n";
				for my $SJ2 (@SJs_current_cluster) {
					my $already_seen = exists $sj_has_been_counted_in_a_group{$SJ2};
					$sj_has_been_counted_in_a_group{$SJ2} = 1;
					push @{$SJ_updated{$SJ2}}, $index_ori;
					my ($SJ_start, $SJ_end, $notSameStrand_SJ, $perfect_count, $all_count, $upstream_2nt, $downstream_2nt, $tag, $isPutative, $isAnno, $isHighConfidence, $cluster_index_ori) = @{$SJ_updated{$SJ2}};
					if (!$already_seen) {
						$summary_data{'perfect_splice_junction_read_count'} += $perfect_count;
						$summary_data{'total_splice_junction_read_count'} += $all_count;
					}
					if ( $tag == 2 or ($tag == 1 and $perfect_count >= $read_num_cutoff and $perfect_count >= $all_count*$read_ratio_cutoff) ) {
						if (!$already_seen) {
							$summary_data{'num_high_confidence_splice_junctions'} ++;
						}
						$SJ_updated{$SJ2}[-2] = 1;
						print SJ_FA ">$SJ2 SJclst:$index_ori: group:$n:\n";
						print SJ_FA substr($chr_seq{$chr}, $SJ_start-$SJFS_dist-$SJFS_dist_add, $SJFS_dist+$SJFS_dist_add).substr($chr_seq{$chr}, $SJ_end, $SJFS_dist+$SJFS_dist_add)."\n";
					}
					print SJ2 "$n\t$SJ2\t";
					print SJ2 "\t$_" for @{$SJ_updated{$SJ2}};
					print SJ2 "\n";
				}

			}

		}
		close SJ2;
	}

	close SJ_FA;


	%chr_seq = ();
	%anno_SJ = ();

	{ # limit scope of variables used to pass arguments to threads
		my %inputs_for_all_threads;
		$inputs_for_all_threads{'scan_number'} = 2;
		$inputs_for_all_threads{'file_ID'} = \%file_ID;
		$inputs_for_all_threads{'out'} = $out;
		$inputs_for_all_threads{'output_titles_ID'} = \%output_titles_ID;
		$inputs_for_all_threads{'group_sep_in_all'} = \%group_sep_in_all;
		$inputs_for_all_threads{'group_all_info'} = \@group_all_info;
		$inputs_for_all_threads{'chr_group'} = \%chr_group;
		$inputs_for_all_threads{'gtf_read_groups'} = $gtf_read_groups;
		$inputs_for_all_threads{'sep_SJ'} = $sep_SJ;
		$inputs_for_all_threads{'keep_tmp'} = $keep_tmp;
		my $serialized_inputs_for_all_threads :shared;
		$serialized_inputs_for_all_threads = freeze(\%inputs_for_all_threads);
		%inputs_for_all_threads = ();

		for my $thread_i (0 .. $#{$sam_distr_ref}) {
			my $files_for_thread = $sam_distr_ref->[$thread_i];
			my $serialized_files_for_thread :shared;
			$serialized_files_for_thread = freeze($files_for_thread);
			my %thread_input :shared;
			$thread_input{'files_for_thread'} = $serialized_files_for_thread;
			$thread_input{'inputs_for_all_threads'} = $serialized_inputs_for_all_threads;
			$thread_input_queues[$thread_i]->enqueue(\%thread_input);
			print "Worker $thread_i begins to scan: \n @{$files_for_thread}\n";
		}
	}

	&wait_for_threads_to_complete_work(scalar(@{$sam_distr_ref}), \@thread_input_queues, \@thread_output_queues, \@worker_threads, \%summary_data);

	&cleanup_threads_and_exit_if_error(\@thread_input_queues, \@worker_threads);
	&write_summary_file($summary_path, \%summary_data);

	print_with_timestamp("ESPRESSO_S finished its work.");

	sub parallel_scan_prep {
		my ($main_input_ref, $summary_data_ref) = @_;
		my $files_ref = thaw($main_input_ref->{'files_for_thread'});
		my $all_input_ref = thaw($main_input_ref->{'inputs_for_all_threads'});
		my $scan_num = $all_input_ref->{'scan_number'};

		for my $file (@{$files_ref}) {
			if ($file =~ /sam$/i) {
				if ($scan_num == 1) {
					&parallel_scan1([$file, 'sam', $all_input_ref, $summary_data_ref]);
				} else {
					&parallel_scan2([$file, 'sam', $all_input_ref, $summary_data_ref]);
				}

			} else {
				if ($scan_num == 1) {
					&parallel_scan1([$file, 'bam', $all_input_ref, $summary_data_ref]);
				} else {
					&parallel_scan2([$file, 'bam', $all_input_ref, $summary_data_ref]);
				}

			}
		}
	}

	sub print_read3_lines {
		my ($pending_lines_ref, $read3_handle, $summary_data_ref) = @_;

		for my $line (@{$pending_lines_ref}) {
			print $read3_handle $line;
			$summary_data_ref->{'number_of_reads_output'} ++;
		}
		@{$pending_lines_ref} = ();
	}

	sub print_read3_lines_for_previous_group {
		my ($group, $pending_group, $pending_lines_ref, $read3_handle, $summary_data_ref) = @_;

		if ((defined $pending_group) and ($pending_group < $group)) {
			&print_read3_lines($pending_lines_ref, $read3_handle, $summary_data_ref);
			$pending_group = undef;
		}
		return $pending_group;
	}

	# With $gtf_read_groups a read could be in multiple groups if the coordinates
	# are between two genes without overlapping either gene.
	# Output a separate line for each group that the read is assigned to, but
	# only write lines for a group after all lines for previous groups are written.
	# This can lead to a read having multiple lines in the --tsv_compt output file.
	# Also a novel intergenic isoform could be reported in both adjacent groups with
	# different novel IDs.
	sub maybe_print_read3_lines {
		my ($group_id_header_i, $chr, $read_start, $chr_group_ref, $group_all_info_ref, $pending_group, $possible_groups_ref, $pending_lines_ref, $columns_ref, $read3_handle, $summary_data_ref) = @_;

		my $group = $possible_groups_ref->[0];
		$columns_ref->[$group_id_header_i] = $group;
		my $out_line = join("\t", @{$columns_ref});
		$out_line .= "\n";

		my $num_possible = scalar(@{$possible_groups_ref});
		if ($num_possible > 2) {
			die "a read is assigned to more than 2 read groups ($num_possible) $out_line";
		}

		if ($num_possible == 2) {
			my $group_2 = $possible_groups_ref->[1];
			$pending_group = &print_read3_lines_for_previous_group($group_2, $pending_group, $pending_lines_ref, $read3_handle, $summary_data_ref);

			print $read3_handle $out_line;
			$summary_data_ref->{'number_of_reads_output'} ++;

			$columns_ref->[$group_id_header_i] = $group_2;
			my $out_line_2 = join("\t", @{$columns_ref});
			$out_line_2 .= "\n";

			push @{$pending_lines_ref}, $out_line_2;
			$pending_group = $group_2;
			return $pending_group;
		}

		my $can_output = 0;
		my $lowest_chr_group = $chr_group_ref->{$chr}[0];
		if ($group == $lowest_chr_group) {
			$can_output = 1;
		} else {
			my $prev_group = $group - 1;
			my $ends_ref = $group_all_info_ref->[$prev_group];
			my $prev_end = $ends_ref->[1];
			if ($read_start > $prev_end) {
				$can_output = 1;
			}
		}

		if ($can_output) {
			# ($group + 1) is used to print any lines for $group
			$pending_group = &print_read3_lines_for_previous_group($group + 1, $pending_group, $pending_lines_ref, $read3_handle, $summary_data_ref);
			print $read3_handle $out_line;
			$summary_data_ref->{'number_of_reads_output'} ++;
		} else {
			&print_read3_lines_for_previous_group($group, $pending_group, $pending_lines_ref, $read3_handle, $summary_data_ref);
			push @{$pending_lines_ref}, $out_line;
			$pending_group = $group;
		}

		return $pending_group;
	}

	sub parallel_scan2 {
		my ($file, $format, $input_ref, $summary_data_ref) = @{$_[0]};
		my $file_ID_ref = $input_ref->{'file_ID'};
		my $out = $input_ref->{'out'};
		my $output_titles_ID_ref = $input_ref->{'output_titles_ID'};
		my $group_sep_in_all_ref = $input_ref->{'group_sep_in_all'};
		my $group_all_info_ref = $input_ref->{'group_all_info'};
		my $chr_group_ref = $input_ref->{'chr_group'};
		my $gtf_read_groups = $input_ref->{'gtf_read_groups'};
		my $sep_SJ = $input_ref->{'sep_SJ'};
		my $keep_tmp = $input_ref->{'keep_tmp'};

		my $num = $file_ID_ref->{$file};
		print "$file\t$num\n";

		my $sort_command = "sort $sort_buffer_size -k 3,3 $out/$num/sam.list > $out/$num/sam.list2";
		my $exit_sort = system($sort_command);
		if ($exit_sort != 0) {
			die "Failed to $sort_command. Exit code is $exit_sort";
		}

		my %read_info;
		my ($last_read, @last_read_alignments);
		my $right_length_index = -1;
		open READ2, "<", "$out/$num/sam.list2" or die "cannot write tmp $out/$num/sam.list2: $!";
		while (<READ2>) {
			chomp;
			my @line = split /\t/;
			next if $line[$output_titles_ID_ref->{'readID'}] eq 'readID';
			if (defined $last_read and $last_read ne $line[$output_titles_ID_ref->{'readID'}]) {
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

			push @last_read_alignments, [$line[$output_titles_ID_ref->{'read_length'}], $line[$output_titles_ID_ref->{'mappedGenome'}], $line[$output_titles_ID_ref->{'mapq'}], $line[$output_titles_ID_ref->{'line_num'}], $line[$output_titles_ID_ref->{'readSeq'}]];
			$right_length_index = $#last_read_alignments if length($last_read_alignments[-1][-1]) == $last_read_alignments[-1][0] and $last_read_alignments[-1][-1] ne 'NA';
			$last_read = $line[$output_titles_ID_ref->{'readID'}];
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
		open(my $read3_handle, '>', "$out/$num/sam.list3") or die "cannot write tmp $out/$num/sam.list3: $!";
		my ($pre_chr, %SJ_updated_all);
		my @pending_read3_lines = ();
		my $pending_read3_lines_group = undef;
		while (<READ>) {
			chomp;
			my @line = split /\t/;
			next if $line[$output_titles_ID_ref->{'readID'}] eq 'readID';
			next if exists $read_info{$line[$output_titles_ID_ref->{'readID'}]} and $read_info{$line[$output_titles_ID_ref->{'readID'}]}[1] != $line[$output_titles_ID_ref->{'line_num'}];
			my $chr = $line[$output_titles_ID_ref->{'chr'}];
			my $orig_group_id = $line[$output_titles_ID_ref->{'group_ID'}];
			my $read_start = $line[$output_titles_ID_ref->{'start'}];
			my $read_end = $line[$output_titles_ID_ref->{'end'}];
			my @possible_groups = &lookup_final_group($num, $orig_group_id, $chr, $read_start, $read_end, $group_sep_in_all_ref, $group_all_info_ref, $chr_group_ref, $gtf_read_groups);
			if (scalar(@possible_groups) == 0) {
				$summary_data_ref->{'num_reads_filtered_no_assigned_read_group'} ++;
				next;
			}

			if ((!defined($pre_chr)) or ($chr ne $pre_chr)) {
				&print_read3_lines(\@pending_read3_lines, $read3_handle, $summary_data_ref);
				$pending_read3_lines_group = undef;

				%SJ_updated_all=();
				if (length($line[$output_titles_ID_ref->{'chr'}]) <= 5) {
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

			if ($line[$output_titles_ID_ref->{'SJcorSeqRef'}] ne 'NA') {
				my @SJcorSeq_strings = split ',', $line[$output_titles_ID_ref->{'SJcorSeqRef'}];
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
				
				$line[$output_titles_ID_ref->{'SJcorSeqRef'}] = join ',', @SJcorSeq_strings;
			}

			if ($line[$output_titles_ID_ref->{'readSeq'}] eq 'NA') {
				if (exists $read_info{$line[$output_titles_ID_ref->{'readID'}]} and defined $read_info{$line[$output_titles_ID_ref->{'readID'}]}[3]) {
					my $read_seq = $read_info{$line[$output_titles_ID_ref->{'readID'}]}[3];
					if ($read_seq eq 'short') {
						$summary_data_ref->{'num_reads_filtered_missing_full_sequence'} ++;
						next;
					}
					$line[$output_titles_ID_ref->{'readSeq'}] = $read_seq;
				} else {
					$line[$output_titles_ID_ref->{'readSeq'}] = 'read_not_recorded';
				}
			}

			$pending_read3_lines_group = &maybe_print_read3_lines($output_titles_ID_ref->{'group_ID'}, $chr, $read_start, $chr_group_ref, $group_all_info_ref, $pending_read3_lines_group, \@possible_groups, \@pending_read3_lines, \@line, $read3_handle, $summary_data_ref);
		}
		close READ;

		&print_read3_lines(\@pending_read3_lines, $read3_handle, $summary_data_ref);
		close $read3_handle;
		$pending_read3_lines_group = undef;

		if (!defined $keep_tmp){
			unlink "$out/$num/sam.list";
			unlink "$out/$num/sam.list2";
		}
	}

	sub parallel_scan1 {
		my ($file, $format, $input_ref, $summary_data_ref) = @{$_[0]};
		my $file_ID_ref = $input_ref->{'file_ID'};
		my $out = $input_ref->{'out'};
		my $mapq_cutoff = $input_ref->{'mapq_cutoff'};
		my $chrM = $input_ref->{'chrM'};
		my $cont_del_max = $input_ref->{'cont_del_max'};
		my $chr_seq_len_ref = $input_ref->{'chr_seq_len'};
		my $inserted_cont_cutoff = $input_ref->{'inserted_cont_cutoff'};
		my $group_extra_range = $input_ref->{'group_extra_range'};
		my $SJFS_dist = $input_ref->{'SJFS_dist'};
		my $chr_seq_ref = $input_ref->{'chr_seq'};
		my $sep_SJ = $input_ref->{'sep_SJ'};
		my $sam_sample_ref = $input_ref->{'sam_sample'};

		my $key_read;
		my $num = $file_ID_ref->{$file};
		my $sample = $sam_sample_ref->{$file};
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

		my (%SJ_read, %group_info, $pre_chr_end, @pre_group);
		my %processed_chrs = ();
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
				my $is_secondary = &ten2b($line[1], 9);
				if ($line[2] eq $chrM) {
					$summary_data_ref->{'number_of_alignments_filtered_for_chrM'} ++;
					next;
				}
				if ($is_secondary) {
					$summary_data_ref->{'number_of_alignments_filtered_for_secondary'} ++;
					next;
				}
				if ($line[4] < $mapq_cutoff) {
					$summary_data_ref->{'number_of_alignments_filtered_for_mapping_quality'} ++;
					next;
				}
				my $msidn = &MSIDN_border($line[5], $cont_del_max);

				$line[11] =~ /NM:i:(\d+)/;
				my $NM_num = $1;
				if (!exists $chr_seq_len_ref->{$line[2]}) {
					die "chr name ($line[2]) not recognized for $line[0] in $file";
				}
				if (${$msidn}{'max_I'} >= $inserted_cont_cutoff){
					$summary_data_ref->{'number_of_alignments_filtered_for_max_insertion'} ++;
					next;
				} elsif ($line[3] - 1 + ${$msidn}{'mappedGenomeN'} > $chr_seq_len_ref->{$line[2]}) {
					$summary_data_ref->{'number_of_alignments_filtered_for_exceeding_chr_coordinates'} ++;
					next;
				}

				my $end = $line[3]-1;
				$end += $_ for @{${$msidn}{'exonIntronRef'}};
				my $out;
				$out .= "$line[0]\t$sample\t${$msidn}{'read_length'}\t$line[1]\t$line[2]\t$line[3]\t$line[4]\t$end\t$notSameStrand\t${$msidn}{'mappedGenome'}\t";
				$out .= "$_," for @{${$msidn}{'clip_ends'}};
				$out .= "\t";
				my $is_chr_switch = 0;
				my $is_switch_to_old_chr = 0;
				my $is_lower_start_coord = 0;
				if (defined $pre_chr_end) {
					$is_chr_switch = $line[2] ne ${$pre_chr_end}[0];
					$is_switch_to_old_chr = ($is_chr_switch and exists $processed_chrs{$line[2]});
					$is_lower_start_coord = (!$is_chr_switch and $line[3] < ${$pre_chr_end}[1]);
				}

				if ($is_switch_to_old_chr or $is_lower_start_coord) {
					die "$file is not sorted";
				}
				$processed_chrs{$line[2]} = 1;

				if ( !defined($pre_chr_end) or $is_chr_switch or $line[3] > ${$pre_chr_end}[2]+$group_extra_range ){
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
					$pre_chr_end = [ $line[2], &max($line[3]-$group_extra_range, 0), &min($end+$group_extra_range, $chr_seq_len_ref->{$line[2]}) ];
				} elsif (${$pre_chr_end}[2] < &min($end+$group_extra_range, $chr_seq_len_ref->{$line[2]})) {
					${$pre_chr_end}[2] = &min($end+$group_extra_range, $chr_seq_len_ref->{$line[2]});
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
									$seq_ref = substr($chr_seq_ref->{$line[2]}, $exon_end2+$SJFS_dist, $substr_length);
									$seq_read = substr($line[9], ${$M_info}[3]+($SJFS_dist-${$M_info}[2]), $substr_length);
								} else {
									print DEFAULT "\tNA";
									next;
								}


							} else {
								if ($substr_length > 0) {
									$seq_ref = substr($chr_seq_ref->{$line[2]}, $exon_end2 + ${$M_info}[2], $substr_length);
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
								$seq_ref = substr($chr_seq_ref->{$line[2]}, $exon_end-$SJFS_dist, $SJFS_dist-${$M_info}[1]);
								$seq_read = substr($line[9], ${$M_info}[3]+${$M_info}[0]-($SJFS_dist-${$M_info}[1]), $SJFS_dist-${$M_info}[1]);
								$substr_current_num += &hamming_distance("\U$seq_read", "\U$seq_ref");
							} else {
								$seq_ref = substr($chr_seq_ref->{$line[2]}, $exon_end - ${$M_info}[1]-${$M_info}[0], ${$M_info}[0]);
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
								$seq_ref = substr($chr_seq_ref->{$line[2]}, $exon_end + ${$M_info}[2], $SJFS_dist-${$M_info}[2]);
								$seq_read = substr($line[9], ${$M_info}[3], $SJFS_dist-${$M_info}[2]);
								$substr_current_num += &hamming_distance("\U$seq_read", "\U$seq_ref");
							} else {
								$seq_ref = substr($chr_seq_ref->{$line[2]}, $exon_end + ${$M_info}[2], ${$M_info}[0]);
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
										print SJ "NA" for @perfect_reads; # TODO no for needed
									}
									
									print SJ "\t";
									print SJ "$_," for keys %{${$SJ_info_ref}[-1]};
									print SJ "\n";
								}
							}
						}

		close GROUP;
		close SJ;
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
		my @counts = split /[MSIDHN=X]/, $_[0];
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
			} elsif (($styles[$i] eq 'M') or ($styles[$i] eq '=') or ($styles[$i] eq 'X')) {
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
			} else {
				die "unexpected cigar operator $styles[$i]: $_[0]";
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

	sub load_annotation_file {
		my ($anno, $fa, $anno_coords_by_chr_by_gene_ref, $anno_SJ_ref, $summary_data_ref, $chr_seq_ref, $strand_digit_ref, $sep_SJ, $gtf_read_groups) = @_;
		my (%anno_exons, %isoform_info, %anno_chr_names);
		open(my $anno_handle, '<', $anno ) or die "cannot open $anno: $!";
		while(<$anno_handle>) {
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
			my $chr = $line[0];
			my $start = $line[3];  # gtf coords are 1-based
			my $end = $line[4];
			my $strand = $line[6];
			$isoform_info{$current_isoform} = [$current_gene, $strand, $chr];
			$anno_exons{$current_isoform}{$start - 1} = $end;
			$anno_chr_names{$chr} = 1;
			if ($gtf_read_groups) {
				if (exists $anno_coords_by_chr_by_gene_ref->{$chr}{$current_gene}) {
					my $gene_coords_ref = $anno_coords_by_chr_by_gene_ref->{$chr}{$current_gene};
					$gene_coords_ref->[0] = &min($gene_coords_ref->[0], $start);
					$gene_coords_ref->[1] = &max($gene_coords_ref->[1], $end);
				} else {
					$anno_coords_by_chr_by_gene_ref->{$chr}{$current_gene} = [$start, $end];
				}
			}
		}
		close $anno_handle;

		my $num_annotated_isoforms = scalar(keys %isoform_info);
		if ($num_annotated_isoforms == 0) {
			die "No isoforms found in $anno";
		}
		$summary_data_ref->{'num_annotated_isoforms'} = $num_annotated_isoforms;

		my $num_chrs_only_in_anno = 0;
		my $num_chrs_only_in_fa = 0;
		my $num_chrs_in_anno_and_fa = 0;
		for my $chr (keys %{$chr_seq_ref}) {
			if (exists $anno_chr_names{$chr}) {
				$num_chrs_in_anno_and_fa ++;
			} else {
				$num_chrs_only_in_fa ++;
			}
		}
		$num_chrs_only_in_anno = scalar(keys %anno_chr_names) - $num_chrs_in_anno_and_fa;

		if ($num_chrs_in_anno_and_fa == 0) {
			die "No overlap in chromosome names in $anno and $fa";
		}
		$summary_data_ref->{'num_chrs_only_in_anno'} = $num_chrs_only_in_anno;
		$summary_data_ref->{'num_chrs_only_in_fa'} = $num_chrs_only_in_fa;
		$summary_data_ref->{'num_chrs_in_anno_and_fa'} = $num_chrs_in_anno_and_fa;

		while (my ($isoform, $exon_start_ref) = each %anno_exons) {
			my @exon_start_sort = sort {$a <=> $b} keys %{$exon_start_ref};
			my $strand = $isoform_info{$isoform}[1];
			my $chr = $isoform_info{$isoform}[2];
			for my $i (1 .. $#exon_start_sort){
				my $prev_exon_end = ${$exon_start_ref}{$exon_start_sort[$i-1]};
				my $exon_start = $exon_start_sort[$i];
				my $SJ2_no_strand = $chr.$sep_SJ.$prev_exon_end.$sep_SJ.$exon_start;
				my $SJ2 = $SJ2_no_strand.$sep_SJ.$strand_digit_ref->{$strand};
				if (!exists $anno_SJ_ref->{$chr}{$SJ2}) {
					$anno_SJ_ref->{$chr}{$SJ2} = [$strand_digit_ref->{$strand},$prev_exon_end,$exon_start];
				}
			}
		}
	}

	sub load_sj_bed_file {
		my ($SJ_bed, $strand_digit_ref, $chr_seq_len_ref, $anno_SJ_ref, $sep_SJ) = @_;
		open(my $bed_handle, '<', $SJ_bed) or die "cannot open $SJ_bed: $!";
		while(<$bed_handle>) {
			chomp;
			my @line = split /\t/;
			my $chr = $line[0];
			my $sj_start = $line[1];
			my $sj_end = $line[2];

			if (!exists $strand_digit_ref->{$line[5]}) {
				die "Strand is not recognized for $_";
			}
			my $strand = $strand_digit_ref->{$line[5]};

			if (!exists $chr_seq_len_ref->{$chr}) {
				die "Chromosome is not found for $_";
			}
			if ($chr_seq_len_ref->{$chr} <= $sj_end) {
				die "Downstream splice site coordinate should be less than chromosome length($chr_seq_len_ref->{$chr} <= $sj_end): $_";
			}
			if ($sj_start <= 0) {
				die "Upstream splice site coordinate should be positive($sj_start <= 0): $_";
			}

			my $SJ2 = $chr.$sep_SJ.$sj_start.$sep_SJ.$sj_end.$sep_SJ.$strand;
			$anno_SJ_ref->{$chr}{$SJ2} = [$strand_digit_ref->{$line[5]},$sj_start,$sj_end] if !exists $anno_SJ_ref->{$chr}{$SJ2};
		}
		close $bed_handle;
	}

	sub create_coordinate_groups_for_genes {
		my ($anno_coords_by_gene_ref, $anno_gene_groups_ref, $gene_group_i_offset) = @_;
		my $gene_group_i = $gene_group_i_offset;
		my @gene_sort = sort {
			my $a_coords = $anno_coords_by_gene_ref->{$a};
			my $a_start = $a_coords->[0];
			my $a_end = $a_coords->[1];
			my $b_coords = $anno_coords_by_gene_ref->{$b};
			my $b_start = $b_coords->[0];
			my $b_end = $b_coords->[1];
			($a_start <=> $b_start) or ($a_end <=> $b_end);
		} keys %{$anno_coords_by_gene_ref};

		for my $gene_id (@gene_sort) {
			my $gene_coords_ref = $anno_coords_by_gene_ref->{$gene_id};
			my $gene_start = $gene_coords_ref->[0];
			my $gene_end = $gene_coords_ref->[1];
			if (!exists $anno_gene_groups_ref->{$gene_group_i}) {
				$anno_gene_groups_ref->{$gene_group_i} = [$gene_start, $gene_end];
			} else {
				my $old_group_coords_ref = $anno_gene_groups_ref->{$gene_group_i};
				if ($gene_start > $old_group_coords_ref->[1]) {
					$gene_group_i ++;
					$anno_gene_groups_ref->{$gene_group_i} = [$gene_start, $gene_end];
				} elsif ($gene_end > $old_group_coords_ref->[1]) {
					$anno_gene_groups_ref->{$gene_group_i}->[1] = $gene_end;
				}
			}
		}
	}

	sub finalize_read_groups {
		my ($gtf_read_groups, $group_sep_ref, $group_all_info_ref, $chr_group_ref, $group_sep_in_all_ref, $chr_sort_ref, $anno_coords_by_chr_by_gene_ref) = @_;
		my $gene_group_i = 0;
		for my $chr (@{$chr_sort_ref}) {
			my $orig_group_ref = $group_sep_ref->{$chr};
			if (keys %{$orig_group_ref} == 0) {
				next;
			}

			my @orig_group_sort = sort {
				my $a_coords = $orig_group_ref->{$a};
				my $a_start = $a_coords->[0];
				my $a_end = $a_coords->[1];
				my $b_coords = $orig_group_ref->{$b};
				my $b_start = $b_coords->[0];
				my $b_end = $b_coords->[1];
				($a_start <=> $b_start) or ($a_end <=> $b_end);
			} keys %{$orig_group_ref};

			if ($gtf_read_groups) {
				my $anno_coords_by_gene_ref = $anno_coords_by_chr_by_gene_ref->{$chr};
				# Only create groups if there are annotated genes for this $chr
				if (keys %{$anno_coords_by_gene_ref} == 0) {
					next;
				}
				my %anno_gene_groups;
				# The offset lets this chr start with the next available group_i.
				my $gene_group_i_offset = $gene_group_i;
				&create_coordinate_groups_for_genes($anno_coords_by_gene_ref, \%anno_gene_groups, $gene_group_i_offset);
				my $highest_coord_used_for_chr = $orig_group_ref->{$orig_group_sort[-1]}[1];

				my $orig_group_i = 0;
				my @sorted_gene_groups = sort {$a <=> $b} keys %anno_gene_groups;
				for my $sort_i (0 .. $#sorted_gene_groups) {
					$gene_group_i = $sorted_gene_groups[$sort_i];
					my $group_coords_ref = $anno_gene_groups{$gene_group_i};
					my $group_start = $group_coords_ref->[0];
					my $group_end = $group_coords_ref->[1];

					# The space between gene groups is added to each adjacent gene group.
					# The first gene group starts at 0.
					# The last gene group ends at $highest_coord_used_for_chr.
					my $is_last_gene_group_for_chr = $sort_i == $#sorted_gene_groups;
					if ($is_last_gene_group_for_chr) {
						$group_end = $highest_coord_used_for_chr;
					} else {
						my $next_gene_group_i = $sorted_gene_groups[$sort_i + 1];
						my $next_group_coords_ref = $anno_gene_groups{$next_gene_group_i};
						my $next_group_start = $next_group_coords_ref->[0];
						$group_end = $next_group_start - 1;
					}
					if ($sort_i == 0) {
						$group_start = 0;
					} else {
						my $prev_gene_group_i = $sorted_gene_groups[$sort_i - 1];
						my $prev_group_coords_ref = $anno_gene_groups{$prev_gene_group_i};
						my $prev_group_end = $prev_group_coords_ref->[1];
						$group_start = $prev_group_end + 1;
					}

					push @{$group_all_info_ref}, [$group_start, $group_end];
					push @{$chr_group_ref->{$chr}}, $gene_group_i;

					while ($orig_group_i < scalar(@orig_group_sort)) {
						my $orig_group_id = $orig_group_sort[$orig_group_i];
						my $orig_group_coord_ref = $orig_group_ref->{$orig_group_id};
						my $orig_group_start = $orig_group_coord_ref->[0];
						my $orig_group_end = $orig_group_coord_ref->[1];
						# Each original group is assigned to the gene-based group
						# that contains the original group start.
						if (($group_start <= $orig_group_start) and ($orig_group_start <= $group_end)) {
							$group_sep_in_all_ref->{$orig_group_id} = $gene_group_i;
							$orig_group_i ++;
						} else {
							last;
						}
					}
				}
				$gene_group_i ++;
			} else {
				my $first_group_id = $orig_group_sort[0];
				my $first_group_ref = ${$orig_group_ref}{$first_group_id};
				push @{$group_all_info_ref}, [$first_group_ref->[0], $first_group_ref->[1]];
				$group_sep_in_all_ref->{$first_group_id} = $gene_group_i;
				push @{$chr_group_ref->{$chr}}, $gene_group_i;
				for my $group_id (@orig_group_sort) {
					my $group_start = ${$orig_group_ref}{$group_id}[0];
					my $group_end = ${$orig_group_ref}{$group_id}[1];
					my $working_group_end = $group_all_info_ref->[-1][1];
					if ($group_start > $working_group_end) {
						$gene_group_i ++;
						push @{$group_all_info_ref}, [$group_start, $group_end];
						push @{$chr_group_ref->{$chr}}, $gene_group_i;
					} elsif ($group_end > $working_group_end) {
						$group_all_info_ref->[-1][1] = $group_end;
					}
					$group_sep_in_all_ref->{$group_id} = $gene_group_i;
				}
				$gene_group_i ++;
			}
		}

		my $total_groups = @{$group_all_info_ref};
		return $total_groups;
	}

	sub lookup_final_group {
		my ($file_num, $orig_group_id, $chr, $start, $end, $group_sep_in_all_ref, $group_all_info_ref, $chr_group_ref, $gtf_read_groups) = @_;

		my @possible_groups = ();
		my $lookup_key = $file_num.'_'.$orig_group_id;
		if (!exists $group_sep_in_all_ref->{$lookup_key}) {
			return @possible_groups;
		}

		my $group_id = $group_sep_in_all_ref->{$lookup_key};
		if (!$gtf_read_groups) {
			push @possible_groups, $group_id;
			return @possible_groups;
		}

		# $group_id corresponds to the gene-based group that contains the
		# original group start. Use $group_id as a starting point to search
		# @group_all_info for a gene-based group that contains [$start, $end].
		my $lower_limit = $group_id;
		my $upper_limit = $chr_group_ref->{$chr}[-1];
		while ($lower_limit <= $upper_limit) {
			my $group_coords_ref = $group_all_info_ref->[$group_id];
			my $group_start = $group_coords_ref->[0];
			my $group_end = $group_coords_ref->[1];
			my $any_check_failed = 0;
			if ($group_start > $start) {
				$any_check_failed = 1;
				$upper_limit = $group_id - 1;
			}
			if ($group_end < $end) {
				$any_check_failed = 1;
				$lower_limit = $group_id + 1;
			}
			if (!$any_check_failed) {
				push @possible_groups, $group_id;
				# [$start, $end] could be in the shared region between two groups.
				my @additional_checks = (($group_id - 1), ($group_id + 1));
				for $group_id (@additional_checks) {
					if ($lower_limit <= $group_id and $group_id <= $upper_limit) {
						$group_coords_ref = $group_all_info_ref->[$group_id];
						$group_start = $group_coords_ref->[0];
						$group_end = $group_coords_ref->[1];
						if ($group_start <= $start and $end <= $group_end) {
							push @possible_groups, $group_id;
						}
					}
				}
				last;
			}
			$group_id = int(($lower_limit + $upper_limit) / 2);
		}

		@possible_groups = sort {$a <=> $b} @possible_groups;
		return @possible_groups;
	}

	sub thread_work_loop {
		my ($input_queue, $output_queue) = @{$_[0]};

		# Allow the main thread to stop other threads with kill('TERM')
		$SIG{'TERM'} = sub { threads->exit(); };

		while (defined(my $work_details = $input_queue->dequeue())) {
			my %summary_data = ();
			&parallel_scan_prep($work_details, \%summary_data);
			my $serialized_summary :shared;
			$serialized_summary = freeze(\%summary_data);
			$output_queue->enqueue($serialized_summary);  # signal that work is done
		}
		$output_queue->end();
	}

	sub wait_for_threads_to_complete_work {
		my ($num_threads, $thread_input_queues_ref, $thread_output_queues_ref, $worker_threads_ref, $summary_data_ref) = @_;
		my @running_threads = 0 .. ($num_threads - 1);
		my @still_running_threads = ();
		my $first_check = 1;
		while (@running_threads) {
			if ($first_check) {
				$first_check = 0;
			} else {
				# Give the threads a chance to work
				my $sleep_seconds = 5;
				sleep $sleep_seconds;
			}
			for my $thread_i (@running_threads) {
				my $thread_is_joinable = $worker_threads_ref->[$thread_i]->is_joinable();
				my $thread_result = $thread_output_queues_ref->[$thread_i]->dequeue_nb();
				if (defined($thread_result)) {
					&merge_summary_results($summary_data_ref, $thread_result);
					print "Worker $thread_i finished reporting.\n";
				} elsif ($thread_is_joinable) {
					# The thread should not be joinable until the work queue is ended.
					# This thread must have had an error.
					&cleanup_threads_and_exit_if_error($thread_input_queues_ref, $worker_threads_ref);
				} else {
					push @still_running_threads, $thread_i;
				}
			}

			@running_threads = @still_running_threads;
			@still_running_threads = ();
		}
	}

	sub cleanup_threads_and_exit_if_error {
		my ($thread_input_queues_ref, $worker_threads_ref) = @_;
		my $sleep_seconds = 1;
		my $any_error = 0;
		for my $thread_i (0 .. $#{$worker_threads_ref}) {
			if ($worker_threads_ref->[$thread_i]->is_joinable()) {
				print "Worker $thread_i terminated early.\n";
				$any_error += 1;
			}
		}

		# After end() is called on an input queue the worker will see
		# 'undef' and then exit.
		for my $thread_i (0 .. $#{$worker_threads_ref}) {
			$thread_input_queues_ref->[$thread_i]->end();
		}
		# Give the threads a chance to exit
		sleep $sleep_seconds;

		for my $thread_i (0 .. $#{$worker_threads_ref}) {
			if (!$worker_threads_ref->[$thread_i]->is_joinable()) {
				print "Terminating worker $thread_i.\n";
				$worker_threads_ref->[$thread_i]->kill('TERM');
				$any_error += 1;
			}
		}
		# Give the threads a chance to exit
		if ($any_error) {
			sleep $sleep_seconds;
		}

		for my $thread_i (0 .. $#{$worker_threads_ref}) {
			if ($worker_threads_ref->[$thread_i]->is_joinable()) {
				$worker_threads_ref->[$thread_i]->join();
			} else {
				print "Worker $thread_i not responding.\n";
				$any_error += 1;
			}
		}
		if ($any_error) {
			die "Exiting due to error in worker thread\n";
		}
	}

	sub merge_summary_results {
		my ($summary_data_ref, $thread_result) = @_;
		my $thread_summary_ref = thaw($thread_result);
		while (my ($key, $value) = each %{$thread_summary_ref}) {
			$summary_data_ref->{$key} += $value;
		}
	}

	sub initialize_summary_data {
		my $summary_data_ref = $_[0];
		$summary_data_ref->{'num_chrs_only_in_anno'} = 0;
		$summary_data_ref->{'num_chrs_only_in_fa'} = 0;
		$summary_data_ref->{'num_chrs_in_anno_and_fa'} = 0;
		$summary_data_ref->{'num_annotated_isoforms'} = 0;
		$summary_data_ref->{'num_annotated_splice_junctions'} = 0;
		$summary_data_ref->{'num_high_confidence_splice_junctions'} = 0;
		$summary_data_ref->{'total_splice_junction_read_count'} = 0;
		$summary_data_ref->{'perfect_splice_junction_read_count'} = 0;
		$summary_data_ref->{'number_of_read_groups'} = 0;
		$summary_data_ref->{'number_of_reads_output'} = 0;
		$summary_data_ref->{'number_of_alignments_filtered_for_chrM'} = 0;
		$summary_data_ref->{'number_of_alignments_filtered_for_secondary'} = 0;
		$summary_data_ref->{'number_of_alignments_filtered_for_mapping_quality'} = 0;
		$summary_data_ref->{'number_of_alignments_filtered_for_max_insertion'} = 0;
		$summary_data_ref->{'number_of_alignments_filtered_for_exceeding_chr_coordinates'} = 0;
		$summary_data_ref->{'num_reads_filtered_missing_full_sequence'} = 0;
		$summary_data_ref->{'num_reads_filtered_no_assigned_read_group'} = 0;
	}

	sub write_summary_file {
		my ($summary_path, $summary_data_ref) = @_;
		open(my $summary_handle, '>', $summary_path) or die "cannot write $summary_path: $!";
		print $summary_handle "$summary_data_ref->{'perl_command'}\n";
		print $summary_handle "number of chromosomes only in input annotation: $summary_data_ref->{'num_chrs_only_in_anno'}\n";
		print $summary_handle "number of chromosomes only in input FASTA: $summary_data_ref->{'num_chrs_only_in_fa'}\n";
		print $summary_handle "number of chromosomes in both annotation and FASTA: $summary_data_ref->{'num_chrs_in_anno_and_fa'}\n";
		print $summary_handle "number of isoforms in input annotation: $summary_data_ref->{'num_annotated_isoforms'}\n";
		print $summary_handle "number of splice junctions in input annotation: $summary_data_ref->{'num_annotated_splice_junctions'}\n";
		print $summary_handle "number of high confidence splice junctions: $summary_data_ref->{'num_high_confidence_splice_junctions'}\n";
		print $summary_handle "total over all splice junctions of supporting reads: $summary_data_ref->{'total_splice_junction_read_count'}\n";
		print $summary_handle "total over all splice junctions of perfect reads: $summary_data_ref->{'perfect_splice_junction_read_count'}\n";
		print $summary_handle "number of read groups: $summary_data_ref->{'number_of_read_groups'}\n";
		print $summary_handle "number of reads in output: $summary_data_ref->{'number_of_reads_output'}\n";
		print $summary_handle "number of chrM alignments filtered: $summary_data_ref->{'number_of_alignments_filtered_for_chrM'}\n";
		print $summary_handle "number of secondary alignments filtered: $summary_data_ref->{'number_of_alignments_filtered_for_secondary'}\n";
		print $summary_handle "number of alignments filtered for mapping quality: $summary_data_ref->{'number_of_alignments_filtered_for_mapping_quality'}\n";
		print $summary_handle "number of alignments filtered for a long insertion: $summary_data_ref->{'number_of_alignments_filtered_for_max_insertion'}\n";
		print $summary_handle "number of alignments filtered for unrecognized coordinates: $summary_data_ref->{'number_of_alignments_filtered_for_exceeding_chr_coordinates'}\n";
		print $summary_handle "number of reads filtered for missing full sequence: $summary_data_ref->{'num_reads_filtered_missing_full_sequence'}\n";
		print $summary_handle "number of reads filtered for not being matched to a read group: $summary_data_ref->{'num_reads_filtered_no_assigned_read_group'}\n";
		close $summary_handle;
	}

	sub print_with_timestamp {
		my ($message) = @_;
		my $time_value = localtime;
		print "[$time_value] $message\n";
	}
}

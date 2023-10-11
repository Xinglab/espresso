use strict;
use threads;
use Thread::Queue;
use Getopt::Long;
use Storable qw(store retrieve);

use File::Basename qw(dirname);
use lib dirname(__FILE__);
use ESPRESSO_Version;
use ESPRESSO_Q_Thread;

my $version_number = ESPRESSO_Version::get_version_number();
my $source_code_commit = ESPRESSO_Version::get_source_code_commit();
my $version_string = "Q_$version_number";
my ($help, $in_dir, $out_dir, $anno, $anno_C, $tsv_compt, $num_thread, $list_samples, $read_num_cutoff, $read_ratio_cutoff, $max_iterate, $SJ_dist, $internal_boundary_limit, $raw, $target_col_index, $tmp_output, $allow_longer_terminal_exons);

my $arguments_before_parsing = "@ARGV";
Getopt::Long::GetOptions (

	'list_samples|L=s'			=>	\$list_samples,
	'out_dir|O=s'				=>	\$out_dir,
	'anno|A=s'					=>	\$anno,
	'anno_C|C=s'				=>	\$anno_C,
	'tsv_compt|V=s'				=>	\$tsv_compt,
	'num_thread|T=i'			=>	\$num_thread,

	'help|H!'					=>	\$help,

	'read_ratio_cutoff|R=f'		=>	\$read_ratio_cutoff,
	'read_num_cutoff|N=i'		=>	\$read_num_cutoff,
	'max_iterate|M=i'			=>	\$max_iterate,
	'SJ_dist|S=i'				=>	\$SJ_dist,
	'internal_boundary_limit=i'	=>	\$internal_boundary_limit,
	'allow_longer_terminal_exons'	=>	\$allow_longer_terminal_exons,
	'raw!'						=>	\$raw,
	'target_col_index=i'		=>	\$target_col_index,
	'tmp_output=s'				=>	\$tmp_output

);

if (defined($help)) {
	print "
Program:  ESPRESSO (Error Statistics PRomoted Evaluator of Splice Site Options)
Version:  $version_string
Contact:  Yuan Gao <gaoy\@email.chop.edu, gy.james\@163.com>

Usage:    perl ESPRESSO_Q.pl -L work_dir/samples.tsv.updated -A anno.gtf

Arguments:

    -L, --list_samples
          tsv list of multiple samples (each bam in a line with 1st column as
          sorted bam file, 2nd column as sample name in output, 3rd column as
          directory of ESPRESSO_C results; this list can be generated by
          ESPRESSO_S according to the initially provided tsv list; required)
    -A, --anno
          input annotation file in GTF format (optional)
    -O, --out_dir
          output directory (default: directory of -L)
    -V, --tsv_compt
          output tsv for compatible isoform(s) of each read (optional)
    -T --num_thread
          how many threads to use (default: 5)

    -H, --help
          show this help information

    -N, --read_num_cutoff
          min perfect read count for all splice junctions of novel isoform
          (default: 2)
    -R, --read_ratio_cutoff
          min perfect read ratio for all splice junctions of novel isoform
          (default: 0)
    -S, --SJ_dist
          max number of bases that an alignment endpoint can extend past the
          start or end of a matched isoform
          (default: 35)
    --internal_boundary_limit
          max number of bases that an alignment endpoint can extend into an
          intron of a matched isoform
          (default: 6)
    --allow_longer_terminal_exons
          allow an alignment to match an isoform even if the alignment endpoint
          extends more than --SJ_dist past the start or end

";
} elsif ( !defined($list_samples) ) {
	print  "The following parameter is required:\n";
	print  "\t--list_samples/-L";
	print  "\nPlease use the --help or -H option to get usage information.\n";
} else {
	my $target_dir;
	if (rindex($list_samples, "/")>=0){
		$target_dir = (substr( $list_samples, 0, rindex($list_samples, "/") ));
	} else {
		$target_dir = '.';
	}

	if (!defined $out_dir) {
		$out_dir = $target_dir;
	} elsif ( !-d $out_dir ) {
		warn " Output directory $out_dir does not exist. Make it by myself\n";
		system("mkdir $out_dir"); 
	}

	if (defined $tsv_compt) {
		if (index($tsv_compt, '/') == -1) {
			$tsv_compt = $out_dir.'/'.$tsv_compt;
		}
	}

	$num_thread = 5 if !defined $num_thread;
	$read_ratio_cutoff = 0 if !defined $read_ratio_cutoff;
	$read_num_cutoff = 2 if !defined $read_num_cutoff;
	$SJ_dist = 35 if !defined $SJ_dist;
	$internal_boundary_limit = 6 if !defined $internal_boundary_limit;
	$max_iterate = 100 if !defined $max_iterate;

	$target_col_index = 7 if !defined $target_col_index;

	my $THREAD_IDLE = 'idle';
	my $THREAD_RUNNING = 'running';

	my @worker_threads;
	my %threads_by_id;
	my @thread_input_queues;
	my @thread_output_queues;
	for my $thread_id (0 .. ($num_thread - 1)) {
		push @thread_input_queues, Thread::Queue->new();
		push @thread_output_queues, Thread::Queue->new();
		my $thread = threads->new({'context' => 'void'}, \&thread_work_loop, [$thread_input_queues[-1], $thread_output_queues[-1]]);
		push @worker_threads, $thread;
		$threads_by_id{$thread_id} = $THREAD_IDLE;
	}

	my $base_in_dir;

	my $dot_index = rindex($list_samples, '.');
	if ($dot_index == 0) {
		$in_dir = substr($list_samples, 1);
	} else {
		$in_dir = substr($list_samples, 0, $dot_index);
	}

	my @in_dir_parts = split /[\/.]/, $in_dir;
	if (length($in_dir_parts[-1]) > 1) {
		$base_in_dir = $in_dir_parts[-2];
	} else {
		$base_in_dir = 'samples';
	}

	my $keep_tmp = 1 if defined $tmp_output;

	# The summary file will be written as a final step
	my $summary_path = "$out_dir/espresso_q_summary.txt";
	my %summary_data = ();
	&initialize_summary_data(\%summary_data);
	$summary_data{'perl_command'} = "$^X " . __FILE__ . " $arguments_before_parsing";

	my %strand_digit = ('+' => 0, '-' => 1);
	my %strand_symbol = ('0' => '+', '1' => '-', 'unknown' => 'unknown');
	my $sep_SJ = ':';

	&print_with_timestamp("Loading annotation");

	my (%isoform_info, %anno_SJ, %anno_SS, %isoform_SJ, %isoform_exon);
	if(defined $anno) {
		my (%anno_exons);
		open ANNO, "<", $anno or die "cannot open $anno: $!";
		while(<ANNO>) {
			chomp;
			my @line = split /\t/;
			next if $line[2] ne 'exon';
			my ($current_isoform, $current_gene);
			my ($gene_name,$isoform_name) = ('NA','NA');
			if($line[-1] =~ /transcript_id \"(\S+)\"/) {
				$current_isoform = $1;
			} else {
				die "no transcript_id found in $_";
			}
			if($line[-1] =~ /gene_id \"(\S+)\"/) {
				$current_gene = $1;
			} else {
				die "no gene_id found in $_";
			}
			if($line[-1] =~ /gene_name \"(\S+)\"/) {
				$gene_name = $1;
			}
			if($line[-1] =~ /transcript_name \"(\S+)\"/) {
				$isoform_name = $1;
			}
			$isoform_info{$current_isoform} = [$current_gene, $line[6], $line[0], $gene_name, $isoform_name];
			$anno_exons{$current_isoform}{$line[3]-1} = $line[4];

			$isoform_exon{$current_isoform}{($line[3]-1).':'.$line[4]} = [$line[3]-1, $line[4], $line[0], $line[6]];

		}
		close ANNO;
		my $num_annotated_isoforms = scalar(keys %isoform_info);
		if ($num_annotated_isoforms == 0) {
			die "No isoforms found in $anno";
		}
		$summary_data{'num_annotated_isoforms'} = $num_annotated_isoforms;

		while (my ($isoform, $exon_start_ref) = each %anno_exons) {
			my @exon_start_sort = sort {$a <=> $b} keys %{$exon_start_ref};
			my $strand = $isoform_info{$isoform}[1];
			my $chr = $isoform_info{$isoform}[2];
			for my $i (1 .. $#exon_start_sort){
				$anno_SS{'0'}{$chr.':'.${$exon_start_ref}{$exon_start_sort[$i-1]}.':'.$strand_digit{$strand}} ++;
				$anno_SS{'1'}{$chr.':'.$exon_start_sort[$i].':'.$strand_digit{$strand}} ++;
				#my $SJ2_no_strand = $chr.$sep_SJ.${$exon_start_ref}{$exon_start_sort[$i-1]}.$sep_SJ.$exon_start_sort[$i];
				#my $SJ2 = $SJ2_no_strand.$sep_SJ.$strand_digit{$strand};
				my $SJ2 = $chr.$sep_SJ.${$exon_start_ref}{$exon_start_sort[$i-1]}.$sep_SJ.$exon_start_sort[$i];
				$isoform_SJ{$chr}{$isoform}{$SJ2} = [${$exon_start_ref}{$exon_start_sort[$i-1]}, $exon_start_sort[$i], $strand_digit{$strand}];
				if (!exists $anno_SJ{$chr}{$SJ2}) {
					$anno_SJ{$chr}{$SJ2} = [$strand_digit{$strand},${$exon_start_ref}{$exon_start_sort[$i-1]},$exon_start_sort[$i],[$isoform]];
				} else {
					push @{$anno_SJ{$chr}{$SJ2}[-1]}, $isoform;
				}
			}
		}
	}
	$summary_data{'num_annotated_splice_junctions'} += (scalar keys %{$anno_SJ{$_}}) for keys %anno_SJ;
	$summary_data{'num_annotated_splice_sites'} += (scalar keys %{$anno_SS{$_}}) for keys %anno_SS;

	my (%isoform_SJ_complete, %anno_SJ_complete);
	if (defined $anno_C) {
		my (%anno_exons);

		open ANNOC, "<", $anno_C or die "cannot open $anno_C: $!";
		while(<ANNOC>) {
			chomp;
			my @line = split /\t/;
			next if $line[2] ne 'exon';
			my ($current_isoform, $current_gene);
			my ($gene_name,$isoform_name) = ('NA','NA');
			if($line[-1] =~ /transcript_id \"(\S+)\"/) {
				$current_isoform = $1;
			} else {
				die "no transcript_id found in $_";
			}
			if($line[-1] =~ /gene_id \"(\S+)\"/) {
				$current_gene = $1;
			} else {
				die "no transcript_id found in $_";
			}
			if($line[-1] =~ /transcript_name \"(\S+)\"/) {
				$isoform_name = $1;
			}
			$isoform_info{$current_isoform} = [$current_gene, $line[6], $line[0], $gene_name, $isoform_name];
			$anno_exons{$current_isoform}{$line[3]-1} = $line[4];

		}
		while (my ($isoform, $exon_start_ref) = each %anno_exons) {
			my @exon_start_sort = sort {$a <=> $b} keys %{$exon_start_ref};
			my $strand = $isoform_info{$isoform}[1];
			my $chr = $isoform_info{$isoform}[2];
			for my $i (1 .. $#exon_start_sort){
				#my $SJ2_no_strand = $chr.$sep_SJ.${$exon_start_ref}{$exon_start_sort[$i-1]}.$sep_SJ.$exon_start_sort[$i];
				#my $SJ2 = $SJ2_no_strand.$sep_SJ.$strand_digit{$strand};
				my $SJ2 = $chr.$sep_SJ.${$exon_start_ref}{$exon_start_sort[$i-1]}.$sep_SJ.$exon_start_sort[$i];
				$isoform_SJ_complete{$chr}{$isoform}{$SJ2} = [${$exon_start_ref}{$exon_start_sort[$i-1]}, $exon_start_sort[$i], $strand_digit{$strand}];
				if (!exists $anno_SJ_complete{$chr}{$SJ2}) {
					$anno_SJ_complete{$chr}{$SJ2} = [$strand_digit{$strand},${$exon_start_ref}{$exon_start_sort[$i-1]},$exon_start_sort[$i],[$isoform]];
				} else {
					push @{$anno_SJ_complete{$chr}{$SJ2}[-1]}, $isoform;
				}
			}
		}
	}

	&print_with_timestamp("Summarizing annotated isoforms");

	my (%multi_exon_isoform_end,%single_exon_isoform_end);
	while (my ($chr, $isoform_ref) = each %isoform_SJ) {
		while (my ($isoform, $SJ_ref) = each %{$isoform_ref}) {
			my @exon_sort = sort {${$isoform_exon{$isoform}}{$a}[0] <=> ${$isoform_exon{$isoform}}{$b}[0]} keys %{$isoform_exon{$isoform}};
			# TODO The use of '0' and '1' in multi_exon_isoform_end{$isoform} could
			# conflict with an actual exon boundary at coordinate 0 or 1.
			# Need to carefully check usage of multi_exon_isoform_end before
			# updating code to use 'start' and 'end' instead of '0' and '1'
			$multi_exon_isoform_end{$isoform}{'0'} = $isoform_exon{$isoform}{$exon_sort[0]}[0];
			$multi_exon_isoform_end{$isoform}{'1'} = $isoform_exon{$isoform}{$exon_sort[-1]}[1];
			$multi_exon_isoform_end{$isoform}{${$isoform_exon{$isoform}{$_}}[1]} = ${$isoform_exon{$isoform}{$_}}[0] for @exon_sort;
			$multi_exon_isoform_end{$isoform}{${$isoform_exon{$isoform}{$_}}[0]} = ${$isoform_exon{$isoform}{$_}}[1] for @exon_sort;

		}
	}

	while (my ($isoform, $exon_ref) = each %isoform_exon) {
		if (!exists $multi_exon_isoform_end{$isoform}){
			while (my ($exon, $exon_info_ref) = each %{$exon_ref}){
				$single_exon_isoform_end{${$exon_info_ref}[2]}{$isoform}{'0'} = ${$exon_info_ref}[0];
				$single_exon_isoform_end{${$exon_info_ref}[2]}{$isoform}{'1'} = ${$exon_info_ref}[1];
				$single_exon_isoform_end{${$exon_info_ref}[2]}{$isoform}{'strand'} = ${$exon_info_ref}[3];
			}
		}
	}

	&print_with_timestamp("Loading corrected splice junctions and alignment information by ESPRESSO");

	my (@samples_sort, %file_sample, %input_file);
	my (%file_input, %sample_file);
	open LIST_S, "<", $list_samples or die "cannot open $list_samples: $!";
	while(<LIST_S>) {
		chomp;
		my @line = split /\t/;
		if ( $_ ne '' and (@line != 3 or $line[0] eq '' or $line[1] eq '' or $line[2] eq '') ){
			die "Please make sure $list_samples has three columns seperated by tab.";
		} elsif ($_ eq '') {
			next;
		}
		if ( (!exists $file_sample{$line[0]} or $file_sample{$line[0]} eq $line[1]) and (!exists $file_input{$line[0]} or $file_input{$line[0]} eq $line[2])) {
			$file_sample{$line[0]} = $line[1];
			$sample_file{$line[1]}{$line[0]} ++;
			$input_file{$line[2]}{$line[0]} ++;
			$file_input{$line[0]} = $line[2];
		} elsif (exists $file_sample{$line[0]}) {
			die "One fastq file should only correspond to one sample. \nThe fastq file \"$line[0]\" corresponds to more than one samples: \"$file_sample{$line[0]}\", \"$line[1]\".";
		} elsif (exists $file_input{$line[0]}) {
			die "One fastq file should only correspond to one input directory. \nThe fastq file \"$line[0]\" corresponds to more than one input directories: \"$file_input{$line[0]}\", \"$line[2]\".";
		}
		
	}
	@samples_sort = sort {$a <=> $b or $a cmp $b} keys %sample_file;

	opendir DIR, $out_dir or die "cannot opendir $out_dir: $!";

	for my $file (readdir DIR) {
		if (-f "$out_dir/$file" and (index($file, 'read_final.tmp')>0 or index($file, 'all_SJ.tmp')>0)) {
			unlink "$out_dir/$file";
		}
	}
	closedir DIR;

	my (%chr_tmp_f, %chr_tmp_sj); #, %input_num

	while (my ($input_dir, $bam_file_ref) = each %input_file) {
		unless ($input_dir =~ /(^\d+$)/) {
			die "unexpected format for $input_dir";
		}
		
		my %read_sample;
		
		while (my ($file, $count) = each %{$bam_file_ref}) {
			my $sample = $file_sample{$file};
			if ($file =~ /\.bam$/) {
				open IN, "samtools view -h $file |" or die "cannot open $file using samtools: $!";
			} else {
				open IN, $file or die "can’t open $file";
			}
			while(<IN>) {
				my ($read_ID, undef) = split /\t/;
				$read_sample{$read_ID} = $sample;
			}
			close IN;
		}
		
		my $valid_read_final_num = 0;
		opendir DIR, "$target_dir/$input_dir" or die "cannot opendir $target_dir/$input_dir: $!";
		for my $file (readdir DIR) {
			my $rindex_read_final = rindex ($file, "_read_final.txt");
			if ( $rindex_read_final > 0 and $rindex_read_final+15==length($file) ) { # 
				my $chr = substr($file, 0, $rindex_read_final);
				my $tmp_f = $out_dir."/$chr.read_final.tmp";
				$valid_read_final_num ++;
				$chr_tmp_f{$chr} ++;
				open TMP_F, ">>", $tmp_f or die "cannot write $tmp_f: $!";
				open IN_F,  "<", "$target_dir/$input_dir".'/'.$file or die "cannot open $target_dir/$input_dir/$file: $!";
				while (<IN_F>) {
					s/\r\n//;
					chomp;
					my @line = split /\t/;
					if ($line[1] eq 'group_ID') {
						if ( exists $read_sample{$line[0]} ) {
							print TMP_F "$line[0]\t$line[1]\t$line[2]\t$read_sample{$line[0]}\t$target_dir/$input_dir\t$input_dir\t$line[3]\n";
						} else {
							die "Don't know which sample $line[0] in $target_dir/$input_dir/$file is from.";
						}
					} else {
						print TMP_F "$_\n";
					}
				}
				close TMP_F;
				close IN_F;
			}
		}
		die "No valid read_final.list can be found in $target_dir/$input_dir.\n" if $valid_read_final_num == 0;
		
		my ($last_chr, $tmp_sj, @last_chr_info);
		if(-s "$target_dir/$input_dir".'/sj.list' > 0) {
			open IN_SJ,  "<", "$target_dir/$input_dir".'/sj.list' or die "cannot open $target_dir/$input_dir/sj.list: $!";
			while (<IN_SJ>) {
				s/\r\n//;
				chomp;
				my @line = split /\t/;
				if (!defined $tmp_sj or $last_chr ne $line[2]) {
					if (defined $tmp_sj) {
						print TMP_SJ "$input_dir\t$_\n" for @last_chr_info;
						@last_chr_info = ();
						close TMP_SJ;
					}
					$tmp_sj = $out_dir."/$line[2].all_SJ.tmp";
					$chr_tmp_sj{$line[2]} ++;
					open TMP_SJ, ">>", $tmp_sj or die "cannot write $tmp_sj: $!";
				}
				push @last_chr_info, $_;
				$last_chr = $line[2];
			}
			print TMP_SJ "$input_dir\t$_\n" for @last_chr_info;
			@last_chr_info = ();
			close TMP_SJ;
			close IN_SJ;
		} else {
			warn "No valid sj.list can be found in $target_dir/$input_dir.\n";
		}
		
	}

	&print_with_timestamp("Categorizing reads according to annotation");

	# Process each chr on a thread.
	# Only run $num_thread at a time.
	my %arguments_shared_for_all_threads;
	my $should_create_default_handle = defined $tmp_output;
	my $should_create_tsv_compt_handle = defined $tsv_compt;
	$arguments_shared_for_all_threads{'out_dir'} = $out_dir;
	$arguments_shared_for_all_threads{'anno_SS'} = \%anno_SS;
	$arguments_shared_for_all_threads{'multi_exon_isoform_end'} = \%multi_exon_isoform_end;
	$arguments_shared_for_all_threads{'samples_sort'} = \@samples_sort;
	$arguments_shared_for_all_threads{'isoform_info'} = \%isoform_info;
	$arguments_shared_for_all_threads{'strand_symbol'} = \%strand_symbol;
	$arguments_shared_for_all_threads{'should_create_default_handle'} = $should_create_default_handle;
	$arguments_shared_for_all_threads{'should_create_tsv_compt_handle'} = $should_create_tsv_compt_handle;
	$arguments_shared_for_all_threads{'keep_tmp'} = $keep_tmp;
	$arguments_shared_for_all_threads{'target_col_index'} = $target_col_index;
	$arguments_shared_for_all_threads{'SJ_dist'} = $SJ_dist;
	$arguments_shared_for_all_threads{'internal_boundary_limit'} = $internal_boundary_limit;
	$arguments_shared_for_all_threads{'allow_longer_terminal_exons'} = $allow_longer_terminal_exons;
	$arguments_shared_for_all_threads{'raw'} = $raw;
	$arguments_shared_for_all_threads{'read_num_cutoff'} = $read_num_cutoff;
	$arguments_shared_for_all_threads{'read_ratio_cutoff'} = $read_ratio_cutoff;
	$arguments_shared_for_all_threads{'max_iterate'} = $max_iterate;
	my $path_to_stored_shared_arguments = $out_dir."/thread_shared_arguments.tmp";
	store \%arguments_shared_for_all_threads, $path_to_stored_shared_arguments;
	%arguments_shared_for_all_threads = ();
	my @argument_temp_files = ($path_to_stored_shared_arguments);

	my @remaining_chrs = keys %chr_tmp_f;
	while ((scalar @remaining_chrs) > 0) {
		&cleanup_finished_threads(\%threads_by_id, \@worker_threads, \@thread_input_queues);

		my @thread_ids = keys %threads_by_id;
		for my $thread_id (@thread_ids) {
			if ((scalar @remaining_chrs) == 0) {
				last;
			}

			my $thread_state = $threads_by_id{$thread_id};
			if ($thread_state ne $THREAD_IDLE) {
				next;
			}
			my $chr = pop @remaining_chrs;
			my $path_to_stored_chr_arguments = $out_dir."/thread_${chr}_arguments.tmp";
			my %arguments_for_chr;
			$arguments_for_chr{'chr'} = $chr;
			$arguments_for_chr{'anno_SJ'} = $anno_SJ{$chr};
			$arguments_for_chr{'exists_chr_tmp_sj'} = exists $chr_tmp_sj{$chr};
			$arguments_for_chr{'single_exon_isoform_end'} = $single_exon_isoform_end{$chr};
			$arguments_for_chr{'isoform_SJ'} = $isoform_SJ{$chr};
			$arguments_for_chr{'anno_SJ_complete'} = $anno_SJ_complete{$chr};
			$arguments_for_chr{'isoform_SJ_complete'} = $isoform_SJ_complete{$chr};
			store \%arguments_for_chr, $path_to_stored_chr_arguments;
			%arguments_for_chr = ();
			push @argument_temp_files, $path_to_stored_chr_arguments;
			my %thread_input;
			$thread_input{'chr'} = $path_to_stored_chr_arguments;
			$thread_input{'shared'} = $path_to_stored_shared_arguments;
			$thread_input_queues[$thread_id]->enqueue(\%thread_input);
			$threads_by_id{$thread_id} = $THREAD_RUNNING;
			&print_with_timestamp("thread_id: $thread_id starting to process chr: $chr");
		}

		if ((scalar @remaining_chrs) > 0) {
			my $sleep_seconds = 10;
			sleep $sleep_seconds;
		}
	}

	# Wait for any remaining threads
	while (1) {
		my $any_running = &cleanup_finished_threads(\%threads_by_id, \@worker_threads, \@thread_input_queues);
		if ($any_running == 1) {
			my $sleep_seconds = 10;
			sleep $sleep_seconds;
		} else {
			last;
		}
	}

	if (!$keep_tmp) {
		for my $argument_temp_file (@argument_temp_files) {
			unlink($argument_temp_file);
		}
	}

	my $gtf_path = $out_dir."/${base_in_dir}_N${read_num_cutoff}_R${read_ratio_cutoff}_updated.gtf";
	open(my $gtf_handle,  ">", $gtf_path) or die "cannot write $gtf_path: $!";
	print $gtf_handle "# ESPRESSO version: $version_number, source code commit: $source_code_commit\n";

	my $abu_path = $out_dir."/${base_in_dir}_N${read_num_cutoff}_R${read_ratio_cutoff}_abundance.esp";
	open(my $abu_handle,  ">", $abu_path) or die "cannot write $abu_path: $!";
	print $abu_handle "transcript_ID\ttranscript_name\tgene_ID";
	print $abu_handle "\t$_" for @samples_sort;
	print $abu_handle "\n";

	open(my $tsv_compt_handle, ">", $tsv_compt) or die "cannot write $tsv_compt: $!" if defined $tsv_compt;

	open(my $default_handle,  ">", $tmp_output) or die "cannot write $tmp_output: $!" if defined $tmp_output;

	for my $chr (keys %chr_tmp_f) {
		my $chr_gtf_path = &ESPRESSO_Q_Thread::gtf_path_for_chr($chr, $out_dir);
		&append_lines_from_path_to_handle($chr_gtf_path, $gtf_handle);
		unlink($chr_gtf_path) if !defined $keep_tmp;
		my $chr_abu_path = &ESPRESSO_Q_Thread::abu_path_for_chr($chr, $out_dir);
		&append_lines_from_path_to_handle($chr_abu_path, $abu_handle);
		unlink($chr_abu_path) if !defined $keep_tmp;
		if (defined $tmp_output) {
			my $chr_default_path = &ESPRESSO_Q_Thread::default_path_for_chr($chr, $out_dir);
			&append_lines_from_path_to_handle($chr_default_path, $default_handle);
			unlink($chr_default_path) if !defined $keep_tmp;
		}
		if (defined $tsv_compt) {
			my $chr_tsv_compt_path = &ESPRESSO_Q_Thread::tsv_compt_path_for_chr($chr, $out_dir);
			&append_lines_from_path_to_handle($chr_tsv_compt_path, $tsv_compt_handle);
			unlink($chr_tsv_compt_path) if !defined $keep_tmp;
		}
		my $chr_summary_path = &ESPRESSO_Q_Thread::summary_path_for_chr($chr, $out_dir);
		&merge_summary_results(\%summary_data, $chr_summary_path);
		unlink($chr_summary_path) if !defined $keep_tmp;
	}

	close $gtf_handle;
	close $abu_handle;
	close $tsv_compt_handle if defined $tsv_compt;
	close $default_handle if defined $tmp_output;

	&print_with_timestamp("cleaning up threads");
	&cleanup_all_threads_and_exit_if_error(\@worker_threads, \@thread_input_queues);
	&write_summary_file($summary_path, \%summary_data);
	&print_with_timestamp("ESPRESSO finished quantification");

	sub cleanup_finished_threads {
		my ($threads_by_id_ref, $worker_threads_ref, $thread_input_queues_ref) = @_;
		my $any_running = 0;
		my @thread_ids = keys %{$threads_by_id_ref};
		for my $thread_id (@thread_ids) {
			my $current_state = $threads_by_id_ref->{$thread_id};
			if ($current_state ne $THREAD_RUNNING) {
				next;
			}
			my $is_joinable = $worker_threads_ref->[$thread_id]->is_joinable();
			my $is_done = defined $thread_output_queues[$thread_id]->dequeue_nb();
			if ($is_done) {
				$threads_by_id_ref->{$thread_id} = $THREAD_IDLE;
				&print_with_timestamp("thread_id: $thread_id finished");
			} elsif ($is_joinable) {
				# The thread should not be joinable until the work queue is ended.
				# This thread must have had an error.
				&cleanup_all_threads_and_exit_if_error($worker_threads_ref, $thread_input_queues_ref);
			} else {
				$any_running = 1;
			}
		}
		return $any_running;
	}

	sub cleanup_all_threads_and_exit_if_error {
		my ($worker_threads_ref, $thread_input_queues_ref) = @_;
		my $num_threads = scalar(@{$worker_threads_ref});
		my $sleep_seconds = 1;
		my $any_error = 0;
		for my $thread_i (0 .. ($num_threads - 1)) {
			if ($worker_threads_ref->[$thread_i]->is_joinable()) {
					&print_with_timestamp("Worker $thread_i terminated early.");
					$any_error += 1;
			}
		}

		# After end() is called on an input queue the worker will see
		# 'undef' and then exit.
		for my $thread_i (0 .. ($num_threads - 1)) {
			$thread_input_queues_ref->[$thread_i]->end();
		}
		# Give the threads a chance to exit
		sleep $sleep_seconds;

		for my $thread_i (0 .. ($num_threads - 1)) {
			if (!$worker_threads_ref->[$thread_i]->is_joinable()) {
				&print_with_timestamp("Terminating worker $thread_i.");
				$worker_threads_ref->[$thread_i]->kill('TERM');
				$any_error += 1;
			}
		}
		# Give the threads a chance to exit
		if ($any_error) {
			sleep $sleep_seconds;
		}

		for my $thread_i (0 .. ($num_threads - 1)) {
			if ($worker_threads_ref->[$thread_i]->is_joinable()) {
				$worker_threads_ref->[$thread_i]->join();
			} else {
				&print_with_timestamp("Worker $thread_i not responding.");
				$any_error += 1;
			}
		}
		if ($any_error) {
			die "Exiting due to error in worker thread\n";
		}
	}

	# It seems that there may be a memory leak in ESPRESSO_Q_Thread.
	# ESPRESSO_Q_Thread is run as a separate process with system().
	# That ensures that memory is properly cleaned up by the OS when
	# the process exists.
	sub thread_work_loop {
		my ($input_queue, $output_queue) = @{$_[0]};

		# Allow the main thread to stop other threads with kill('TERM')
		$SIG{'TERM'} = sub { threads->exit(); };

		my $path_to_perl_executable = $^X;
		my $path_to_code_file = dirname(__FILE__)."/ESPRESSO_Q_Thread.pm";
		while (defined(my $work_details = $input_queue->dequeue())) {
			my $path_of_shared_arguments = $work_details->{'shared'};
			my $path_of_chr_arguments = $work_details->{'chr'};
			my $command = "$path_to_perl_executable $path_to_code_file $path_of_shared_arguments $path_of_chr_arguments";
			print "running command: $command\n";
			my $return_code = system($command);
			print "finished command: $command\n";
			if ($return_code != 0) {
				die "call for $command exited with return_code: $return_code\n";
			}
			$output_queue->enqueue(1);  # signal that work is done
		}
		$output_queue->end();
	}

	sub append_lines_from_path_to_handle {
		my ($read_path, $write_handle) = @_;
		open(my $read_handle, "<", $read_path) or die "cannot read $read_path: $!";
		while (<$read_handle>) {
			print $write_handle $_;
		}
		close $read_handle;
	}

	sub print_with_timestamp {
		my ($message) = @_;
		my $time_value = localtime;
		print "[$time_value] $message\n";
	}

	sub merge_summary_results {
		my ($summary_data_ref, $thread_result_path) = @_;
		my $thread_summary_ref = retrieve($thread_result_path);
		while (my ($key, $value) = each %{$thread_summary_ref}) {
			$summary_data_ref->{$key} += $value;
		}
	}

	sub initialize_summary_data {
		my $summary_data_ref = $_[0];
		$summary_data_ref->{'num_annotated_isoforms'} = 0;
		$summary_data_ref->{'num_annotated_splice_junctions'} = 0;
		$summary_data_ref->{'num_annotated_splice_sites'} = 0;
		$summary_data_ref->{'num_reads_assigned'} = 0;
		$summary_data_ref->{'num_reads_not_assigned'} = 0;
		$summary_data_ref->{'num_fsm_reads'} = 0;
		$summary_data_ref->{'num_fsm_perfect_match_reads'} = 0;
		$summary_data_ref->{'num_ism_reads'} = 0;
		$summary_data_ref->{'num_nic_reads'} = 0;
		$summary_data_ref->{'num_nnc_reads'} = 0;
		$summary_data_ref->{'num_single_exon_reads'} = 0;
		$summary_data_ref->{'num_reads_with_failed_SJ'} = 0;
		$summary_data_ref->{'num_fsm_chains'} = 0;
		$summary_data_ref->{'num_ism_chains'} = 0;
		$summary_data_ref->{'num_nic_chains'} = 0;
		$summary_data_ref->{'num_validated_nic_chains'} = 0;
		$summary_data_ref->{'num_nnc_chains'} = 0;
		$summary_data_ref->{'num_validated_nnc_chains'} = 0;
		$summary_data_ref->{'num_chains_with_failed_SJ'} = 0;
		$summary_data_ref->{'total_fsm_abundance'} = 0;
		$summary_data_ref->{'total_novel_ism_abundance'} = 0;
		$summary_data_ref->{'total_nic_abundance'} = 0;
		$summary_data_ref->{'total_nnc_abundance'} = 0;
		$summary_data_ref->{'total_single_exon_abundance'} = 0;
		$summary_data_ref->{'num_fsm_isoforms_detected'} = 0;
		$summary_data_ref->{'num_novel_ism_isoforms_detected'} = 0;
		$summary_data_ref->{'num_nic_isoforms_detected'} = 0;
		$summary_data_ref->{'num_nnc_isoforms_detected'} = 0;
		$summary_data_ref->{'num_single_exon_isoforms_detected'} = 0;
		$summary_data_ref->{'num_internal_exon_boundary_check_fails'} = 0;
		$summary_data_ref->{'num_terminal_exon_boundary_check_fails'} = 0;
	}

	sub write_summary_file {
		my ($summary_path, $summary_data_ref) = @_;
		open(my $summary_handle, '>', $summary_path) or die "cannot write $summary_path: $!";
		print $summary_handle "$summary_data_ref->{'perl_command'}\n";
		print $summary_handle "number of isoforms in input annotation: $summary_data_ref->{'num_annotated_isoforms'}\n";
		print $summary_handle "number of splice junctions in input annotation: $summary_data_ref->{'num_annotated_splice_junctions'}\n";
		print $summary_handle "number of splice sites in input annotation: $summary_data_ref->{'num_annotated_splice_sites'}\n";
		print $summary_handle "number of reads assigned: $summary_data_ref->{'num_reads_assigned'}\n";
		print $summary_handle "number of reads not assigned: $summary_data_ref->{'num_reads_not_assigned'}\n";
		print $summary_handle "number of full splice match reads: $summary_data_ref->{'num_fsm_reads'}\n";
		print $summary_handle "number of FSM reads with 'perfect match' endpoints: $summary_data_ref->{'num_fsm_perfect_match_reads'}\n";
		print $summary_handle "number of incomplete splice match reads: $summary_data_ref->{'num_ism_reads'}\n";
		print $summary_handle "number of novel in catalog reads: $summary_data_ref->{'num_nic_reads'}\n";
		print $summary_handle "number of novel not in catalog reads: $summary_data_ref->{'num_nnc_reads'}\n";
		print $summary_handle "number of single exon reads: $summary_data_ref->{'num_single_exon_reads'}\n";
		print $summary_handle "number of reads with a failed splice junction: $summary_data_ref->{'num_reads_with_failed_SJ'}\n";
		print $summary_handle "number of FSM splice junction chains: $summary_data_ref->{'num_fsm_chains'}\n";
		print $summary_handle "number of ISM splice junction chains: $summary_data_ref->{'num_ism_chains'}\n";
		print $summary_handle "number of NIC splice junction chains: $summary_data_ref->{'num_nic_chains'}\n";
		print $summary_handle "number of validated NIC chains: $summary_data_ref->{'num_validated_nic_chains'}\n";
		print $summary_handle "number of NNC splice junction chains: $summary_data_ref->{'num_nnc_chains'}\n";
		print $summary_handle "number of validated NNC chains: $summary_data_ref->{'num_validated_nnc_chains'}\n";
		print $summary_handle "number of splice junction chains with a failed junction: $summary_data_ref->{'num_chains_with_failed_SJ'}\n";
		print $summary_handle "total FSM abundance: $summary_data_ref->{'total_fsm_abundance'}\n";
		print $summary_handle "total novel ISM abundance: $summary_data_ref->{'total_novel_ism_abundance'}\n";
		print $summary_handle "total NIC abundance: $summary_data_ref->{'total_nic_abundance'}\n";
		print $summary_handle "total NNC abundance: $summary_data_ref->{'total_nnc_abundance'}\n";
		print $summary_handle "total single exon abundance: $summary_data_ref->{'total_single_exon_abundance'}\n";
		print $summary_handle "number of detected FSM isoforms: $summary_data_ref->{'num_fsm_isoforms_detected'}\n";
		print $summary_handle "number of detected novel ISM isoforms: $summary_data_ref->{'num_novel_ism_isoforms_detected'}\n";
		print $summary_handle "number of detected NIC isoforms: $summary_data_ref->{'num_nic_isoforms_detected'}\n";
		print $summary_handle "number of detected NNC isoforms: $summary_data_ref->{'num_nnc_isoforms_detected'}\n";
		print $summary_handle "number of detected single exon isoforms: $summary_data_ref->{'num_single_exon_isoforms_detected'}\n";
		print $summary_handle "number of internal exon boundary check failures: $summary_data_ref->{'num_internal_exon_boundary_check_fails'}\n";
		print $summary_handle "number of terminal exon boundary check failures: $summary_data_ref->{'num_terminal_exon_boundary_check_fails'}\n";
		close $summary_handle;
	}
}

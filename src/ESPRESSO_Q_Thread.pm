package ESPRESSO_Q_Thread;

use strict;

use Storable qw(retrieve);

if ((scalar @ARGV) == 2) {
	# Run as a script
	my $stored_shared_arguments = $ARGV[0];
	my $stored_chr_arguments = $ARGV[1];
	&main($stored_shared_arguments, $stored_chr_arguments);
} else {
	# Just for imports
}

sub main {
	my ($stored_shared_arguments, $stored_chr_arguments) = @_;
	my $shared_arguments_ref = retrieve($stored_shared_arguments);
	my $chr_arguments_ref = retrieve($stored_chr_arguments);

	my $out_dir = $shared_arguments_ref->{'out_dir'};
	my $anno_SS_ref = $shared_arguments_ref->{'anno_SS'};
	my $multi_exon_isoform_end_ref = $shared_arguments_ref->{'multi_exon_isoform_end'};
	my $samples_sort_ref = $shared_arguments_ref->{'samples_sort'};
	my $isoform_info_ref = $shared_arguments_ref->{'isoform_info'};
	my $strand_symbol_ref = $shared_arguments_ref->{'strand_symbol'};
	my $should_create_default_handle = $shared_arguments_ref->{'should_create_default_handle'};
	my $should_create_tsv_compt_handle = $shared_arguments_ref->{'should_create_tsv_compt_handle'};
	my $keep_tmp = $shared_arguments_ref->{'keep_tmp'};
	my $target_col_index = $shared_arguments_ref->{'target_col_index'};
	my $SJ_dist = $shared_arguments_ref->{'SJ_dist'};
	my $raw = $shared_arguments_ref->{'raw'};
	my $read_num_cutoff = $shared_arguments_ref->{'read_num_cutoff'};
	my $read_ratio_cutoff = $shared_arguments_ref->{'read_ratio_cutoff'};
	my $max_iterate = $shared_arguments_ref->{'max_iterate'};

	my $chr = $chr_arguments_ref->{'chr'};
	my $anno_SJ_ref = $chr_arguments_ref->{'anno_SJ'};
	my $exists_chr_tmp_sj = $chr_arguments_ref->{'exists_chr_tmp_sj'};
	my $single_exon_isoform_end_ref = $chr_arguments_ref->{'single_exon_isoform_end'};
	my $isoform_SJ_ref = $chr_arguments_ref->{'isoform_SJ'};
	my $anno_SJ_complete_ref = $chr_arguments_ref->{'anno_SJ_complete'};
	my $isoform_SJ_complete_ref = $chr_arguments_ref->{'isoform_SJ_complete'};

	&process_results_for_chr($chr, $out_dir, $anno_SJ_ref, $exists_chr_tmp_sj, $single_exon_isoform_end_ref, $anno_SS_ref, $isoform_SJ_ref, $multi_exon_isoform_end_ref, $anno_SJ_complete_ref, $isoform_SJ_complete_ref, $samples_sort_ref, $isoform_info_ref, $strand_symbol_ref, $should_create_default_handle, $should_create_tsv_compt_handle, $keep_tmp, $target_col_index, $SJ_dist, $raw, $read_num_cutoff, $read_ratio_cutoff, $max_iterate);
}

sub process_results_for_chr {
	my($chr, $out_dir, $anno_SJ_ref, $exists_chr_tmp_sj, $single_exon_isoform_end_ref, $anno_SS_ref, $isoform_SJ_ref, $multi_exon_isoform_end_ref, $anno_SJ_complete_ref, $isoform_SJ_complete_ref, $samples_sort_ref, $isoform_info_ref, $strand_symbol_ref, $should_create_default_handle, $should_create_tsv_compt_handle, $keep_tmp, $target_col_index, $SJ_dist, $raw, $read_num_cutoff, $read_ratio_cutoff, $max_iterate) = @_;
	my %chr_type;
	my $tmp_f = $out_dir."/$chr.read_final.tmp";
	my $tmp_sj = $out_dir."/$chr.all_SJ.tmp";
	my $default_path = &default_path_for_chr($chr, $out_dir);
	open(my $default_handle, ">", $default_path) or die "cannot write $default_path: $!" if $should_create_default_handle;
	my $tsv_compt_path = &tsv_compt_path_for_chr($chr, $out_dir);
	open(my $tsv_compt_handle, ">", $tsv_compt_path) or die "cannot write $tsv_compt_path: $!" if $should_create_tsv_compt_handle;
	my $abu_path = &abu_path_for_chr($chr, $out_dir);
	open(my $abu_handle, ">", $abu_path) or die "cannot write $abu_path: $!";
	my $gtf_path = &gtf_path_for_chr($chr, $out_dir);
	open(my $gtf_handle, ">", $gtf_path) or die "cannot write $gtf_path: $!";

	my %SJ_read_improved;
	my %read_filtered;
	my %read_info;
	my %SJ_read_order;
	my %unanno_SJ;

	my %all_SJ_reads;
	my ($read_input_ID,$strand_isoform);
	open(my $tmp_f_handle, "<", $tmp_f) or die "cannot open $tmp_f: $!";
	while (<$tmp_f_handle>) {
		s/\r\n//;
		chomp;
		my @line = split /\t/;
		if ($line[1] eq 'SJ'){

			if ($line[$target_col_index] ne 'NA') {
				my ($chr, $start, $end, $strand) = split ':', $line[$target_col_index];
				$line[$target_col_index] = join ':', ($chr, $start, $end);
				$SJ_read_improved{$read_input_ID}{$line[$target_col_index]} = [$chr, $start, $end, $strand_isoform, 0];
				if (exists $anno_SJ_ref->{$line[$target_col_index]}){
					$SJ_read_improved{$read_input_ID}{$line[$target_col_index]}[-1] = 1;
					$all_SJ_reads{$line[$target_col_index]}{$line[0]} = 0;
				} else {
					$unanno_SJ{$line[$target_col_index]} = [$strand, $start, $end];
					$all_SJ_reads{$line[$target_col_index]}{$line[0]} = 0;
				}
				$SJ_read_order{$read_input_ID}{$line[2]} = $line[$target_col_index];
			} else {
				$read_filtered{$read_input_ID} ++ if $target_col_index == 7;
				$SJ_read_order{$read_input_ID}{$line[2]} = 'x' if $target_col_index == 7;
			}
		} elsif ($line[1] eq 'start' or $line[1] eq 'end') {
			$read_info{$read_input_ID}{$line[1]} = $line[$target_col_index];
		} elsif ($line[1] eq 'group_ID') {
			$read_input_ID = $line[5].'_'.$line[6];
			$read_info{$read_input_ID}{'sample'} = $line[3];
			$read_info{$read_input_ID}{'ori_ID'} = $line[0];
		} elsif ($line[1] eq 'strand_isoform') {
			$strand_isoform = $line[2];
			$read_info{$read_input_ID}{'strand_isoform'} = $strand_isoform;
		}
	}
	close $tmp_f_handle;

	if ($exists_chr_tmp_sj) {
		open(my $tmp_sj_handle,  "<", $tmp_sj) or die "cannot open $tmp_sj: $!";
		while (<$tmp_sj_handle>) {
			s/\r\n//;
			chomp;
			my @line = split /\t/;
			next if $line[-2] eq 'NA';
			my @perfect_reads = split ',', $line[-2];
			if (exists $all_SJ_reads{$line[2]}) {
				$all_SJ_reads{$line[2]}{$_} ++ for @perfect_reads;
			}
		}
		close $tmp_sj_handle;
	}

	unlink($tmp_f) if !defined $keep_tmp;
	unlink($tmp_sj) if !defined $keep_tmp;
	
	my %group_reads;
	my %group_ends;
	
	my @group_sort;
#need update

	while (my ($read, $info_ref) = each %read_info){
		my ($sample_ID, $group_ID, $read_ID) = split '_', $read;
		push @{$group_reads{$group_ID}}, $read;
		if(!exists $group_ends{$group_ID}) {
			$group_ends{$group_ID} = [${$info_ref}{'start'}, ${$info_ref}{'end'}];
		} else {
			if (${$info_ref}{'start'} < $group_ends{$group_ID}[0]) {
				$group_ends{$group_ID}[0] = ${$info_ref}{'start'};
			} if (${$info_ref}{'end'} > $group_ends{$group_ID}[1]) {
				$group_ends{$group_ID}[1] = ${$info_ref}{'end'};
			}
		}
	}
	# The groups originally had non-overlapping coordinates, but
	# since ESPRESSO_C can correct the end points of each read,
	# the coordinates may now overlap.
	# Sort by both the start coordinate and end coordinate.
	@group_sort = sort {
		my $a_info = $group_ends{$a};
		my $a_start = $a_info->[0];
		my $a_end = $a_info->[1];
		my $b_info = $group_ends{$b};
		my $b_start = $b_info->[0];
		my $b_end = $b_info->[1];
		($a_start <=> $b_start) or ($a_end <=> $b_end);
	} keys %group_reads;

	# The single exon reads in a group will be checked against the
	# single exon isoforms.
	# Sort the isoforms and determine which isoforms overlap with
	# each group. An extra distance of $SJ_dist is allowed
	my @single_exon_isoforms_chr_sort = sort {
		my $a_info = $single_exon_isoform_end_ref->{$a};
		my $a_start = $a_info->{'0'};
		my $a_end = $a_info->{'1'};
		my $b_info = $single_exon_isoform_end_ref->{$b};
		my $b_start = $b_info->{'0'};
		my $b_end = $b_info->{'1'};
		($a_start <=> $b_start) or ($a_end <=> $b_end);
	} keys %{$single_exon_isoform_end_ref};

	my %single_exon_isoforms_by_group;
	{
		my $extra_length = $SJ_dist;
		# Since the isoforms and groups are sorted by coordinate, keep track
		# of which groups end before the remaining isoforms start
		my $min_group_i = 0;
		for my $isoform (@single_exon_isoforms_chr_sort) {
			my $isoform_info = $single_exon_isoform_end_ref->{$isoform};
			my $isoform_start = $isoform_info->{'0'};
			my $isoform_end = $isoform_info->{'1'};
			while ($min_group_i < (scalar @group_sort)) {
				my $min_group_id = $group_sort[$min_group_i];
				my $min_group_info = $group_ends{$min_group_id};
				my $min_group_end = $min_group_info->[1] + $extra_length;
				if ($isoform_start > $min_group_end) {
					# all remaninig isoforms have a start past this
					# group end so stop checking this group
					$min_group_i += 1;
					next;
				}
				last;
			}
			for my $group_i ($min_group_i .. $#group_sort) {
				my $group_id = $group_sort[$group_i];
				my $group_info = $group_ends{$group_id};
				my $group_start = $group_info->[0] - $extra_length;
				my $group_end = $group_info->[1] + $extra_length;
				if ($isoform_start <= $group_end and $isoform_end >= $group_start) {
					push @{$single_exon_isoforms_by_group{$group_id}}, $isoform;
					next;
				}
				if ($isoform_end < $group_start) {
					# all remaining groups have a start past this
					# isoform end so stop checking this isoform
					last;
				}
			}
		}
	}

	my @annotated_SJ_chr = sort {$anno_SJ_ref->{$a}[1] <=> $anno_SJ_ref->{$b}[1] or $anno_SJ_ref->{$a}[2] <=> $anno_SJ_ref->{$b}[2]} keys %{$anno_SJ_ref};
	my @unannotated_SJ_chr = sort {$unanno_SJ{$a}[1] <=> $unanno_SJ{$b}[1] or $unanno_SJ{$a}[2] <=> $unanno_SJ{$b}[2]} keys %unanno_SJ;

	my ($start_SJ_index_anno, $start_SJ_index_unanno) = (0,0);
	my $chr_SJ_reads_ref = \%all_SJ_reads;
	my %single_exon_isoform_total;
	
	for my $n (@group_sort) {

		my (%fsm, %ism, %nic, %nnc, %isoform_perfect_match);
		my @single_exon_reads;
		for my $read (@{$group_reads{$n}}) {
			if (exists $read_filtered{$read}){
				print $default_handle "$read\t$read_info{$read}{'ori_ID'}\tNA\tfiltered\tNA\n" if $should_create_default_handle; # and exists $SJ_read_improved{$read};
			} elsif (exists $SJ_read_improved{$read}){
				
				my @SJ_correct = grep {$SJ_read_improved{$read}{$_}[-1] == 1} keys %{$SJ_read_improved{$read}};
				if ( scalar(keys %{$SJ_read_improved{$read}}) != scalar(@SJ_correct) ){
					my @same_SSs = grep {exists $anno_SS_ref->{'0'}{${$_}[0].':'.${$_}[1].':'.${$_}[3]} and exists $anno_SS_ref->{'1'}{${$_}[0].':'.${$_}[2].':'.${$_}[3]}} values %{$SJ_read_improved{$read}};
					if ( @same_SSs == scalar(keys %{$SJ_read_improved{$read}}) ){
						print $default_handle "$read\t$read_info{$read}{'ori_ID'}\t",scalar(@same_SSs), "\tNIC\tNA\t" if $should_create_default_handle;
						$nic{$read} ++;
						$chr_type{$chr}{'nic'} ++;
					} else {
						print $default_handle "$read\t$read_info{$read}{'ori_ID'}\t",scalar(@same_SSs), "\tNNC\tNA\t" if $should_create_default_handle;
						$nnc{$read} ++;
						$chr_type{$chr}{'nnc'} ++;
					}
				} else {
					my @fsm_isoforms;
					my @ism_isoforms;

					my $read_info_ref = $read_info{$read};
					my %possible_isoforms;
					while (my ($SJ2, undef) = each %{$SJ_read_improved{$read}}) {
						die "not_recorded:$chr\t$read\t$SJ2\t\n" if !exists $anno_SJ_ref->{$SJ2};
						die "not_array:$chr\t$read\t$SJ2\t\n" if !defined $anno_SJ_ref->{$SJ2}[-1];
						$possible_isoforms{$_}++ for @{$anno_SJ_ref->{$SJ2}[-1]};
					}
					
					while (my ($isoform, undef) = each %possible_isoforms){
						my $SJ_isoform_ref = $isoform_SJ_ref->{$isoform};
						my $comp_result = &comp_minimap_reference_SJ($SJ_read_improved{$read}, $SJ_isoform_ref, $default_handle);
						if ($comp_result == 1){
							push @fsm_isoforms, $isoform;

						} elsif ($comp_result == 2) {
							push @ism_isoforms, $isoform;
						}
					}
					if (@fsm_isoforms >= 1 ) {
						if ($should_create_default_handle) {
							print $default_handle "$read\t$read_info{$read}{'ori_ID'}\tall\tFSM\t";
							print $default_handle "$_," for @fsm_isoforms;
							print $default_handle "\t";
						}
						push @{$fsm{$read}}, $_ for @fsm_isoforms;
						$chr_type{$chr}{'fsm'} ++;
						for my $fsm_isoform (@fsm_isoforms) {
							if ( abs(${$read_info_ref}{'start'}-$multi_exon_isoform_end_ref->{$fsm_isoform}{'0'}) <= $SJ_dist and abs(${$read_info_ref}{'end'}-$multi_exon_isoform_end_ref->{$fsm_isoform}{'1'}) <= $SJ_dist ) {
								push @{$isoform_perfect_match{$fsm_isoform}{${$read_info_ref}{'sample'}}}, $read;
							}
						}

					} elsif (@ism_isoforms >= 1){
						if ($should_create_default_handle) {
							print $default_handle "$read\t$read_info{$read}{'ori_ID'}\tall\tISM\t";
							print $default_handle "$_," for @ism_isoforms;
							print $default_handle "\t";
						}
						push @{$ism{$read}}, $_ for @ism_isoforms;
						$chr_type{$chr}{'ism'} ++;
					} else {
						print $default_handle "$read\tall\tNIC\tNA\t" if $should_create_default_handle;
						$nic{$read}++;
						$chr_type{$chr}{'nic'} ++;
					}

				}
				print $default_handle "$chr\n" if $should_create_default_handle;
			} else {
				push @single_exon_reads, $read;
				print $default_handle "$read\t$read_info{$read}{'ori_ID'}\tall\tone-exon\tNA\t$chr\n" if $should_create_default_handle;
			}
		}

		my %all_SJ_group;
		
		my ($group_start, $group_end) = @{$group_ends{$n}};
		print $default_handle "group $n: ($group_start, $group_end)\n" if $should_create_default_handle;
		for my $i ($start_SJ_index_anno .. $#annotated_SJ_chr) {
			if ($anno_SJ_ref->{$annotated_SJ_chr[$i]}[1] > $group_end) {
				last;
			} elsif ($anno_SJ_ref->{$annotated_SJ_chr[$i]}[2] >= $group_start and $anno_SJ_ref->{$annotated_SJ_chr[$i]}[1] <= $group_end) {
				$all_SJ_group{$annotated_SJ_chr[$i]} = [@{$anno_SJ_ref->{$annotated_SJ_chr[$i]}}, 0];
			} elsif ($anno_SJ_ref->{$annotated_SJ_chr[$i]}[2] < $group_start){
				$start_SJ_index_anno = $i;
			}
		}
		for my $i ($start_SJ_index_unanno .. $#unannotated_SJ_chr) {
			if ($unanno_SJ{$unannotated_SJ_chr[$i]}[1] > $group_end) {
				last;
			} elsif ($unanno_SJ{$unannotated_SJ_chr[$i]}[2] >= $group_start and $unanno_SJ{$unannotated_SJ_chr[$i]}[1] <= $group_end) {

				$all_SJ_group{$unannotated_SJ_chr[$i]} = [@{$unanno_SJ{$unannotated_SJ_chr[$i]}}, 1];
			} elsif ($unanno_SJ{$unannotated_SJ_chr[$i]}[2] < $group_start){
				$start_SJ_index_unanno = $i;
			}
		}

		my %annotated_isoforms_included;
		while (my ($SJ, $SJ_info) = each %all_SJ_group) {
			if (${$SJ_info}[-1] == 0) {
				for my $isoform (@{${$SJ_info}[-2]}) {
					$annotated_isoforms_included{$isoform} ++;
				}
			}
		}

		while (my ($isoform, undef) = each %annotated_isoforms_included) {
			while (my ($SJ, undef) = each %{$isoform_SJ_ref->{$isoform}}){
				$all_SJ_group{$SJ} = [@{$anno_SJ_ref->{$SJ}}, 0] if !exists $all_SJ_group{$SJ};
			}
		}

		my @all_SJ_group_sort = sort {$all_SJ_group{$a}[1] <=> $all_SJ_group{$b}[1] or $all_SJ_group{$a}[2] <=> $all_SJ_group{$b}[2]} keys %all_SJ_group;
		my %all_SJ_group_ID;
		$all_SJ_group_ID{$all_SJ_group_sort[$_]} = $_ for 0 .. $#all_SJ_group_sort;

		my %fsm_chain;
		while (my ($isoform, undef) = each %annotated_isoforms_included) {
			my @SJ_in_group_ID_sort = sort {$a <=> $b} (map {$all_SJ_group_ID{$_}} keys %{$isoform_SJ_ref->{$isoform}});
			$fsm_chain{$isoform} = join '_', @SJ_in_group_ID_sort;
			print $default_handle "annotated_isoform: $isoform\t$fsm_chain{$isoform}\n" if $should_create_default_handle;
		}

		my (%SJ_chain_read, %SJ_chain_cat);
		for my $read (@{$group_reads{$n}}) {
			if (exists $SJ_read_order{$read}) {
				my $SJ_count = scalar (keys %{$SJ_read_order{$read}});
				if (exists $SJ_read_order{$read}{'-1'}){
					$SJ_read_order{$read}{$SJ_count} = $SJ_read_order{$read}{'-1'};
					delete $SJ_read_order{$read}{'-1'};
				}
				my @SJ_order_sort = sort {$a <=> $b} keys %{$SJ_read_order{$read}};
				my @SJ_IDs;
				for my $SJ_order (@SJ_order_sort) {
					if ($SJ_read_order{$read}{$SJ_order} ne 'x'){
						if (exists $all_SJ_group_ID{$SJ_read_order{$read}{$SJ_order}}) {
							push @SJ_IDs, $all_SJ_group_ID{$SJ_read_order{$read}{$SJ_order}}
						} else {
							die "$read has no SJ ID for: $SJ_read_order{$read}{$SJ_order}";
						}
						
					} else {
						push @SJ_IDs, 'x';
					}
				}
				my $SJ_chain = join '_', @SJ_IDs;
				
				if (!exists $SJ_chain_read{$SJ_chain}) {
					my $cat;
					if (exists $fsm{$read}) {
						$SJ_chain_cat{'fsm'}{$SJ_chain} = [scalar(@SJ_IDs), $fsm{$read}];
						$cat = 'fsm';

					} elsif (exists $ism{$read}) {
						$SJ_chain_cat{'ism'}{$SJ_chain} = [scalar(@SJ_IDs), $ism{$read}];
						$cat = 'ism';
					} elsif (exists $nic{$read}) {
						$SJ_chain_cat{'nic'}{$SJ_chain} = [scalar(@SJ_IDs)];
						$cat = 'nic';
					} elsif (exists $nnc{$read}) {
						$SJ_chain_cat{'nnc'}{$SJ_chain} = [scalar(@SJ_IDs)];
						$cat = 'nnc';
					} else {
						$SJ_chain_cat{'filtered'}{$SJ_chain} = [scalar(@SJ_IDs)];
						$cat = 'filtered';
					}
					$SJ_chain_read{$SJ_chain} = [$cat, $read_info{$read}{'start'}, $read_info{$read}{'end'}, $read_info{$read}{'strand_isoform'}, scalar(@SJ_IDs), {}, [[$read, $read_info{$read}{'start'}, $read_info{$read}{'end'}]]];
				} else {
					push @{$SJ_chain_read{$SJ_chain}[-1]}, [$read, $read_info{$read}{'start'}, $read_info{$read}{'end'}];

				}
			}
		}


		my (%SJ_chain_read_exon_end, %SJ_included_chain);
		while (my ($SJ_chain, $SJ_chain_info) = each %SJ_chain_read) {
			#
			my @SJ_IDs = split '_', $SJ_chain;
			my @start_reads_sort = sort {${$a}[1] <=> ${$b}[1]} @{${$SJ_chain_info}[-1]};
			my @end_reads_sort = sort {${$a}[2] <=> ${$b}[2]} @{${$SJ_chain_info}[-1]};
			$SJ_chain_read{$SJ_chain}[1] = ${$start_reads_sort[int(@start_reads_sort/2)]}[1];
			$SJ_chain_read{$SJ_chain}[2] = ${$end_reads_sort[int(@end_reads_sort/2)]}[2];
			next if defined $raw;
			for my $SJ_ID_i (0 .. $#SJ_IDs) {
				push @{$SJ_included_chain{$SJ_IDs[$SJ_ID_i]}}, $SJ_chain if $SJ_IDs[$SJ_ID_i] ne 'x'; #record chain in to hash with keys as SJ
				next if exists $SJ_chain_cat{'filtered'}{$SJ_chain};
				my ($start, $end);
				if ($SJ_ID_i == 0) {
					$start = ${$SJ_chain_info}[1];
				}
				if ($SJ_ID_i == $#SJ_IDs) {
					$end = ${$SJ_chain_info}[2];
				}
				if(!defined $start) {
					my $pre_SJ = $all_SJ_group_sort[$SJ_IDs[$SJ_ID_i-1]];
					my ($pre_SJ_chr, $pre_SJ_start, $pre_SJ_end, $pre_SJ_strand) = split ':', $pre_SJ;
					$start = $pre_SJ_end;
				} if (!defined $end) {
					my $next_SJ = $all_SJ_group_sort[$SJ_IDs[$SJ_ID_i+1]];
					my ($next_SJ_chr, $next_SJ_start, $next_SJ_end, $next_SJ_strand) = split ':', $next_SJ;
					$end = $next_SJ_start;
				}
				$SJ_chain_read_exon_end{$SJ_chain}{'start'}{$SJ_IDs[$SJ_ID_i]} = $start;
				$SJ_chain_read_exon_end{$SJ_chain}{'end'}{$SJ_IDs[$SJ_ID_i]} = $end;
				print $default_handle "$chr:exon_end: $SJ_chain\t$SJ_IDs[$SJ_ID_i]\t$start\t$end\n" if $should_create_default_handle;
			}
		}

		while (my ($SJ_ID, $SJ_chain_ref) = each %SJ_included_chain) {
			my @chain_length_sort = sort {$SJ_chain_read{$a}[4] <=> $SJ_chain_read{$b}[4]} @{$SJ_chain_ref}; #record chain in to hash with keys as SJ
			for my $i( 0 .. $#chain_length_sort-1 ){        # record possible longer chains of a shorter(filtered) chain
				my $chain_short = $chain_length_sort[$i];
				for my $j ($i+1 .. $#chain_length_sort){
					my $chain_long = $chain_length_sort[$j];
					if ($SJ_chain_read{$chain_short}[4] == $SJ_chain_read{$chain_long}[4]){
						if (exists $SJ_chain_cat{'filtered'}{$chain_long} and !exists $SJ_chain_cat{'filtered'}{$chain_short}) {
							${$SJ_chain_read{$chain_long}[-2]}{$chain_short}++;
						} elsif (exists $SJ_chain_cat{'filtered'}{$chain_short} and !exists $SJ_chain_cat{'filtered'}{$chain_long}) {
							${$SJ_chain_read{$chain_short}[-2]}{$chain_long}++;
						}
					} elsif (!exists $SJ_chain_cat{'filtered'}{$chain_long}) {
						${$SJ_chain_read{$chain_short}[-2]}{$chain_long}++;
					}
				}
			}
		}

		my (%longest_nc_valid, @nc_valid_longest, %substring_nc, @nic_valid_all, @nnc_valid_all);
		my (%read_count_isoform, %cat_SJ_chains_sort, %possible_fsm_filtered, %unannotated_isoforms);

		if (defined $raw) {
			for my $cat ('ism', 'nic', 'nnc') {
				if (exists $SJ_chain_cat{$cat}){
					while (my ($chain, undef) = each %{$SJ_chain_cat{$cat}}) {
						my $validation_result = &valid_nc($default_handle, $chain, $chr_SJ_reads_ref, $SJ_chain_read{$chain}[-1], \%read_info, \@all_SJ_group_sort, $anno_SJ_complete_ref, $isoform_SJ_complete_ref, $read_num_cutoff, $read_ratio_cutoff);
						$unannotated_isoforms{$chain} = $validation_result;
						$longest_nc_valid{$chain} = -1 if ${$validation_result}[0] > 0;
					}
				}
			}
		} else {


		
		if (exists $SJ_chain_cat{'nnc'}) {
			my @chain_length_sort = sort {$SJ_chain_cat{'nnc'}{$a}[0] <=> $SJ_chain_cat{'nnc'}{$b}[0]} keys %{$SJ_chain_cat{'nnc'}};
			for my $chain(@chain_length_sort) {
				my $validation_result = &valid_nc($default_handle, $chain, $chr_SJ_reads_ref, $SJ_chain_read{$chain}[-1], \%read_info, \@all_SJ_group_sort, $anno_SJ_complete_ref, $isoform_SJ_complete_ref, $read_num_cutoff, $read_ratio_cutoff);
				if (${$validation_result}[0] >= 1) { #chr_anno_SJ_ref
					push @nnc_valid_all, $chain;
					$unannotated_isoforms{$chain} = $validation_result; # if @{${$validation_result}[1]} > 0;
				}
			}
			for my $i( 0 .. $#nnc_valid_all-1 ){
				my @SJs_i = split '_', $nnc_valid_all[$i];
				for my $j ($i+1 .. $#nnc_valid_all){
					next unless exists ${$SJ_chain_read{$nnc_valid_all[$i]}[-2]}{$nnc_valid_all[$j]};
					my $comp_result = &comp_SJ_chain($nnc_valid_all[$i], $nnc_valid_all[$j]);
					print $default_handle "$nnc_valid_all[$i]\t$nnc_valid_all[$j]\t$comp_result\t$SJ_chain_read_exon_end{$nnc_valid_all[$j]}{'start'}{$SJs_i[0]}\t$SJ_chain_read_exon_end{$nnc_valid_all[$j]}{'end'}{$SJs_i[-1]}\t$SJ_chain_read{$nnc_valid_all[$i]}[1]\t$SJ_chain_read{$nnc_valid_all[$i]}[2]\n" if $should_create_default_handle;
					if ($comp_result == 2) {
						push @{$substring_nc{$nnc_valid_all[$i]}}, $nnc_valid_all[$j] if $SJ_chain_read_exon_end{$nnc_valid_all[$j]}{'start'}{$SJs_i[0]} - $SJ_dist <= $SJ_chain_read{$nnc_valid_all[$i]}[1] and $SJ_chain_read_exon_end{$nnc_valid_all[$j]}{'end'}{$SJs_i[-1]} + $SJ_dist >= $SJ_chain_read{$nnc_valid_all[$i]}[2];
					}
				}
			}
			for my $chain(@nnc_valid_all) {
				if (!exists $substring_nc{$chain}) {
					push @nc_valid_longest, $chain;
					$longest_nc_valid{$chain} = $#nc_valid_longest;
				}
			}
		}
		if (exists $SJ_chain_cat{'nic'}) {
			my @chain_length_sort = sort {$SJ_chain_cat{'nic'}{$a}[0] <=> $SJ_chain_cat{'nic'}{$b}[0]} keys %{$SJ_chain_cat{'nic'}};
			for my $chain(@chain_length_sort) {
				my $validation_result = &valid_nc($default_handle, $chain, $chr_SJ_reads_ref, $SJ_chain_read{$chain}[-1], \%read_info, \@all_SJ_group_sort, $anno_SJ_complete_ref, $isoform_SJ_complete_ref, $read_num_cutoff, $read_ratio_cutoff);
				if (${$validation_result}[0] >= 1) {
					push @nic_valid_all, $chain;
					$unannotated_isoforms{$chain} = $validation_result; 
				}
			}

			for my $i( 0 .. $#nic_valid_all-1 ){
				my @SJs_i = split '_', $nic_valid_all[$i];
				for my $j ($i+1 .. $#nic_valid_all){
					next unless exists ${$SJ_chain_read{$nic_valid_all[$i]}[-2]}{$nic_valid_all[$j]};
					my $comp_result = &comp_SJ_chain($nic_valid_all[$i], $nic_valid_all[$j]);
					print $default_handle "$nic_valid_all[$i]\t$nic_valid_all[$j]\t$comp_result\t$SJ_chain_read_exon_end{$nic_valid_all[$j]}{'start'}{$SJs_i[0]}\t$SJ_chain_read_exon_end{$nic_valid_all[$j]}{'end'}{$SJs_i[-1]}\t$SJ_chain_read{$nic_valid_all[$i]}[1]\t$SJ_chain_read{$nic_valid_all[$i]}[2]\n" if $should_create_default_handle;
					if ($comp_result == 2) {
						push @{$substring_nc{$nic_valid_all[$i]}}, $nic_valid_all[$j] if $SJ_chain_read_exon_end{$nic_valid_all[$j]}{'start'}{$SJs_i[0]}- $SJ_dist <= $SJ_chain_read{$nic_valid_all[$i]}[1] and $SJ_chain_read_exon_end{$nic_valid_all[$j]}{'end'}{$SJs_i[-1]} + $SJ_dist >= $SJ_chain_read{$nic_valid_all[$i]}[2];
					}
				}
			}
			my %longest_nic_valid;
			for my $chain(@nic_valid_all) {
				my @SJs_i = split '_', $chain;
				while (my ($nc_chain, undef) = each %longest_nc_valid) {
					next unless exists ${$SJ_chain_read{$chain}[-2]}{$nc_chain};
					my $comp_result = &comp_SJ_chain($chain, $nc_chain);
					if ($comp_result == 2) {
						push @{$substring_nc{$chain}}, $nc_chain if $SJ_chain_read_exon_end{$nc_chain}{'start'}{$SJs_i[0]} - $SJ_dist <= $SJ_chain_read{$chain}[1] and $SJ_chain_read_exon_end{$nc_chain}{'end'}{$SJs_i[-1]} + $SJ_dist >= $SJ_chain_read{$chain}[2];
					}
				}
				if (!exists $substring_nc{$chain}) {
					push @nc_valid_longest, $chain;
					$longest_nic_valid{$chain} = $#nc_valid_longest;
				}
			}
			$longest_nc_valid{$_} = $longest_nic_valid{$_} for keys %longest_nic_valid;
		}

		if (exists $SJ_chain_cat{'ism'}) {
			my @chain_length_sort = sort {$SJ_chain_cat{'ism'}{$a}[0] <=> $SJ_chain_cat{'ism'}{$b}[0]} keys %{$SJ_chain_cat{'ism'}};
			my (@possible_nc_ism, %chain_w_compatible_annotated_fsm);
			for my $chain(@chain_length_sort) {
				my (@isoforms_fsm_valid);
				my @SJs_i = split '_', $chain;
				if ($should_create_default_handle) {
					print $default_handle "ism\t$chain\t";
					print $default_handle "$_," for @{$SJ_chain_cat{'ism'}{$chain}[1]};
				}
				my $compatible_annotated_fsm_num = 0;
				for my $isoform_fsm(@{$SJ_chain_cat{'ism'}{$chain}[1]}) {
					if (exists $multi_exon_isoform_end_ref->{$isoform_fsm}){
						if( exists $multi_exon_isoform_end_ref->{$isoform_fsm}{$anno_SJ_ref->{$all_SJ_group_sort[$SJs_i[0]]}[1]} and exists $multi_exon_isoform_end_ref->{$isoform_fsm}{$anno_SJ_ref->{$all_SJ_group_sort[$SJs_i[-1]]}[2]} and $multi_exon_isoform_end_ref->{$isoform_fsm}{$anno_SJ_ref->{$all_SJ_group_sort[$SJs_i[0]]}[1]} - $SJ_dist <= $SJ_chain_read{$chain}[1] and $multi_exon_isoform_end_ref->{$isoform_fsm}{$anno_SJ_ref->{$all_SJ_group_sort[$SJs_i[-1]]}[2]} + $SJ_dist >= $SJ_chain_read{$chain}[2] ){
							$compatible_annotated_fsm_num ++;
							push @isoforms_fsm_valid, $isoform_fsm if exists $SJ_chain_read{$fsm_chain{$isoform_fsm}};	#check if there is any full length annotated isoform for this chain
						}
					}
				}
				if ($should_create_default_handle) {
					print $default_handle "\t";
					print $default_handle "$_," for @isoforms_fsm_valid;
					print $default_handle "\n";
				}
				$SJ_chain_cat{'ism'}{$chain}[1] = \@isoforms_fsm_valid;
				if (@isoforms_fsm_valid == 0) { #check if the ism chain is valid if no full length annotated isoform contains the chain
					$chain_w_compatible_annotated_fsm{$chain} = 1 if $compatible_annotated_fsm_num > 0; #check if completely compatible with annotated isoform
					my $validation_result = &valid_nc($default_handle, $chain, $chr_SJ_reads_ref, $SJ_chain_read{$chain}[-1], \%read_info, \@all_SJ_group_sort, $anno_SJ_complete_ref, $isoform_SJ_complete_ref, $read_num_cutoff, $read_ratio_cutoff);
					if (${$validation_result}[0] >= 1) {
						push @possible_nc_ism, $chain;
						$unannotated_isoforms{$chain} = $validation_result; # if @{${$validation_result}[1]} > 0;
					}
				}
			}
			for my $i ( 0 .. $#possible_nc_ism-1 ){ #pairwise comaprison of valid ism chains without full length annotated isoform to find the longest
				my @SJs_i = split '_', $possible_nc_ism[$i];
				for my $j ($i+1 .. $#possible_nc_ism){
					next unless exists ${$SJ_chain_read{$possible_nc_ism[$i]}[-2]}{$possible_nc_ism[$j]};
					my $comp_result = &comp_SJ_chain($possible_nc_ism[$i], $possible_nc_ism[$j]);
					print $default_handle "$possible_nc_ism[$i]\t$possible_nc_ism[$j]\t$comp_result\t$SJ_chain_read_exon_end{$possible_nc_ism[$j]}{'start'}{$SJs_i[0]}\t$SJ_chain_read_exon_end{$possible_nc_ism[$j]}{'end'}{$SJs_i[-1]}\t$SJ_chain_read{$possible_nc_ism[$i]}[1]\t$SJ_chain_read{$possible_nc_ism[$i]}[2]\n" if $should_create_default_handle;
					if ($comp_result == 2) {
						push @{$substring_nc{$possible_nc_ism[$i]}}, $possible_nc_ism[$j] if $SJ_chain_read_exon_end{$possible_nc_ism[$j]}{'start'}{$SJs_i[0]} - $SJ_dist <= $SJ_chain_read{$possible_nc_ism[$i]}[1] and $SJ_chain_read_exon_end{$possible_nc_ism[$j]}{'end'}{$SJs_i[-1]} + $SJ_dist >= $SJ_chain_read{$possible_nc_ism[$i]}[2];
					}
				}
			}

			for my $chain (@possible_nc_ism) {
				my @SJs_i = split '_', $chain;
				for my $nc_chain (@nc_valid_longest){   #compare valid ism chains without full length annotated isoform to valid longest nc chians
					next unless exists ${$SJ_chain_read{$chain}[-2]}{$nc_chain};
					my $comp_result = &comp_SJ_chain($chain, $nc_chain);
					if ($comp_result == 2) {
						push @{$substring_nc{$chain}}, $nc_chain if $SJ_chain_read_exon_end{$nc_chain}{'start'}{$SJs_i[0]} - $SJ_dist <= $SJ_chain_read{$chain}[1] and $SJ_chain_read_exon_end{$nc_chain}{'end'}{$SJs_i[-1]} + $SJ_dist >= $SJ_chain_read{$chain}[2];
					}
				}
				if (exists $substring_nc{$chain}) {
					my @nc_isoforms;
					for my $nc_chain (@{$substring_nc{$chain}}) {
						if (exists $longest_nc_valid{$nc_chain}) {
							push @nc_isoforms, $longest_nc_valid{$nc_chain};
						}
					}
					$SJ_chain_cat{'ism'}{$chain}[1] = \@nc_isoforms;
				} elsif (exists $chain_w_compatible_annotated_fsm{$chain}) {
					push @nc_valid_longest, $chain;
					$longest_nc_valid{$chain} = $#nc_valid_longest;
				}
			}
		}

		}
		
		while (my ($SJ_chain, $chain_info_ref) = each %SJ_chain_read) {
			my ($cat, $chain_start, $chain_end, $strand, $SJ_count, $related_chain_ref, $read_info) = @{$chain_info_ref};
			my %read_count_sample;
			for my $read_info_ref (@{$read_info}){
				$read_count_sample{$read_info{${$read_info_ref}[0]}{'sample'}} ++;
			}
			for my $sample (@{$samples_sort_ref}) {
				if (exists $read_count_sample{$sample}) {
					$read_count_isoform{$SJ_chain}{$sample} = [$read_count_sample{$sample},$read_count_sample{$sample},0,0,0]; #before_assign, direct_assign, EM_assign, current_increase, last_increase
				} else {
					$read_count_isoform{$SJ_chain}{$sample} = [0,0,0,0,0]; 
				}				
			}

		}
		
		my %ism_mt1_chain_count;
		if (!defined $raw) {
		
		for my $cat ('ism', 'fsm', 'nic', 'nnc', 'filtered') {
			@{$cat_SJ_chains_sort{$cat}} = sort {$SJ_chain_cat{$cat}{$b}[0] <=> $SJ_chain_cat{$cat}{$a}[0]} keys %{$SJ_chain_cat{$cat}};
		}

		for my $cat ('nic', 'nnc') {
			for my $chain (@{$cat_SJ_chains_sort{$cat}}) {
				my $type;
				if (exists $substring_nc{$chain}) {
					my @nc_isoforms;
					for my $nc_chain (@{$substring_nc{$chain}}) {
						if (exists $longest_nc_valid{$nc_chain}) {
							push @nc_isoforms, $longest_nc_valid{$nc_chain};
						}
					}
					
					if (@nc_isoforms>0) {
						push @{$SJ_chain_cat{$cat}{$chain}}, \@nc_isoforms;
						if (scalar(@nc_isoforms) == 1){
							my $read_count_total=0;
							for my $sample (@{$samples_sort_ref}) {
								$read_count_isoform{ $nc_valid_longest[$nc_isoforms[0]] }{$sample}[1] += $read_count_isoform{$chain}{$sample}[0];
								$read_count_total += $read_count_isoform{$chain}{$sample}[0];
							}
							$type = "direct";
							printf $default_handle "!!read_count_summary\t$chr\t$cat\t$n:$chain\t$type\t%g!=$read_count_total\n", scalar(@{$SJ_chain_read{$chain}[-1]}) if $should_create_default_handle and $read_count_total!=@{$SJ_chain_read{$chain}[-1]};
						} else {
							$type = "EM";
						}
						if ($should_create_tsv_compt_handle){
							for my $read_info_ref (@{$SJ_chain_read{$chain}[-1]}) {
								my $read_ID = ${$read_info_ref}[0];
								print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tNIC/NNC\t";
								print $tsv_compt_handle "ESPRESSO:$chr:$n:$_," for @{$SJ_chain_cat{$cat}{$chain}[1]};
								print $tsv_compt_handle "\n";
							}
						}
						printf $default_handle "read_count_summary\t$chr\t$cat\t$n:$chain\t$type\t%g\n", scalar(@{$SJ_chain_read{$chain}[-1]}) if $should_create_default_handle;
					} else {
						if ($should_create_tsv_compt_handle){
							for my $read_info_ref (@{$SJ_chain_read{$chain}[-1]}) {
								my $read_ID = ${$read_info_ref}[0];
								print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tNIC/NNC\tNA\n";
							}
						}
						printf $default_handle "read_count_summary\t$chr\t$cat\t$n:$chain\tNA1\t%g\n", scalar(@{$SJ_chain_read{$chain}[-1]}) if $should_create_default_handle;
					}
				} elsif (exists $longest_nc_valid{$chain}) {
					if ($should_create_tsv_compt_handle){
						for my $read_info_ref (@{$SJ_chain_read{$chain}[-1]}){
							my $read_ID = ${$read_info_ref}[0];
							print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tNIC/NNC\tESPRESSO:$chr:$n:$longest_nc_valid{$chain}\n";
						}
					}
					
					printf $default_handle "read_count_summary\t$chr\t$cat\t$n:$chain\titself\t%g\n", scalar(@{$SJ_chain_read{$chain}[-1]}) if $should_create_default_handle;
				} else {
					if ($should_create_tsv_compt_handle){
						for my $read_info_ref (@{$SJ_chain_read{$chain}[-1]}) {
							my $read_ID = ${$read_info_ref}[0];
							print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tNIC/NNC\tNA\n";
						}
					}
					printf $default_handle "read_count_summary\t$chr\t$cat\t$n:$chain\tNA2\t%g\n", scalar(@{$SJ_chain_read{$chain}[-1]}) if $should_create_default_handle;
					#print "$cat\t$chain\t",scalar(@{$SJ_chain_read{$chain}[-1]}),"\n";
				}
				
			}
		}

		
		for my $ism_SJ_chain (@{$cat_SJ_chains_sort{'ism'}}) {
			if (exists $longest_nc_valid{$ism_SJ_chain}){
				if ($should_create_tsv_compt_handle){
					for my $read_info_ref (@{$SJ_chain_read{$ism_SJ_chain}[-1]}){
						my $read_ID = ${$read_info_ref}[0];
						print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tNIC/NNC\t";
						print $tsv_compt_handle "ESPRESSO:$chr:$n:$longest_nc_valid{$ism_SJ_chain},";
						print $tsv_compt_handle "\n";
					}
				}
				printf $default_handle "read_count_summary\t$chr\tism\t$n:$ism_SJ_chain\titself\t%g\n", scalar(@{$SJ_chain_read{$ism_SJ_chain}[-1]}) if $should_create_default_handle;
			} else {
				if (!exists $substring_nc{$ism_SJ_chain}) {
					my @SJs_i = split '_', $ism_SJ_chain;
					while (my ($nc_chain, $nc_chain_ID) = each %longest_nc_valid) {
						next unless exists ${$SJ_chain_read{$ism_SJ_chain}[-2]}{$nc_chain};
						my $comp_result = &comp_SJ_chain($ism_SJ_chain, $nc_chain);
						if ($comp_result == 2) {
							push @{$SJ_chain_cat{'ism'}{$ism_SJ_chain}[1]}, $nc_chain_ID if $SJ_chain_read_exon_end{$nc_chain}{'start'}{$SJs_i[0]} - $SJ_dist <= $SJ_chain_read{$ism_SJ_chain}[1] and $SJ_chain_read_exon_end{$nc_chain}{'end'}{$SJs_i[-1]} + $SJ_dist >= $SJ_chain_read{$ism_SJ_chain}[2];
						}
					}
				}
				
				if (scalar(@{$SJ_chain_cat{'ism'}{$ism_SJ_chain}[1]}) == 1) {
					my @isoforms_only = @{$SJ_chain_cat{'ism'}{$ism_SJ_chain}[1]};
					my $isoform_fsm = $SJ_chain_cat{'ism'}{$ism_SJ_chain}[1][0];
					my ($read_compatible_isoform_1);
					if (exists $fsm_chain{$isoform_fsm}) {
						for my $sample (@{$samples_sort_ref}) {
							$read_count_isoform{$fsm_chain{$isoform_fsm}}{$sample}[1] += $read_count_isoform{$ism_SJ_chain}{$sample}[0];
						}
						$read_compatible_isoform_1 = $isoform_fsm;
					} elsif ($isoform_fsm =~ /^\d+$/ and defined $nc_valid_longest[$isoform_fsm]){
						for my $sample (@{$samples_sort_ref}) {
							$read_count_isoform{$nc_valid_longest[$isoform_fsm]}{$sample}[1] += $read_count_isoform{$ism_SJ_chain}{$sample}[0];
						}
						$read_compatible_isoform_1 = "ESPRESSO:$chr:$n:$isoform_fsm";
					} else {
						die "don't know what fsm is for $ism_SJ_chain in group $n($group_start, $group_end)\t$isoform_fsm\n";
					}
					if ($should_create_tsv_compt_handle){
						for my $read_info_ref (@{$SJ_chain_read{$ism_SJ_chain}[-1]}){
							my $read_ID = ${$read_info_ref}[0];
							print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tISM\t";
							print $tsv_compt_handle "$read_compatible_isoform_1,\n";
						}
					}
					printf $default_handle "read_count_summary\t$chr\tism\t$n:$ism_SJ_chain\tdirect\t%g\n", scalar(@{$SJ_chain_read{$ism_SJ_chain}[-1]}) if $should_create_default_handle;
				} elsif (scalar(@{$SJ_chain_cat{'ism'}{$ism_SJ_chain}[1]}) > 1) {
					my @ism_SJ_order = split '_', $ism_SJ_chain;
					my @SJ_chain_pos = map {$all_SJ_group_sort[$_]} @ism_SJ_order;
					my (%chains_string_freq);
					for my $read_info_ref (@{$SJ_chain_read{$ism_SJ_chain}[-1]}){
						my ($read_ID, $read_start, $read_end, $strand) = @{$read_info_ref};
						my $sample = $read_info{$read_ID}{'sample'};
						my (%fsm_SJ_chain_right, %read_compatible_isoform);
						for my $isoform_fsm (@{$SJ_chain_cat{'ism'}{$ism_SJ_chain}[1]}){
							if (exists $multi_exon_isoform_end_ref->{$isoform_fsm}){
								my $fsm_chain_from_isoform = $fsm_chain{$isoform_fsm};
								# Multiple isoforms could have the same chain.
								if (!exists $fsm_SJ_chain_right{$fsm_chain_from_isoform}) {
									$fsm_SJ_chain_right{$fsm_chain_from_isoform} = 0;
								}
								$read_compatible_isoform{$isoform_fsm} = 0;
								if( exists $multi_exon_isoform_end_ref->{$isoform_fsm}{$anno_SJ_ref->{$SJ_chain_pos[0]}[1]} and exists $multi_exon_isoform_end_ref->{$isoform_fsm}{$anno_SJ_ref->{$SJ_chain_pos[-1]}[2]} and $multi_exon_isoform_end_ref->{$isoform_fsm}{$anno_SJ_ref->{$SJ_chain_pos[0]}[1]} - $SJ_dist <= $read_start and $multi_exon_isoform_end_ref->{$isoform_fsm}{$anno_SJ_ref->{$SJ_chain_pos[-1]}[2]} + $SJ_dist >= $read_end ){
									$fsm_SJ_chain_right{$fsm_chain_from_isoform} ++;
									$read_compatible_isoform{$isoform_fsm} ++;
								}
							} elsif ($isoform_fsm =~ /^\d+$/ and defined $nc_valid_longest[$isoform_fsm]){
								$fsm_SJ_chain_right{$nc_valid_longest[$isoform_fsm]} ++;
								$read_compatible_isoform{"ESPRESSO:$chr:$n:$isoform_fsm"} ++;

							} else {
								die "don't know what fsm is for2 $ism_SJ_chain in group $n($group_start, $group_end): $isoform_fsm";
							}
						}
						my @sort_chains;
						my @valid_SJ_chain = grep {$fsm_SJ_chain_right{$_}>0} keys %fsm_SJ_chain_right;
						if (scalar(@valid_SJ_chain)>=1) {
							@sort_chains = sort {$a cmp $b} @valid_SJ_chain;
							print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tISM\t" if $should_create_tsv_compt_handle;
							my @valid_isoforms = grep {$read_compatible_isoform{$_}>0} keys %read_compatible_isoform;
							if ($should_create_tsv_compt_handle) {
								print $tsv_compt_handle "$_," for @valid_isoforms;
								print $tsv_compt_handle "\n";
							}
						} else {
							@sort_chains = sort {$a cmp $b} keys %fsm_SJ_chain_right;
							if ($should_create_tsv_compt_handle) {
								print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tISM\t";
								print $tsv_compt_handle "$_," for keys %read_compatible_isoform;
								print $tsv_compt_handle "\n";
							}
						}

						my $chains_compt_read_string = join(':', @sort_chains);
						$chains_string_freq{$chains_compt_read_string}{$sample} ++;
					}
					my $total_count = 0;
					while (my ($chains_string, $sample_ref) = each %chains_string_freq) {
						while (my ($sample, $count) = each %{$sample_ref}) {
							push @{$ism_mt1_chain_count{$ism_SJ_chain}{$sample}}, [$chains_string, $count];
							$total_count += $count;
						}
					}
					die "$total_count!=@{$SJ_chain_read{$ism_SJ_chain}[-1]}" if $total_count!=@{$SJ_chain_read{$ism_SJ_chain}[-1]};
					printf $default_handle "read_count_summary\t$chr\tism\t$n:$ism_SJ_chain\tEM\t%g\n", scalar(@{$SJ_chain_read{$ism_SJ_chain}[-1]}) if $should_create_default_handle;
				} else {
					if ($should_create_tsv_compt_handle){
						for my $read_info_ref (@{$SJ_chain_read{$ism_SJ_chain}[-1]}){
							my $read_ID = ${$read_info_ref}[0];
							print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tISM\t";
							print $tsv_compt_handle "NA\n";
						}
					}
					
					printf $default_handle "read_count_summary\t$chr\tism\t$n:$ism_SJ_chain\tNA0\t%g\n", scalar(@{$SJ_chain_read{$ism_SJ_chain}[-1]}) if $should_create_default_handle;
				}
			}
		}
		
		for my $filtered_SJ_chain (@{$cat_SJ_chains_sort{'filtered'}}) {
			my @SJs_i = split '_', $filtered_SJ_chain;
			my (@SJ_chains, @isoforms_compatible);
			for my $fsm_SJ_chain (@{$cat_SJ_chains_sort{'fsm'}}) {
				next unless exists ${$SJ_chain_read{$filtered_SJ_chain}[-2]}{$fsm_SJ_chain};
				my $comp_result = &comp_SJ_chain($filtered_SJ_chain, $fsm_SJ_chain);
				if ($comp_result == 3) {
					if ($SJ_chain_read_exon_end{$fsm_SJ_chain}{'start'}{$SJs_i[0]} - $SJ_dist <= $SJ_chain_read{$filtered_SJ_chain}[1] and $SJ_chain_read_exon_end{$fsm_SJ_chain}{'end'}{$SJs_i[-1]} + $SJ_dist >= $SJ_chain_read{$filtered_SJ_chain}[2]){
						push @SJ_chains, $fsm_SJ_chain;
						my @isoforms = @{$SJ_chain_cat{'fsm'}{$fsm_SJ_chain}[1]};
						push @isoforms_compatible, @isoforms;;
					}
				}
			}
			while (my ($nc_chain, $nc_chain_ID) = each %longest_nc_valid) {
				next unless exists ${$SJ_chain_read{$filtered_SJ_chain}[-2]}{$nc_chain};
				my $comp_result = &comp_SJ_chain($filtered_SJ_chain, $nc_chain);
				if ($comp_result == 3) {
					if ($SJ_chain_read_exon_end{$nc_chain}{'start'}{$SJs_i[0]} - $SJ_dist <= $SJ_chain_read{$filtered_SJ_chain}[1] and $SJ_chain_read_exon_end{$nc_chain}{'end'}{$SJs_i[-1]} + $SJ_dist >= $SJ_chain_read{$filtered_SJ_chain}[2]){
						push @SJ_chains, $nc_chain;
						push @isoforms_compatible, "ESPRESSO:$chr:$n:$longest_nc_valid{$nc_chain}";
					}
				}
			}
			if (@SJ_chains >= 1) {
				$possible_fsm_filtered{scalar(@SJ_chains)}{$filtered_SJ_chain} = \@SJ_chains;
				if ($should_create_tsv_compt_handle){
					for my $read_info_ref (@{$SJ_chain_read{$filtered_SJ_chain}[-1]}) {
						my $read_ID = ${$read_info_ref}[0];
						print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tNCD\t";
						for my $isoform (@isoforms_compatible) {
							#my @isoforms = @{$SJ_chain_cat{'fsm'}{$chain}[1]};
							print $tsv_compt_handle "$isoform,";
						}
						print $tsv_compt_handle "\n";
					}
				}
				printf $default_handle "read_count_summary\t$chr\tfiltered\t$n:$filtered_SJ_chain\tEM\t%g\n", scalar(@{$SJ_chain_read{$filtered_SJ_chain}[-1]}) if $should_create_default_handle and @SJ_chains > 1;
			} else {
				for my $read_info_ref (@{$SJ_chain_read{$filtered_SJ_chain}[-1]}) {
					my $read_ID = ${$read_info_ref}[0];
					print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tNCD\tNA\n" if $should_create_tsv_compt_handle;
				}
				printf $default_handle "read_count_summary\t$chr\tfiltered\t$n:$filtered_SJ_chain\tNA\t%g\n", scalar(@{$SJ_chain_read{$filtered_SJ_chain}[-1]}) if $should_create_default_handle;
			}
		}
		while(my ($filtered_SJ_chain, $chain_ref) = each %{$possible_fsm_filtered{'1'}}){
			for my $sample (@{$samples_sort_ref}) {
				$read_count_isoform{${$chain_ref}[0]}{$sample}[1] += $read_count_isoform{$filtered_SJ_chain}{$sample}[0];
			}
			printf $default_handle "read_count_summary\t$chr\tfiltered\t$n:$filtered_SJ_chain\tdirect\t%g\n", scalar(@{$SJ_chain_read{$filtered_SJ_chain}[-1]}) if $should_create_default_handle;
		}
		while (my ($chain, $sample_ref) = each %read_count_isoform) {
			while (my ($sample, $count_ref) = each %{$sample_ref}) {
				$read_count_isoform{$chain}{$sample}[2]=$read_count_isoform{$chain}{$sample}[1];
			}
		}



		}

		if ($should_create_default_handle) {
			while (my ($SJ_chain, $chain_info_ref) = each %SJ_chain_read) {
				my @SJ_IDs = split '_', $SJ_chain;
				my ($cat, $chain_start, $chain_end, $strand, $SJ_count, $related_chain_ref, $read_info) = @{$chain_info_ref};
				print $default_handle "isoform_chain\t$chr\t$SJ_chain\t$cat\t", scalar(@{$read_info}), "\t$SJ_chain_read{$SJ_chain}[1]\t$SJ_chain_read{$SJ_chain}[2]\t";

				for my $SJ_ID (@SJ_IDs){
					if ($SJ_ID eq 'x'){
						print $default_handle "x,";
					} else {
						print $default_handle "$all_SJ_group_sort[$SJ_ID],";
					}
				}
				print $default_handle "\t";
				if (exists $SJ_chain_cat{$cat}{$SJ_chain} and defined $SJ_chain_cat{$cat}{$SJ_chain}[1] and @{$SJ_chain_cat{$cat}{$SJ_chain}[1]} > 0) {
					print $default_handle "$_," for @{$SJ_chain_cat{$cat}{$SJ_chain}[1]};
				} else {
					print $default_handle "NA";
				}
				print $default_handle "\t";
				if ($cat eq 'fsm') {
					print $default_handle "NA";
				} elsif (exists ($longest_nc_valid{$SJ_chain})) {
					print $default_handle "valid($longest_nc_valid{$SJ_chain})";
				} else {
					print $default_handle "non_valid";
				}
				print $default_handle "\t";
				if (exists ($substring_nc{$SJ_chain})) {
					print $default_handle "$_," for @{$substring_nc{$SJ_chain}};
				} else {
					print $default_handle "NA";
				}
				print $default_handle "\t";
				if (exists $unannotated_isoforms{$SJ_chain} and @{$unannotated_isoforms{$SJ_chain}[2]}>0) {
					print $default_handle "$_," for @{$unannotated_isoforms{$SJ_chain}[1]};
				} else {
					print $default_handle "NA";
				}
				print $default_handle "\t";
				if (exists $unannotated_isoforms{$SJ_chain} and @{$unannotated_isoforms{$SJ_chain}[2]}>0) {
					print $default_handle "$_," for @{$unannotated_isoforms{$SJ_chain}[2]};
				} else {
					print $default_handle "NA";
				}
				if (exists $unannotated_isoforms{$SJ_chain}) {
					print $default_handle "\t$unannotated_isoforms{$SJ_chain}[3]\t$unannotated_isoforms{$SJ_chain}[4]";
				} else {
					print $default_handle "\tNA\tNA";
				}

				for my $sample (@{$samples_sort_ref}) {
					print $default_handle "\t";
					print $default_handle "$_," for @{$read_count_isoform{$SJ_chain}{$sample}};
				}
				print $default_handle "\n";
			}
		}
		

		if (!defined $raw) {


		my $m = 0;
		my $largest_diff;
		# The EM algorithm loops until the largest_diff drops below a threshold.
		# The outcome of the largest_diff check can change based on the order of
		# floating point operations. For example the largest_diff could be calculated as
		# 0.09999 instead of 0.1 which would lead to an extra loop iteration.
		# Sort the data structures to ensure a consistent result.
		my @sorted_ism_mt1_chain_count_keys;
		my %sorted_samples_by_ism_mt1_chain;
		while (my ($chain, $sample_ref) = each %ism_mt1_chain_count) {
			push @sorted_ism_mt1_chain_count_keys, $chain;
			my @sorted_samples_for_chain;
			while (my ($sample, $chains_count_ref) = each %{$sample_ref}) {
				push @sorted_samples_for_chain, $sample;
				@{$chains_count_ref} = sort {$a->[0] cmp $b->[0]} @{$chains_count_ref};
			}
			@sorted_samples_for_chain = sort {$a cmp $b} @sorted_samples_for_chain;
			$sorted_samples_by_ism_mt1_chain{$chain} = \@sorted_samples_for_chain;
		}
		@sorted_ism_mt1_chain_count_keys = sort {$a cmp $b} @sorted_ism_mt1_chain_count_keys;

		for my $cat ('nic', 'nnc') {
			for my $chain (@{$cat_SJ_chains_sort{$cat}}) {
				if (defined $SJ_chain_cat{$cat}{$chain}[1] and scalar(@{$SJ_chain_cat{$cat}{$chain}[1]}) > 1) {
					@{$SJ_chain_cat{$cat}{$chain}[1]} = sort {$a cmp $b} @{$SJ_chain_cat{$cat}{$chain}[1]};
				}
			}
		}

		my @sorted_possible_fsm_num;
		my %sorted_chains_by_possible_fsm_num;
		while (my ($possible_fsm_num, $isoforms_by_chain_ref) = each %possible_fsm_filtered) {
			push @sorted_possible_fsm_num, $possible_fsm_num;
			my @sorted_chains_for_fsm_num;
			while (my ($filtered_SJ_chain, $fsm_chain_ref) = each %{$isoforms_by_chain_ref}) {
				push @sorted_chains_for_fsm_num, $filtered_SJ_chain;
				@{$fsm_chain_ref} = sort {$a cmp $b} @{$fsm_chain_ref};
			}
			@sorted_chains_for_fsm_num = sort {$a cmp $b} @sorted_chains_for_fsm_num;
			$sorted_chains_by_possible_fsm_num{$possible_fsm_num} = \@sorted_chains_for_fsm_num;
		}
		@sorted_possible_fsm_num = sort {$a <=> $b} @sorted_possible_fsm_num;
		while (1){
			$m ++;
			$largest_diff = 0;

			for my $ism_SJ_chain (@sorted_ism_mt1_chain_count_keys) {
				my $sample_ref = $ism_mt1_chain_count{$ism_SJ_chain};
				my $total_count = 0;
				for my $sample (@{$sorted_samples_by_ism_mt1_chain{$ism_SJ_chain}}) {
					my $chains_count_ref = $sample_ref->{$sample};
					for my $chains_string_freq_ref (@{$chains_count_ref}) {
						my ($chains_string, $count) = @{$chains_string_freq_ref};
						$total_count += $count;
						my @chains = split ':', $chains_string;
						my $total_read_count_fsm = 0;
						$total_read_count_fsm += $read_count_isoform{$_}{$sample}[2] for @chains;
						if ($total_read_count_fsm > 0) {
							$read_count_isoform{$_}{$sample}[3] += sprintf("%.2f",$count/$total_read_count_fsm*$read_count_isoform{$_}{$sample}[2]) for @chains;
						} else {
							$read_count_isoform{$_}{$sample}[3] += sprintf("%.2f",$count/@chains) for @chains;
						}
						if ($should_create_default_handle) {
							print $default_handle "special_ISM_read_assigned:$ism_SJ_chain\t$sample\t",scalar(@chains), "\t";
							print $default_handle "$_($count)" for @chains;
							print $default_handle "\n";
						}
					}
					
				}
				printf $default_handle "read_count_summary_in_EM\t$chr\tISM\t$n:$ism_SJ_chain\tEM\t%g\t$m\n", $total_count if $should_create_default_handle;
			}

			for my $cat ('nic', 'nnc') {
				for my $chain (@{$cat_SJ_chains_sort{$cat}}) {
					if (defined $SJ_chain_cat{$cat}{$chain}[1] and scalar(@{$SJ_chain_cat{$cat}{$chain}[1]}) > 1) {
						my @chains_nc = map {$nc_valid_longest[$_]} @{$SJ_chain_cat{$cat}{$chain}[1]};
						my $total_count = 0;
						for my $sample (@{$samples_sort_ref}) {
							next if $read_count_isoform{$chain}{$sample}[0] == 0;
							$total_count += $read_count_isoform{$chain}{$sample}[0];
							my $total_read_count_nc = 0;
							$total_read_count_nc += $read_count_isoform{$_}{$sample}[2] for @chains_nc;
							if ($total_read_count_nc > 0) {
								$read_count_isoform{$_}{$sample}[3] += sprintf("%.2f",$read_count_isoform{$chain}{$sample}[0]/$total_read_count_nc*$read_count_isoform{$_}{$sample}[2]) for @chains_nc;
							} else {
								$read_count_isoform{$_}{$sample}[3] += sprintf("%.2f",$read_count_isoform{$chain}{$sample}[0]/@chains_nc) for @chains_nc;
							}
						}
						printf $default_handle "read_count_summary_in_EM\t$chr\t$cat\t$n:$chain\tEM\t%g\t$m\n", $total_count if $should_create_default_handle;
					}
				}
			}
			for my $possible_fsm_num(@sorted_possible_fsm_num) {
				if ($possible_fsm_num>1){
					my $isoforms_by_chain_ref = $possible_fsm_filtered{$possible_fsm_num};
					for my $filtered_SJ_chain (@{$sorted_chains_by_possible_fsm_num{$possible_fsm_num}}) {
						my $fsm_chain_ref = $isoforms_by_chain_ref->{$filtered_SJ_chain};
						my $total_count = 0;
						for my $sample (@{$samples_sort_ref}) {
							next if $read_count_isoform{$filtered_SJ_chain}{$sample}[0] == 0;
							$total_count += $read_count_isoform{$filtered_SJ_chain}{$sample}[0];
							my $total_read_count_fsm = 0;
							$total_read_count_fsm += $read_count_isoform{$_}{$sample}[2] for @{$fsm_chain_ref};
							if ($total_read_count_fsm > 0) {
								$read_count_isoform{$_}{$sample}[3] += sprintf("%.2f",$read_count_isoform{$filtered_SJ_chain}{$sample}[0]/$total_read_count_fsm*$read_count_isoform{$_}{$sample}[2]) for @{$fsm_chain_ref};
							} else {
								$read_count_isoform{$_}{$sample}[3] += sprintf("%.2f",$read_count_isoform{$filtered_SJ_chain}{$sample}[0]/@{$fsm_chain_ref}) for @{$fsm_chain_ref};
							}
						}
						printf $default_handle "read_count_summary_in_EM\t$chr\tfiltered\t$n:$filtered_SJ_chain\tEM\t%g\t$m\n", $total_count if $should_create_default_handle;
					}
					
				}
			}
			for my $chain (keys %read_count_isoform){
				for my $sample (@{$samples_sort_ref}) {
					$read_count_isoform{$chain}{$sample}[2] = $read_count_isoform{$chain}{$sample}[1] + $read_count_isoform{$chain}{$sample}[3];
					$largest_diff = abs($read_count_isoform{$chain}{$sample}[3] - $read_count_isoform{$chain}{$sample}[4]) if $largest_diff < abs($read_count_isoform{$chain}{$sample}[3] - $read_count_isoform{$chain}{$sample}[4]);
					$read_count_isoform{$chain}{$sample}[4] = $read_count_isoform{$chain}{$sample}[3];
					$read_count_isoform{$chain}{$sample}[3] = 0;
				}
			}

			if ($largest_diff < .1){
				print $default_handle "converge: $n\t$largest_diff\t$m\n" if $should_create_default_handle;
				last;
			} elsif ($m>=$max_iterate){
				print $default_handle "not_converge: $n\t$largest_diff\t$m\n" if $should_create_default_handle;
				last;
			}
		}
		@sorted_ism_mt1_chain_count_keys = ();
		%sorted_samples_by_ism_mt1_chain = ();
		@sorted_possible_fsm_num = ();
		%sorted_chains_by_possible_fsm_num = ();

		while(my ($fsm_SJ_chain, $chain_info_ref) = each %{$SJ_chain_cat{'fsm'}}){
			my $isSameStrand = $SJ_chain_read{$fsm_SJ_chain}[3];
			my @SJ_IDs = split '_', $fsm_SJ_chain;
			my @SJs = map {$all_SJ_group_sort[$_]} @SJ_IDs;
			if ($should_create_default_handle) {
				print $default_handle "assign_reads\t$fsm_SJ_chain\t";
				print $default_handle "$read_count_isoform{$fsm_SJ_chain}{$_}[0]," for @{$samples_sort_ref};
				print $default_handle "\t";
				printf $default_handle "%.2f,", $read_count_isoform{$fsm_SJ_chain}{$_}[1] for @{$samples_sort_ref};
				print $default_handle "\t";
				printf $default_handle "%.2f,", $read_count_isoform{$fsm_SJ_chain}{$_}[2] for @{$samples_sort_ref};
				print $default_handle "\t";
				
				print $default_handle "$_," for @{${$chain_info_ref}[-1]};
				print $default_handle "\t";
				print $default_handle "$SJ_chain_read{$fsm_SJ_chain}[1],";
				print $default_handle "$_," for @SJs;
				print $default_handle "$SJ_chain_read{$fsm_SJ_chain}[2]";
				print $default_handle "\n";
			}
			
			#my @SJ0_info = split ':', $SJs[0];

			my (%possible_genes, %possible_isoform_names, %total_perfect_match);
			for my $isoform (@{${$chain_info_ref}[-1]}) {
				for my $sample (@{$samples_sort_ref}){
					$total_perfect_match{$sample} += @{$isoform_perfect_match{$isoform}{$sample}} if exists $isoform_perfect_match{$isoform} and exists $isoform_perfect_match{$isoform}{$sample};
				}
			}

			for my $isoform ( @{${$chain_info_ref}[-1]} ){
				my $abu_output = "$isoform\t$isoform_info_ref->{$isoform}[-1]\t$isoform_info_ref->{$isoform}[0]";
				my $tag = 0;
				for my $sample (@{$samples_sort_ref}){
					if ($total_perfect_match{$sample} > 0){
						if ( exists $isoform_perfect_match{$isoform} and exists $isoform_perfect_match{$isoform}{$sample} ) {
							$abu_output .= sprintf ("\t%.2f", $read_count_isoform{$fsm_SJ_chain}{$sample}[2]*@{$isoform_perfect_match{$isoform}{$sample}}/$total_perfect_match{$sample});
							$tag = 1;
						} else {
							$abu_output .= "\t0";
						}
					} else {
						if ($read_count_isoform{$fsm_SJ_chain}{$sample}[2]>0) {
							$abu_output .= sprintf ("\t%.2f", $read_count_isoform{$fsm_SJ_chain}{$sample}[2]/@{${$chain_info_ref}[-1]});
							$tag = 1;
						} else {
							$abu_output .= "\t0";
						}
					}
				}
				if ($tag == 1){
					print $abu_handle $abu_output."\n";
					&gtf_output($gtf_handle, $chr, $isoform, [$multi_exon_isoform_end_ref->{$isoform}{'0'},$multi_exon_isoform_end_ref->{$isoform}{'1'}], \@SJs, $strand_symbol_ref->{$isSameStrand}, 'annotated_isoform');
				}
				
			}
			if ($should_create_tsv_compt_handle){
				for my $read_info_ref (@{$SJ_chain_read{$fsm_SJ_chain}[-1]}) {
					my $read_ID = ${$read_info_ref}[0];
					print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tFSM\t";
					print $tsv_compt_handle "$_," for @{${$chain_info_ref}[-1]};
					print $tsv_compt_handle "\n";
				}
			}
			printf $default_handle "read_count_summary\t$chr\tfsm\t$n:$fsm_SJ_chain\titself\t%g\n", scalar(@{$SJ_chain_read{$fsm_SJ_chain}[-1]}) if $should_create_default_handle;
		}
		while(my ($nc_SJ_chain, $nc_SJ_chain_ID) = each %longest_nc_valid){
			my $isSameStrand = $SJ_chain_read{$nc_SJ_chain}[3];
			my @SJ_IDs = split '_', $nc_SJ_chain;
			my @SJs = map {$all_SJ_group_sort[$_]} @SJ_IDs;
			if ($should_create_default_handle) {
				print $default_handle "assign_reads\t$nc_SJ_chain\t";
				print $default_handle "$read_count_isoform{$nc_SJ_chain}{$_}[0]," for @{$samples_sort_ref};
				print $default_handle "\t";
				printf $default_handle "%.2f,", $read_count_isoform{$nc_SJ_chain}{$_}[1] for @{$samples_sort_ref};
				print $default_handle "\t";
				printf $default_handle "%.2f,", $read_count_isoform{$nc_SJ_chain}{$_}[2] for @{$samples_sort_ref};
				print $default_handle "\t";
				
				if (exists $unannotated_isoforms{$nc_SJ_chain} and @{$unannotated_isoforms{$nc_SJ_chain}[1]} > 0){
					print $default_handle "$_," for @{$unannotated_isoforms{$nc_SJ_chain}[1]};
				} else {
					print $default_handle "$n:$nc_SJ_chain_ID\t";
				}
				if (exists $unannotated_isoforms{$nc_SJ_chain}){
					print $default_handle "$unannotated_isoforms{$nc_SJ_chain}[3]\t$unannotated_isoforms{$nc_SJ_chain}[4]\t";
				} else {
					print $default_handle "0\t0\t";
				}
				print $default_handle "$SJ_chain_read{$nc_SJ_chain}[1],";
				print $default_handle "$_," for @SJs;
				print $default_handle "$SJ_chain_read{$nc_SJ_chain}[2]";
				print $default_handle "\n";
			}
			#my @SJ0_info = split ':', $SJs[0];
			

			if (exists $unannotated_isoforms{$nc_SJ_chain} and @{$unannotated_isoforms{$nc_SJ_chain}[1]} > 0){
				my $isoform_ID = join '/', @{$unannotated_isoforms{$nc_SJ_chain}[1]};

				my (%possible_genes, %possible_isoform_names);
				for my $isoform (@{$unannotated_isoforms{$nc_SJ_chain}[1]}) {
					$possible_genes{$isoform_info_ref->{$isoform}[0]} ++;
					$possible_isoform_names{$isoform_info_ref->{$isoform}[-1]} ++;
				}
				my @possible_genes_sort = sort {$possible_genes{$b} <=> $possible_genes{$a}} keys %possible_genes;
				my @possible_isoform_names_sort = sort {$possible_isoform_names{$b} <=> $possible_isoform_names{$a}} keys %possible_isoform_names;
				my $possible_genes_ID = join ',', @possible_genes_sort;
				my $possible_isoform_name = join ',', @possible_isoform_names_sort;
				print $abu_handle "$isoform_ID\t$possible_isoform_name\t$possible_genes_ID";

				printf $abu_handle "\t%.2f", $read_count_isoform{$nc_SJ_chain}{$_}[2] for @{$samples_sort_ref};
				print $abu_handle "\n";
				&gtf_output($gtf_handle, $chr, $isoform_ID, [$SJ_chain_read{$nc_SJ_chain}[1],$SJ_chain_read{$nc_SJ_chain}[2]], \@SJs, $strand_symbol_ref->{$isSameStrand}, 'annotated_isoform');
			} else {
				my $isoform_ID = "ESPRESSO:$chr:$n:$nc_SJ_chain_ID";
				print $abu_handle "$isoform_ID\tNA";
				my %possible_genes;
				for my $SJ(@SJs) {
					if (exists $anno_SJ_ref->{$SJ}){
						$possible_genes{$isoform_info_ref->{$_}[0]} ++ for @{$anno_SJ_ref->{$SJ}[-1]};
					} elsif (exists $anno_SJ_complete_ref->{$SJ}){
						$possible_genes{$isoform_info_ref->{$_}[0]} ++ for @{$anno_SJ_complete_ref->{$SJ}[-1]};
					}
				}
				my @possible_genes_sort = sort {$possible_genes{$b} <=> $possible_genes{$a}} keys %possible_genes;
				if (@possible_genes_sort > 0) {
					print $abu_handle "\t$possible_genes_sort[0]";
					print $abu_handle ",$possible_genes_sort[$_]" for 1 .. $#possible_genes_sort;
				} else {
					print $abu_handle "\tNA";
				}
				
				printf $abu_handle "\t%.2f", $read_count_isoform{$nc_SJ_chain}{$_}[2] for @{$samples_sort_ref};
				print $abu_handle "\n";
				&gtf_output($gtf_handle, $chr, "$isoform_ID", [$SJ_chain_read{$nc_SJ_chain}[1],$SJ_chain_read{$nc_SJ_chain}[2]], \@SJs, $strand_symbol_ref->{$isSameStrand}, 'novel_isoform');
			}
			
		}

		my %single_exon_isoform_read_count;
		my @single_exon_isoforms_group_sort = sort {
			my $a_info = $single_exon_isoform_end_ref->{$a};
			my $a_start = $a_info->{'0'};
			my $a_end = $a_info->{'1'};
			my $a_length = $a_end - $a_start;
			my $b_info = $single_exon_isoform_end_ref->{$b};
			my $b_start = $b_info->{'0'};
			my $b_end = $b_info->{'1'};
			my $b_length = $b_end - $b_start;
			$a_length <=> $b_length;
		} @{$single_exon_isoforms_by_group{$n}};

		for my $single_exon_isoform(@single_exon_isoforms_group_sort){
			$single_exon_isoform_read_count{$single_exon_isoform}{$_} = [0,0,0] for @{$samples_sort_ref}; #direct_assign, last_EM_assign, current_EM_assign
		}
		my %possible_single_exon_isoform;
		my ($total_count1, $total_count2) = (0,0);
		for my $read_ID (@single_exon_reads){
			my $sample = $read_info{$read_ID}{'sample'};
			if (exists $read_info{$read_ID}{'start'} and exists $read_info{$read_ID}{'end'}) {
				for my $single_exon_isoform(@single_exon_isoforms_group_sort){
					if ( $single_exon_isoform_end_ref->{$single_exon_isoform}{'0'}-$SJ_dist<=$read_info{$read_ID}{'start'} and abs($single_exon_isoform_end_ref->{$single_exon_isoform}{'1'}+$SJ_dist >=$read_info{$read_ID}{'end'}) ){
						push @{$possible_single_exon_isoform{$read_ID}}, $single_exon_isoform;
					}
				}
			}
			if(exists $possible_single_exon_isoform{$read_ID}) {
				$total_count1 ++;
				print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tsingle-exon\t" if $should_create_tsv_compt_handle;
				for my $single_exon_isoform (@{$possible_single_exon_isoform{$read_ID}}) {
					$single_exon_isoform_read_count{$single_exon_isoform}{$sample}[0] += sprintf("%.2f", 1/@{$possible_single_exon_isoform{$read_ID}});
					print $tsv_compt_handle "$single_exon_isoform," if $should_create_tsv_compt_handle;
				}
				print $tsv_compt_handle "\n" if $should_create_tsv_compt_handle;
			} else {
				$total_count2 ++;
				print $tsv_compt_handle "$read_info{$read_ID}{'ori_ID'}\t$read_info{$read_ID}{'sample'}\tsingle-exon\tNA\n" if $should_create_tsv_compt_handle;
				
			}
		}
		if ($should_create_default_handle) {
			printf $default_handle "read_count_summary\t$chr\tone_exon\t$n:total\tEM\t%g\n", $total_count1;
			printf $default_handle "read_count_summary\t$chr\tone_exon\t$n:total\tNA\t%g\n", $total_count2;
			printf $default_handle "read_count_summary_one_exon\t$chr\t$n\t%g\n", scalar(@single_exon_reads);
		}
		
		$m = 0;
		for my $single_exon_isoform(@single_exon_isoforms_group_sort){
			$single_exon_isoform_read_count{$single_exon_isoform}{$_}[1] = $single_exon_isoform_read_count{$single_exon_isoform}{$_}[0] for @{$samples_sort_ref};
		}
		# Sort to ensure consistent result of EM algorithm
		my @sorted_possible_single_exon_isoform_keys;
		while (my ($read_ID, $isoforms_ref) = each %possible_single_exon_isoform) {
			push @sorted_possible_single_exon_isoform_keys, $read_ID;
			@{$isoforms_ref} = sort {$a cmp $b} @{$isoforms_ref};
		}
		@sorted_possible_single_exon_isoform_keys = sort {$a cmp $b} @sorted_possible_single_exon_isoform_keys;
		while (1){
			$m ++;
			$largest_diff = 0;
			for my $read_ID (@sorted_possible_single_exon_isoform_keys) {
				my $isoforms_ref = $possible_single_exon_isoform{$read_ID};
				my $sample = $read_info{$read_ID}{'sample'};
				my $total_read_count_isoforms= 0;
				$total_read_count_isoforms += $single_exon_isoform_read_count{$_}{$sample}[1] for @{$isoforms_ref};
				if ($total_read_count_isoforms > 0) {
					for my $single_exon_isoform (@{$isoforms_ref}) {
						$single_exon_isoform_read_count{$single_exon_isoform}{$sample}[2] += sprintf("%.2f", 1/$total_read_count_isoforms*$single_exon_isoform_read_count{$single_exon_isoform}{$sample}[1]);
					}
				} else {
					for my $single_exon_isoform (@{$isoforms_ref}) {
						$single_exon_isoform_read_count{$single_exon_isoform}{$sample}[2] += sprintf("%.2f", 1/@{$isoforms_ref});
					}
				}
			}
			printf $default_handle "read_count_summary_in_EM\t$chr\tone_exon\t$n:total\tEM\t%g\t$m\n", scalar(keys %possible_single_exon_isoform) if $should_create_default_handle;
			for my $single_exon_isoform (@single_exon_isoforms_group_sort){
				for my $sample (@{$samples_sort_ref}) {
					$largest_diff = abs($single_exon_isoform_read_count{$single_exon_isoform}{$sample}[2] - $single_exon_isoform_read_count{$single_exon_isoform}{$sample}[1]) if $largest_diff < abs($single_exon_isoform_read_count{$single_exon_isoform}{$sample}[2] - $single_exon_isoform_read_count{$single_exon_isoform}{$sample}[1]);
					$single_exon_isoform_read_count{$single_exon_isoform}{$sample}[1] = $single_exon_isoform_read_count{$single_exon_isoform}{$sample}[2];
					$single_exon_isoform_read_count{$single_exon_isoform}{$sample}[2] = 0;
				}
				
			}


			if ($largest_diff < .1){
				print $default_handle "single_converge: $n\t$largest_diff\t$m\n" if $should_create_default_handle;
				last;
			} elsif ($m>=$max_iterate){
				print $default_handle "single_not_converge: $n\t$largest_diff\t$m\n" if $should_create_default_handle;
				last;
			}
		}
		@sorted_possible_single_exon_isoform_keys = ();

		while (my ($single_exon_isoform, $sample_ref) = each %single_exon_isoform_read_count){
			my $total_read_count_EM = 0;
			$total_read_count_EM += ${$sample_ref}{$_}[1] for @{$samples_sort_ref};
			if($total_read_count_EM > 0) {				

				if ($should_create_default_handle) {
					print $default_handle "assign_reads\tNA\tNA\t";
					print $default_handle "${$sample_ref}{$_}[0]," for @{$samples_sort_ref};
					print $default_handle "\t";
					printf $default_handle "%.2f,", ${$sample_ref}{$_}[1] for @{$samples_sort_ref};
					print $default_handle "\t";
					print $default_handle "$single_exon_isoform,\t";
					print $default_handle "$single_exon_isoform_end_ref->{$single_exon_isoform}{'0'},";
					print $default_handle "$single_exon_isoform_end_ref->{$single_exon_isoform}{'1'}";
					print $default_handle "\n";
				}
				$single_exon_isoform_total{$single_exon_isoform}{$_} += ${$sample_ref}{$_}[1] for @{$samples_sort_ref};
			}
		}



		}
	}

	while (my ($single_exon_isoform, $sample_ref) = each %single_exon_isoform_total){
		print $abu_handle "$single_exon_isoform\t$isoform_info_ref->{$single_exon_isoform}[-1]\t$isoform_info_ref->{$single_exon_isoform}[0]";
		printf $abu_handle "\t%.2f", ${$sample_ref}{$_} for @{$samples_sort_ref};
		print $abu_handle "\n";
		&gtf_output($gtf_handle, $chr, $single_exon_isoform, [$single_exon_isoform_end_ref->{$single_exon_isoform}{'0'},$single_exon_isoform_end_ref->{$single_exon_isoform}{'1'}], [], $single_exon_isoform_end_ref->{$single_exon_isoform}{'strand'}, 'annotated_isoform');

	}

	if ($should_create_default_handle) {
		my %total_type;
		while (my ($chr, $type_ref) = each %chr_type) {
			print $default_handle "$chr";
			for my $type ("fsm", "ism", "nic", "nnc") {
				if (exists $chr_type{$chr}{$type}) {
					print $default_handle " $chr_type{$chr}{$type}";
					$total_type{$type} += $chr_type{$chr}{$type};
				} else {
					print $default_handle " NA";
				}
			}
			print $default_handle "\n";
		}
		for my $type ("fsm", "ism", "nic", "nnc") {
			if (exists $total_type{$type}) {
				print $default_handle " $total_type{$type}";
			} else {
				print $default_handle " NA";
			}
		}
		print $default_handle "\n";
	}

	close $default_handle if $should_create_default_handle;
	close $tsv_compt_handle if $should_create_tsv_compt_handle;
	close $abu_handle;
	close $gtf_handle;
}

sub gtf_output {
	my ($gtf_handle, $chr, $isoform_ID, $ends_ref, $SJs_ref, $strand, $type) = @_;
	print $gtf_handle "$chr\t$type\ttranscript\t",${$ends_ref}[0]+1,"\t${$ends_ref}[1]\t.\t$strand\t.\ttranscript_id \"$isoform_ID\";\n";
	if (@{$SJs_ref} > 0) {
		my @exons = ([${$ends_ref}[0]+1]);
		for my $SJ(@{$SJs_ref}) {
			my @SJ_info = split ':', $SJ;
			push @{$exons[-1]}, $SJ_info[1];
			push @exons, [$SJ_info[2]+1];
		}
		push @{$exons[-1]}, ${$ends_ref}[1];
		for my $exon_ID (0 .. $#exons) {
			printf $gtf_handle "$chr\t$type\texon\t${$exons[$exon_ID]}[0]\t${$exons[$exon_ID]}[1]\t.\t$strand\t.\ttranscript_id \"$isoform_ID\"; exon_number \"%g\";\n", $exon_ID+1;
		}
	} else {
		print $gtf_handle "$chr\t$type\texon\t",${$ends_ref}[0]+1,"\t${$ends_ref}[1]\t.\t$strand\t.\ttranscript_id \"$isoform_ID\"; exon_number \"1\";\n";
	}
}

sub valid_nc {
	my ($default_handle, $chain, $chr_SJ_reads_ref, $SJ_chain_read_ref, $read_info_ref, $all_SJ_group_sort_ref, $chr_anno_SJ_ref, $SJ_isoform_SJ_ref, $read_num_cutoff, $read_ratio_cutoff) = @_;
	my @SJ_IDs = split '_', $chain;
	print $default_handle "validation_for_nc_chain: $chain" if defined $default_handle;

	my (@SJs, @SJs_valid);
	my (%possible_isoforms, @fsm_unanno, @ism_unanno, %SJ_info_nc);
	for my $SJ_ID (@SJ_IDs) {
		my $SJ = ${$all_SJ_group_sort_ref}[$SJ_ID];
		my ($chr, $start, $end, $strand) = split ':', $SJ;
		$SJ_info_nc{$SJ} = [$start, $end, $strand];
		if(exists ${$chr_anno_SJ_ref}{$SJ}) {
			$possible_isoforms{$_}++ for @{${$chr_anno_SJ_ref}{$SJ}[-1]};
		}
		push @SJs, $SJ;
		print $default_handle "\t$SJ" if defined $default_handle;
		if ( !exists ${$chr_SJ_reads_ref}{$SJ} ) {
			die "0:$chain\t$SJ_ID\t$SJ";
		} 
	}
	print $default_handle "\n" if defined $default_handle;
	
	for my $possible_fsm_complete (keys %possible_isoforms) {
		my $comp_result = &comp_minimap_reference_SJ(\%SJ_info_nc, ${$SJ_isoform_SJ_ref}{$possible_fsm_complete}, $default_handle);
		if ($comp_result == 1) {
			push @fsm_unanno, $possible_fsm_complete;
		} elsif ($comp_result == 2) {
			push @ism_unanno, $possible_fsm_complete;
		}
	}
	if(@fsm_unanno >= 1) {
		if (defined $default_handle) {
			print $default_handle "unannotated_isoforms $chain: ";
			print $default_handle "$_," for @fsm_unanno;
			print $default_handle "\n";
		}
	} elsif(@ism_unanno >= 1) {
		if (defined $default_handle) {
			print $default_handle "unannotated_isoforms_fragment $chain: ";
			print $default_handle "$_," for @ism_unanno;
			print $default_handle "\n";
		}
	}
	
	my @reads_info = @{$SJ_chain_read_ref};
	die "no reads loaded for $chain" if scalar(@reads_info) == 0;
	my @valid_read_SJs;
	for my $SJ (@SJs) {
		my ($chr,$start,$end,$strand) = split ':', $SJ;
		
			my @valid_reads = grep {${$chr_SJ_reads_ref}{$SJ}{${$read_info_ref}{${$_}[0]}{'ori_ID'}} >= 1} @reads_info;
			if (@valid_reads >= $read_num_cutoff and @valid_reads/@reads_info >= $read_ratio_cutoff) {
				push @SJs_valid, $SJ;
			}
			printf $default_handle "\t$SJ(%g/%g)", scalar(@valid_reads), scalar(@reads_info) if defined $default_handle;
			push @valid_read_SJs, scalar(@valid_reads);

	}
	print $default_handle "\n" if defined $default_handle;
	my @valid_read_SJs_sort = sort {$a <=> $b} @valid_read_SJs;
	if (@SJs_valid == @SJs) {
		print $default_handle "yes_validated: $chain\n" if defined $default_handle;
		[1,\@fsm_unanno,\@ism_unanno, $valid_read_SJs_sort[0], scalar(@reads_info)];
	} else {
		print $default_handle "not_validated: $chain\n" if defined $default_handle;
		[0,\@fsm_unanno,\@ism_unanno, $valid_read_SJs_sort[0], scalar(@reads_info)];
	}
}

sub comp_SJ_chain{
	my ($SJ_chain_short, $SJ_chain_long) = @_;
	my @SJs_short = split '_', $SJ_chain_short;
	my @SJs_long = split '_', $SJ_chain_long;
	if (@SJs_short > @SJs_long or @SJs_short == 0 or @SJs_long == 0) {
		return 0;
	}
	my @non_x_SJ_IDs = grep {$SJs_short[$_] ne 'x'} (0 .. $#SJs_short);
	my ($matched_1st_non_x_SJ_ID, $match_SJ_num);

	for my $XJ_ID_long ($non_x_SJ_IDs[0] .. $#SJs_long) {
		if (defined $matched_1st_non_x_SJ_ID and defined $SJs_short[$XJ_ID_long-$matched_1st_non_x_SJ_ID+$non_x_SJ_IDs[0]] and $SJs_short[$XJ_ID_long-$matched_1st_non_x_SJ_ID+$non_x_SJ_IDs[0]] eq $SJs_long[$XJ_ID_long]) {
			$match_SJ_num ++;
		} elsif ($SJs_short[$non_x_SJ_IDs[0]] eq $SJs_long[$XJ_ID_long]) {
			$matched_1st_non_x_SJ_ID = $XJ_ID_long;
			if ($#SJs_long-$matched_1st_non_x_SJ_ID < $#SJs_short-$non_x_SJ_IDs[0]) {
				return 0;
			}
			$match_SJ_num ++;
		}
	}
	if ($match_SJ_num == @SJs_long){
		1;
	} elsif ($match_SJ_num == @SJs_short){
		2;
	} elsif ($match_SJ_num == @non_x_SJ_IDs){
		3;
	} else {
		0;
	}
}

sub comp_minimap_reference_SJ{
	my ($minimap_SJ_ref, $ref_SJ_ref, $default_handle) = @_;
	my @ref_SJ_sort = sort {${$ref_SJ_ref}{$a}[0] <=> ${$ref_SJ_ref}{$b}[0] or ${$ref_SJ_ref}{$a}[1] <=> ${$ref_SJ_ref}{$b}[1]} keys %{$ref_SJ_ref};
	my @minimap_SJ_sort = sort {${$minimap_SJ_ref}{$a}[1] <=> ${$minimap_SJ_ref}{$b}[1] or ${$minimap_SJ_ref}{$a}[2] <=> ${$minimap_SJ_ref}{$b}[2]} keys %{$minimap_SJ_ref};
	my ($ref_SJ_sort_string, $minimap_SJ_sort_string) = (join('_',@ref_SJ_sort), join('_',@minimap_SJ_sort) );
	my $result = &comp_SJ_chain($minimap_SJ_sort_string, $ref_SJ_sort_string);
	print $default_handle "$minimap_SJ_sort_string\t$ref_SJ_sort_string\t$result\n" if defined $default_handle;
	return $result;
}

sub gtf_path_for_chr {
	my ($chr, $out_dir) = @_;
	return $out_dir."/$chr.updated.gtf.tmp";
}

sub abu_path_for_chr {
	my ($chr, $out_dir) = @_;
	return $out_dir."/$chr.abundance.esp.tmp";
}

sub tsv_compt_path_for_chr {
	my ($chr, $out_dir) = @_;
	return $out_dir."/$chr.tsv_compt.tmp";
}

sub default_path_for_chr {
	my ($chr, $out_dir) = @_;
	return $out_dir."/$chr.default_output.tmp";
}

# package must return true
1;

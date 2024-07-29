package ESPRESSO_Q_Thread;

use strict;
use warnings;

use Fcntl qw(SEEK_SET);
use Storable qw(store retrieve);

if ($0 =~ /ESPRESSO_Q_Thread/) {
  # Run as a script
  my $stored_shared_arguments = $ARGV[0];
  my $stored_chr_arguments = $ARGV[1];
  main($stored_shared_arguments, $stored_chr_arguments);
} else {
  # Just for imports
}

sub split_sj_string {
  my ($sj_string) = @_;

  my ($chr, $start, $end) = split(':', $sj_string);
  my %result = ();
  $result{'chr'} = $chr;
  $result{'start'} = $start;
  $result{'end'} = $end;
  return \%result;
}

sub sort_by_numeric_first_column {
  my ($unsorted_path, $sorted_path, $buffer_size) = @_;

  my $sort_command = (
    "sort $buffer_size -k1,1n $unsorted_path --output $sorted_path");
  my $exit_code = system($sort_command);
  if ($exit_code != 0) {
    die "'$sort_command' exited with code $exit_code";
  }
}

sub write_tsv_columns {
  my ($handle, $columns_ref) = @_;

  my $with_tabs = join("\t", @{$columns_ref});
  print $handle "$with_tabs\n";
}

sub read_tsv_columns {
  my ($line, $columns_ref) = @_;

  $line =~ s/\r\n//;
  chomp($line);
  @{$columns_ref} = split(/\t/, $line);
}

# Trailing commas are used to be consistent with previous versions
sub join_with_trailing {
  my ($sep, $values_ref) = @_;

  my $result = join($sep, @{$values_ref});
  $result = $result . $sep;
  return $result;
}

sub update_group_ends {
  my ($group, $start, $end, $group_ends_ref) = @_;

  if (!exists $group_ends_ref->{$group}) {
    $group_ends_ref->{$group} = [$start, $end];
    return;
  }

  my $ends_ref = $group_ends_ref->{$group};
  if ($start < $ends_ref->[0]) {
    $ends_ref->[0] = $start;
  }
  if ($end > $ends_ref->[1]) {
    $ends_ref->[1] = $end;
  }
}

sub load_read_final {
  my ($out_dir, $chr, $keep_tmp, $target_col_index, $sort_buffer_size,
      $group_info_path, $read_final_paths_ref, $anno_SJ_ref, $unanno_SJ_ref,
      $group_ends_ref, $group_sort_ref, $sj_to_groups_ref,
      $summary_data_ref) = @_;

  my $use_corrected = ($target_col_index == 7);
  my $num_filtered = 0;
  my %read_info = ();

  open(my $group_info_handle, '>', $group_info_path)
    or die "cannot open $group_info_path: $!";
  for my $read_final_path_info (@{$read_final_paths_ref}) {
    my $read_final_path = $read_final_path_info->{'path'};
    my $input_dir = $read_final_path_info->{'input_dir'};

    open(my $read_final_handle, '<', $read_final_path)
      or die "cannot open $read_final_path: $!";
    while (<$read_final_handle>) {
      s/\r\n//;
      chomp;

      my @line = ();
      read_tsv_columns($_, \@line);

      my $read_id = $line[0];
      my $feature = $line[1];
      if ($feature eq 'group_ID') {
        # group_ID will be the first line (feature) for each read.
        # All lines for a previous read have been read.
        if (scalar(keys %read_info) != 0) {
          my $group = $read_info{'group'};
          my $start = $read_info{'start'};
          my $end = $read_info{'end'};
          update_group_ends($group, $start, $end, $group_ends_ref);
          write_tsv_columns($group_info_handle,
                            [$group, 'read', $read_info{'id'},
                             $read_info{'ori_id'}, $read_info{'sample'},
                             $read_info{'strand'}, $start, $end]);
        }
        %read_info = ();

        my $group_id = $line[2];
        my $updated_group_id = $line[3];
        my $sample = $line[4];
        my $read_input_ID = $input_dir . '_' . $updated_group_id;
        $read_info{'group'} = $group_id;
        $read_info{'id'} = $read_input_ID;
        $read_info{'ori_id'} = $read_id;
        $read_info{'sample'} = $sample;
      } elsif ($feature eq 'strand_isoform') {
        my $strand_isoform = $line[2];
        if (!defined $strand_isoform) {
          # TODO seems that this should be '0', '1', or 'unknown'
          # but some reads have no value
          $strand_isoform = '';
        }
        $read_info{'strand'} = $strand_isoform;
      } elsif ($feature eq 'SJ') {
        my $read_input_ID = $read_info{'id'};
        my $group = $read_info{'group'};

        my $sj_string_with_strand = $line[$target_col_index];
        my $sj_number = $line[2];
        if ($sj_string_with_strand eq 'NA') {
          if ($use_corrected) {
            $num_filtered ++;
            write_tsv_columns($group_info_handle,
                              [$group, 'order', $read_input_ID,
                               $sj_number, 'x']);
          }
          next;
        }
        my ($chr, $start, $end, $strand) = split(':', $sj_string_with_strand);
        my $sj_string = join(':', ($chr, $start, $end));
        my $is_annotated = 0;
        write_tsv_columns($group_info_handle, [$group, 'order', $read_input_ID,
                                               $sj_number, $sj_string]);
        if (exists $anno_SJ_ref->{$sj_string}) {
          $is_annotated = 1;
        } else {
          my %details = ();
          $details{'strand'} = $strand;
          $details{'start'} = $start;
          $details{'end'} = $end;
          $unanno_SJ_ref->{$sj_string} = \%details;
        }
        write_tsv_columns($group_info_handle,
                          [$group, 'sj', $read_input_ID, $read_info{'ori_id'},
                           $sj_string, $read_info{'strand'}, $is_annotated]);
        $sj_to_groups_ref->{$sj_string}{$group} = 1;
      } elsif ($feature eq 'start' or $feature eq 'end') {
        my $coordinate = $line[$target_col_index];
        $read_info{$feature} = $coordinate;
      }
    }
    close $read_final_handle;
  }

  if (scalar(keys %read_info) != 0) {
    my $group = $read_info{'group'};
    my $start = $read_info{'start'};
    my $end = $read_info{'end'};
    update_group_ends($group, $start, $end, $group_ends_ref);
    write_tsv_columns($group_info_handle,
                      [$group, 'read', $read_info{'id'}, $read_info{'ori_id'},
                       $read_info{'sample'}, $read_info{'strand'}, $start,
                       $end]);
  }
  %read_info = ();
  close $group_info_handle;
  $summary_data_ref->{'num_reads_with_failed_SJ'} = $num_filtered;

  # The groups originally had non-overlapping coordinates, but
  # since ESPRESSO_C can correct the end points of each read,
  # the coordinates may now overlap.
  # Sort by both the start coordinate and end coordinate.
  @{$group_sort_ref} = sort {
    my $a_info = $group_ends_ref->{$a};
    my $a_start = $a_info->[0];
    my $a_end = $a_info->[1];
    my $b_info = $group_ends_ref->{$b};
    my $b_start = $b_info->[0];
    my $b_end = $b_info->[1];
    ($a_start <=> $b_start) or ($a_end <=> $b_end);
  } keys %{$group_ends_ref};
}

sub load_group_info {
  my ($n, $group_info_path, $group_info_offsets_ref, $SJ_read_improved_ref,
      $SJ_read_order_ref, $all_SJ_reads_ref, $read_info_ref, $read_filtered_ref,
      $group_reads_ref) = @_;

  # %group_info_offsets has the starting offset for groups already seen.
  # If the offset for $n is not known, then use the offset of a lower group,
  # or start from the beginning of the file.
  my $curr_offset = 0;
  if (exists $group_info_offsets_ref->{$n}) {
    $curr_offset = $group_info_offsets_ref->{$n};
  } else {
    while (my ($group, $offset) = each %{$group_info_offsets_ref}) {
      if (($group < $n) and ($offset > $curr_offset)) {
        $curr_offset = $offset;
      }
    }
  }

  open(my $group_info_handle, '<', $group_info_path)
    or die "cannot read $group_info_path: $!";
  seek($group_info_handle, $curr_offset, SEEK_SET)
    or die "failed to seek in $group_info_path";
  while (my $line = <$group_info_handle>) {
    my $next_offset = tell($group_info_handle);
    if ($next_offset == -1) {
      die "failed to get offset in $group_info_path";
    }

    my @columns = ();
    read_tsv_columns($line, \@columns);
    my $group = $columns[0];
    if (!exists $group_info_offsets_ref->{$group}) {
      $group_info_offsets_ref->{$group} = $curr_offset;
    }
    $curr_offset = $next_offset;
    if ($group > $n) {
      last;
    }
    if ($group < $n) {
      next;
    }

    my $feature = $columns[1];
    if ($feature eq 'read') {
      my $read_input_ID = $columns[2];
      my $ori_ID = $columns[3];
      my $sample = $columns[4];
      my $strand_isoform = $columns[5];
      my $start = $columns[6];
      my $end = $columns[7];
      push @{$group_reads_ref}, $read_input_ID;
      $read_info_ref->{$read_input_ID}{'sample'} = $sample;
      $read_info_ref->{$read_input_ID}{'ori_ID'} = $ori_ID;
      $read_info_ref->{$read_input_ID}{'strand_isoform'} = $strand_isoform;
      $read_info_ref->{$read_input_ID}{'start'} = $start;
      $read_info_ref->{$read_input_ID}{'end'} = $end;
    } elsif ($feature eq 'sj') {
      my $read_input_ID = $columns[2];
      my $ori_ID = $columns[3];
      my $sj_string = $columns[4];
      my $strand_isoform = $columns[5];
      my $is_annotated = $columns[6];
      my ($chr, $start, $end) = split(':', $sj_string);
      my %details = ();
      $details{'chr'} = $chr;
      $details{'start'} = $start;
      $details{'end'} = $end;
      $details{'strand'} = $strand_isoform;
      $details{'is_annotated'} = $is_annotated;
      $SJ_read_improved_ref->{$read_input_ID}{$sj_string} = \%details;
      if (!exists $all_SJ_reads_ref->{$sj_string}{$ori_ID}) {
        $all_SJ_reads_ref->{$sj_string}{$ori_ID} = 0;
      }
    } elsif ($feature eq 'order') {
      my $read_input_ID = $columns[2];
      my $sj_number = $columns[3];
      my $sj_string = $columns[4];
      $SJ_read_order_ref->{$read_input_ID}{$sj_number} = $sj_string;
      if ($sj_string eq 'x') {
        $read_filtered_ref->{$read_input_ID} ++;
      }
    } elsif ($feature eq 'perfect') {
      my $sj_string = $columns[2];
      my $perfect_reads_string = $columns[3];
      my @perfect_reads = split(',', $perfect_reads_string);
      for my $read (@perfect_reads) {
        $all_SJ_reads_ref->{$sj_string}{$read} ++;
      }
    }
  }
  close $group_info_handle;
}

sub add_perfect_reads_to_group_info {
  my ($out_dir, $chr, $keep_tmp, $has_tmp_sj, $sj_to_groups_ref,
      $group_info_path) = @_;

  if (!$has_tmp_sj) {
    return;
  }

  open(my $group_info_handle, '>>', $group_info_path)
    or die "cannot write $group_info_path";

  my $sj_path = "$out_dir/$chr.all_SJ.tmp";
  open(my $sj_handle,  '<', $sj_path) or die "cannot open $sj_path: $!";
  while (<$sj_handle>) {
    my @line = ();
    read_tsv_columns($_, \@line);

    my $perfect_reads_string = $line[-2];
    if (($perfect_reads_string eq 'NA') or ($perfect_reads_string eq '')) {
      next;
    }

    my $sj_string = $line[2];
    if (!exists $sj_to_groups_ref->{$sj_string}) {
      next;
    }
    for my $group (keys %{$sj_to_groups_ref->{$sj_string}}) {
      write_tsv_columns($group_info_handle, [$group, 'perfect', $sj_string,
                                             $perfect_reads_string]);
    }
  }
  close $sj_handle;
  unlink($sj_path) if !$keep_tmp;
}

sub sort_single_exon_isoforms {
  my ($SJ_dist, $single_exon_isoform_end_ref, $group_sort_ref, $group_ends_ref,
      $single_exon_isoforms_by_group_ref) = @_;

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

  my $extra_length = $SJ_dist;
  my $num_groups = scalar(@{$group_sort_ref});
  # Since the isoforms and groups are sorted by coordinate, keep track
  # of which groups end before the remaining isoforms start
  my $min_group_i = 0;
  for my $isoform (@single_exon_isoforms_chr_sort) {
    my $isoform_info = $single_exon_isoform_end_ref->{$isoform};
    my $isoform_start = $isoform_info->{'0'};
    my $isoform_end = $isoform_info->{'1'};
    while ($min_group_i < $num_groups) {
      my $min_group_id = $group_sort_ref->[$min_group_i];
      my $min_group_info = $group_ends_ref->{$min_group_id};
      my $min_group_end = $min_group_info->[1] + $extra_length;
      if ($isoform_start > $min_group_end) {
        # all remaninig isoforms have a start past this
        # group end so stop checking this group
        $min_group_i += 1;
        next;
      }
      last;
    }
    for my $group_i ($min_group_i .. ($num_groups - 1)) {
      my $group_id = $group_sort_ref->[$group_i];
      my $group_info = $group_ends_ref->{$group_id};
      my $group_start = $group_info->[0] - $extra_length;
      my $group_end = $group_info->[1] + $extra_length;
      if ($isoform_start <= $group_end and $isoform_end >= $group_start) {
        push @{$single_exon_isoforms_by_group_ref->{$group_id}}, $isoform;
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

sub sort_splice_junctions {
  my ($anno_SJ_ref, $unanno_SJ_ref, $annotated_SJ_chr_ref,
      $unannotated_SJ_chr_ref) = @_;

  @{$annotated_SJ_chr_ref} = sort {
    my $a_info = $anno_SJ_ref->{$a};
    my $a_start = $a_info->{'start'};
    my $a_end = $a_info->{'end'};
    my $b_info = $anno_SJ_ref->{$b};
    my $b_start = $b_info->{'start'};
    my $b_end = $b_info->{'end'};
    ($a_start <=> $b_start) or ($a_end <=> $b_end);
  } keys %{$anno_SJ_ref};

  @{$unannotated_SJ_chr_ref} = sort {
    my $a_info = $unanno_SJ_ref->{$a};
    my $a_start = $a_info->{'start'};
    my $a_end = $a_info->{'end'};
    my $b_info = $unanno_SJ_ref->{$b};
    my $b_start = $b_info->{'start'};
    my $b_end = $b_info->{'end'};
    ($a_start <=> $b_start) or ($a_end <=> $b_end);
  } keys %{$unanno_SJ_ref};
}

# returns 2 values:
# The 1st indicates the comparison result.
# The 2nd is $num_leading_unpaired_SJs which is the number of SJs from
# the long chain which are not paired to an SJ in the short chain (x or non-x).
sub comp_SJ_chain {
  my ($SJ_chain_short, $SJ_chain_long) = @_;

  my $num_leading_unpaired_SJs = 0;
  my @SJs_short = split('_', $SJ_chain_short);
  my @SJs_long = split('_', $SJ_chain_long);
  my $num_short = scalar(@SJs_short);
  my $num_long = scalar(@SJs_long);
  if (($num_short > $num_long) or ($num_short == 0) or ($num_long == 0)) {
    return (0, $num_leading_unpaired_SJs);
  }

  my $num_non_x = 0;
  my $first_non_x_index;
  for my $short_i (0 .. $#SJs_short) {
    if ($SJs_short[$short_i] eq 'x') {
      next;
    }
    $num_non_x ++;
    if (!defined $first_non_x_index) {
      $first_non_x_index = $short_i;
    }
  }

  my $match_SJ_num = 0;
  my $short_i;
  for my $long_i ($first_non_x_index .. $#SJs_long) {
    if (($match_SJ_num > 0) and (defined $SJs_short[$short_i])
        and ($SJs_short[$short_i] eq $SJs_long[$long_i])) {
      $match_SJ_num ++;
      $short_i ++;
    } elsif ($SJs_short[$first_non_x_index] eq $SJs_long[$long_i]) {
      $short_i = $first_non_x_index;
      my $remaining_long = $#SJs_long - $long_i;
      my $remaining_short = $#SJs_short - $short_i;
      if ($remaining_long < $remaining_short) {
        return (0, $num_leading_unpaired_SJs);
      }
      $match_SJ_num ++;
      $short_i ++;
    } elsif ($match_SJ_num == 0) {
      $num_leading_unpaired_SJs ++;
    } else {
      $short_i ++;
    }
  }

  if ($match_SJ_num == $num_long) {
    return (1, $num_leading_unpaired_SJs);
  }
  if ($match_SJ_num == $num_short) {
    return (2, $num_leading_unpaired_SJs);
  }
  if ($match_SJ_num == $num_non_x) {
    return (3, $num_leading_unpaired_SJs);
  }

  return (0, $num_leading_unpaired_SJs);
}

sub compare_SJs_to_reference {
  my ($SJ_ref, $reference_SJ_ref, $default_handle) = @_;

  my @SJ_sort = sort {
    my $a_info = $SJ_ref->{$a};
    my $a_start = $a_info->{'start'};
    my $a_end = $a_info->{'end'};
    my $b_info = $SJ_ref->{$b};
    my $b_start = $b_info->{'start'};
    my $b_end = $b_info->{'end'};
    ($a_start <=> $b_start) or ($a_end <=> $b_end);
  } keys %{$SJ_ref};

  my @ref_SJ_sort = sort {
    my $a_info = $reference_SJ_ref->{$a};
    my $a_start = $a_info->{'start'};
    my $a_end = $a_info->{'end'};
    my $b_info = $reference_SJ_ref->{$b};
    my $b_start = $b_info->{'start'};
    my $b_end = $b_info->{'end'};
    ($a_start <=> $b_start) or ($a_end <=> $b_end);
  } keys %{$reference_SJ_ref};

  my $SJ_sort_string = join('_', @SJ_sort);
  my $ref_SJ_sort_string = join('_', @ref_SJ_sort);
  my ($result, $num_leading_unpaired_SJs) = comp_SJ_chain($SJ_sort_string,
                                                          $ref_SJ_sort_string);
  if (defined $default_handle) {
    write_tsv_columns($default_handle, [$SJ_sort_string, $ref_SJ_sort_string,
                                        $result]);

  }

  return $result;
}

sub check_fsm_perfect_match_endpoints {
  my ($expected_start, $expected_end, $actual_start, $actual_end, $SJ_dist) = @_;
  return ((abs($actual_start - $expected_start) <= $SJ_dist)
          and (abs($actual_end - $expected_end) <= $SJ_dist));
}

sub classify_reads_and_match_to_annotated_isoforms {
  my ($chr, $n, $SJ_dist, $should_create_default_handle, $anno_SJ_ref,
      $anno_SS_ref, $isoform_SJ_ref, $multi_exon_isoform_end_ref,
      $group_reads_ref, $read_filtered_ref, $SJ_read_improved_ref,
      $read_info_ref, $chr_type_ref, $fsm_ref, $ism_ref, $nic_ref, $nnc_ref,
      $isoform_perfect_match_ref, $single_exon_reads_ref, $summary_data_ref,
      $default_handle) = @_;

  for my $read (@{$group_reads_ref}) {
    if (exists $read_filtered_ref->{$read}) {
      if ($should_create_default_handle) {
        write_tsv_columns($default_handle,
                          [$read, $read_info_ref->{$read}{'ori_ID'}, 'NA',
                           'filtered', 'NA']);
      }
    } elsif (exists $SJ_read_improved_ref->{$read}) {
      my $num_sjs = scalar(keys %{$SJ_read_improved_ref->{$read}});
      my @anno_sjs = grep {
        $SJ_read_improved_ref->{$read}{$_}{'is_annotated'} == 1;
      } keys %{$SJ_read_improved_ref->{$read}};

      if (scalar(@anno_sjs) != $num_sjs) {
        my @same_SSs = grep {
          my $start_string = ${$_}{'start'} . ':' . ${$_}{'strand'};
          my $end_string = ${$_}{'end'} . ':' . ${$_}{'strand'};
          my $start_is_anno = exists $anno_SS_ref->{'0'}{$start_string};
          my $end_is_anno = exists $anno_SS_ref->{'1'}{$end_string};
          $start_is_anno and $end_is_anno;
        } values %{$SJ_read_improved_ref->{$read}};

        if (scalar(@same_SSs) == $num_sjs) {
          if ($should_create_default_handle) {
            write_tsv_columns($default_handle,
                              [$read, $read_info_ref->{$read}{'ori_ID'},
                               scalar(@same_SSs), 'NIC', 'NA', $chr]);
          }
          $nic_ref->{$read} ++;
          $chr_type_ref->{$chr}{'nic'} ++;
        } else {
          if ($should_create_default_handle) {
            write_tsv_columns($default_handle,
                              [$read, $read_info_ref->{$read}{'ori_ID'},
                               scalar(@same_SSs), 'NNC', 'NA', $chr]);
          }
          $nnc_ref->{$read} ++;
          $chr_type_ref->{$chr}{'nnc'} ++;
        }
      } else {
        my @fsm_isoforms = ();
        my @ism_isoforms = ();

        my $info_ref = $read_info_ref->{$read};
        my %possible_isoforms = ();
        for my $SJ2 (keys %{$SJ_read_improved_ref->{$read}}) {
          if (!exists $anno_SJ_ref->{$SJ2}) {
            die "not_recorded:$chr\t$read\t$SJ2";
          }
          if (!defined $anno_SJ_ref->{$SJ2}{'isoforms'}) {
            die "not_array:$chr\t$read\t$SJ2";
          }
          for my $isoform (@{$anno_SJ_ref->{$SJ2}{'isoforms'}}) {
            $possible_isoforms{$isoform} ++;
          }
        }

        for my $isoform (keys %possible_isoforms) {
          my $SJ_isoform_ref = $isoform_SJ_ref->{$isoform};
          my $comp_result = compare_SJs_to_reference(
            $SJ_read_improved_ref->{$read}, $SJ_isoform_ref, $default_handle);
          if ($comp_result == 1) {
            push @fsm_isoforms, $isoform;
          } elsif ($comp_result == 2) {
            push @ism_isoforms, $isoform;
          }
        }
        if (@fsm_isoforms >= 1) {
          if ($should_create_default_handle) {
            my $isoforms_string = join_with_trailing(',', \@fsm_isoforms);
            write_tsv_columns($default_handle,
                              [$read, $info_ref->{$read}{'ori_ID'}, 'all',
                               'FSM', $isoforms_string, $chr]);
          }
          $chr_type_ref->{$chr}{'fsm'} ++;
          for my $isoform (@fsm_isoforms) {
            push @{$fsm_ref->{$read}}, $isoform;

            my $expected_start = $multi_exon_isoform_end_ref->{$isoform}{'0'};
            my $expected_end = $multi_exon_isoform_end_ref->{$isoform}{'1'};
            my $actual_start = $info_ref->{'start'};
            my $actual_end = $info_ref->{'end'};
            if (check_fsm_perfect_match_endpoints($expected_start, $expected_end,
                                                  $actual_start, $actual_end,
                                                  $SJ_dist)) {
              my $sample = $info_ref->{'sample'};
              push @{$isoform_perfect_match_ref->{$isoform}{$sample}}, $read;
              $summary_data_ref->{'num_fsm_perfect_match_reads'} ++;
            }
          }
        } elsif (@ism_isoforms >= 1) {
          if ($should_create_default_handle) {
            my $isoforms_string = join_with_trailing(',', \@ism_isoforms);
            write_tsv_columns($default_handle,
                              [$read, $info_ref->{$read}{'ori_ID'}, 'all',
                               'ISM', $isoforms_string, $chr]);
          }
          for my $isoform (@ism_isoforms) {
            push @{$ism_ref->{$read}}, $isoform;
          }
          $chr_type_ref->{$chr}{'ism'} ++;
        } else {
          if ($should_create_default_handle) {
            write_tsv_columns($default_handle,
                              [$read, 'all', 'NIC', 'NA', $chr]);
          }
          $nic_ref->{$read} ++;
          $chr_type_ref->{$chr}{'nic'} ++;
        }
      }
    } else {
      push @{$single_exon_reads_ref}, $read;
      if ($should_create_default_handle) {
        write_tsv_columns($default_handle,
                          [$read, $read_info_ref->{$read}{'ori_ID'},
                           'all', 'one-exon', 'NA', $chr]);
      }
    }
  }
  $summary_data_ref->{'num_fsm_reads'} += scalar(keys %{$fsm_ref});
  $summary_data_ref->{'num_ism_reads'} += scalar(keys %{$ism_ref});
  $summary_data_ref->{'num_nic_reads'} += scalar(keys %{$nic_ref});
  $summary_data_ref->{'num_nnc_reads'} += scalar(keys %{$nnc_ref});
  $summary_data_ref->{'num_single_exon_reads'} += (
    scalar(@{$single_exon_reads_ref}));
}

sub get_sjs_for_group {
  my ($start_i, $group_start, $group_end, $is_annotated, $sorted_sjs_ref,
      $sj_details_ref, $out_sj_ref) = @_;

  my $num_sjs = scalar(@{$sorted_sjs_ref});
  for my $i ($start_i .. ($num_sjs - 1)) {
    my $sj = $sorted_sjs_ref->[$i];
    my $sj_details = $sj_details_ref->{$sj};
    if ($sj_details->{'start'} > $group_end) {
      last;
    } elsif (($sj_details->{'end'} >= $group_start)
             and ($sj_details->{'start'} <= $group_end)) {
      my %out_details = ();
      $out_details{'start'} = $sj_details->{'start'};
      $out_details{'end'} = $sj_details->{'end'};
      $out_details{'isoforms'} = $sj_details->{'isoforms'};
      $out_details{'is_anno'} = $is_annotated;
      $out_sj_ref->{$sj} = \%out_details;
    } elsif ($sj_details->{'end'} < $group_start) {
      $start_i = $i;
    }
  }

  return $start_i;
}

sub add_other_sjs_for_annotated_isoforms {
  my ($anno_SJ_ref, $isoform_SJ_ref, $annotated_isoforms_included_ref,
      $out_sj_ref) = @_;

  while (my ($sj, $details) = each %{$out_sj_ref}) {
    if ($details->{'is_anno'}) {
      for my $isoform (@{$details->{'isoforms'}}) {
        $annotated_isoforms_included_ref->{$isoform} ++;
      }
    }
  }

  for my $isoform (keys %{$annotated_isoforms_included_ref}) {
    for my $sj (keys %{$isoform_SJ_ref->{$isoform}}) {
      if (exists $out_sj_ref->{$sj}) {
        next;
      }

      my $sj_details = $anno_SJ_ref->{$sj};
      my %new_details = ();
      $new_details{'start'} = $sj_details->{'start'};
      $new_details{'end'} = $sj_details->{'end'};
      $new_details{'isoforms'} = $sj_details->{'isoforms'};
      $new_details{'is_anno'} = 1;
      $out_sj_ref->{$sj} = \%new_details;
    }
  }
}

sub create_sorted_sj_mapping {
  my ($sj_ref, $sorted_sj_ref, $sj_map_ref) = @_;

  @{$sorted_sj_ref} = sort {
    my $a_info = $sj_ref->{$a};
    my $a_start = $a_info->{'start'};
    my $a_end = $a_info->{'end'};
    my $b_info = $sj_ref->{$b};
    my $b_start = $b_info->{'start'};
    my $b_end = $b_info->{'end'};
    ($a_start <=> $b_start) or ($a_end <=> $b_end);
  } keys %{$sj_ref};

  my $num_sjs = scalar(@{$sorted_sj_ref});
  for my $i (0 .. ($num_sjs - 1)) {
    my $sj = $sorted_sj_ref->[$i];
    $sj_map_ref->{$sj} = $i;
  }
}

sub create_sj_id_chains {
  my ($should_create_default_handle, $isoform_SJ_ref, $isoforms_included_ref,
      $sj_id_ref, $out_chain_ref, $default_handle) = @_;

  for my $isoform (keys %{$isoforms_included_ref}) {
    my @sj_ids = map {$sj_id_ref->{$_}} keys %{$isoform_SJ_ref->{$isoform}};
    my @sj_id_sort = sort {$a <=> $b} (@sj_ids);
    my $chain = join('_', @sj_id_sort);
    $out_chain_ref->{$isoform} = $chain;
    if ($should_create_default_handle) {
      print $default_handle "annotated_isoform: $isoform\t$chain\n";
    }
  }
}

sub create_sj_chain_mappings {
  my ($n, $group_reads_ref, $SJ_read_order_ref, $all_SJ_group_ID_ref,
      $read_info_ref, $fsm_ref, $ism_ref, $nic_ref, $nnc_ref, $SJ_chain_read_ref,
      $SJ_chain_cat_ref, $summary_data_ref) = @_;

  for my $read (@{$group_reads_ref}) {
    if (!exists $SJ_read_order_ref->{$read}) {
      next;
    }

    my $num_sj = scalar (keys %{$SJ_read_order_ref->{$read}});
    if (exists $SJ_read_order_ref->{$read}{'-1'}) {
      $SJ_read_order_ref->{$read}{$num_sj} = $SJ_read_order_ref->{$read}{'-1'};
      delete $SJ_read_order_ref->{$read}{'-1'};
    }

    my @SJ_order_sort = sort {$a <=> $b} keys %{$SJ_read_order_ref->{$read}};
    my @SJ_IDs = ();
    for my $SJ_order (@SJ_order_sort) {
      if ($SJ_read_order_ref->{$read}{$SJ_order} eq 'x') {
        push @SJ_IDs, 'x';
        next;
      }
      my $sj_string = $SJ_read_order_ref->{$read}{$SJ_order};
      if (!exists $all_SJ_group_ID_ref->{$sj_string}) {
        die "$read has no SJ ID for: $sj_string";
      }
      push @SJ_IDs, $all_SJ_group_ID_ref->{$sj_string};
    }
    my $SJ_chain = join('_', @SJ_IDs);

    my %read_details = ();
    $read_details{'id'} = $read;
    $read_details{'start'} = $read_info_ref->{$read}{'start'};
    $read_details{'end'} = $read_info_ref->{$read}{'end'};
    if (!exists $SJ_chain_read_ref->{$SJ_chain}) {
      my $cat;
      my %chain_cat_details = ();
      $chain_cat_details{'num_sj'} = $num_sj;
      if (exists $fsm_ref->{$read}) {
        $chain_cat_details{'isoforms'} = $fsm_ref->{$read};
        $SJ_chain_cat_ref->{'fsm'}{$SJ_chain} = \%chain_cat_details;
        $cat = 'fsm';
      } elsif (exists $ism_ref->{$read}) {
        $chain_cat_details{'isoforms'} = $ism_ref->{$read};
        $SJ_chain_cat_ref->{'ism'}{$SJ_chain} = \%chain_cat_details;
        $cat = 'ism';
      } elsif (exists $nic_ref->{$read}) {
        $SJ_chain_cat_ref->{'nic'}{$SJ_chain} = \%chain_cat_details;
        $cat = 'nic';
      } elsif (exists $nnc_ref->{$read}) {
        $SJ_chain_cat_ref->{'nnc'}{$SJ_chain} = \%chain_cat_details;
        $cat = 'nnc';
      } else {
        $SJ_chain_cat_ref->{'filtered'}{$SJ_chain} = \%chain_cat_details;
        $cat = 'filtered';
      }
      my %chain_read_details = ();
      $chain_read_details{'cat'} = $cat;
      $chain_read_details{'start'} = $read_info_ref->{$read}{'start'};
      $chain_read_details{'end'} = $read_info_ref->{$read}{'end'};
      $chain_read_details{'strand_isoform'} = (
        $read_info_ref->{$read}{'strand_isoform'});
      $chain_read_details{'num_sj'} = $num_sj;
      $chain_read_details{'compat_chains'} = {};
      $chain_read_details{'reads'} = [\%read_details];
      $SJ_chain_read_ref->{$SJ_chain} = \%chain_read_details;
    } else {
      push @{$SJ_chain_read_ref->{$SJ_chain}{'reads'}}, \%read_details;
    }
  }
  $summary_data_ref->{'num_fsm_chains'} += scalar(
    keys %{$SJ_chain_cat_ref->{'fsm'}});
  $summary_data_ref->{'num_ism_chains'} += scalar(
    keys %{$SJ_chain_cat_ref->{'ism'}});
  $summary_data_ref->{'num_nic_chains'} += scalar(
    keys %{$SJ_chain_cat_ref->{'nic'}});
  $summary_data_ref->{'num_nnc_chains'} += scalar(
    keys %{$SJ_chain_cat_ref->{'nnc'}});
  $summary_data_ref->{'num_chains_with_failed_SJ'} += scalar(
    keys %{$SJ_chain_cat_ref->{'filtered'}});
}

sub create_endpoint_mapping_by_chain {
  my ($chr, $raw, $should_create_default_handle, $all_SJ_group_sort_ref,
      $SJ_chain_cat_ref, $SJ_chain_read_ref, $SJ_chain_read_exon_end_ref,
      $SJ_included_chain_ref, $default_handle) = @_;

  while (my ($SJ_chain, $SJ_chain_info) = each %{$SJ_chain_read_ref}) {
    my @SJ_IDs = split('_', $SJ_chain);
    # Use the median start and end coordinates as the chain start and end.
    # For an even number of reads there are two middle reads.
    # Choose the middle read separately for start and end,
    # and pick the middle read that makes the exon longer.
    my @start_reads_sort = sort {
      $a->{'start'} <=> $b->{'start'};
    } @{$SJ_chain_info->{'reads'}};
    my @end_reads_sort = sort {
      $a->{'end'} <=> $b->{'end'};
    } @{$SJ_chain_info->{'reads'}};

    my $middle_index = int(scalar(@start_reads_sort) / 2);
    my $middle_start_index = $middle_index;
    my $middle_end_index = $middle_index;
    my $has_even_read_count = (scalar(@start_reads_sort) % 2) == 0;
    if ($has_even_read_count) {
      $middle_start_index--;
    }

    my $chain_start = $start_reads_sort[$middle_start_index]->{'start'};
    my $chain_end = $end_reads_sort[$middle_end_index]->{'end'};
    $SJ_chain_info->{'start'} = $chain_start;
    $SJ_chain_info->{'end'} = $chain_end;
    if (defined $raw) {
      next;
    }

    my $is_filtered_chain = exists $SJ_chain_cat_ref->{'filtered'}{$SJ_chain};
    for my $SJ_ID_i (0 .. $#SJ_IDs) {
      my $SJ_ID = $SJ_IDs[$SJ_ID_i];
      if ($SJ_ID ne 'x') {
        push @{$SJ_included_chain_ref->{$SJ_ID}}, $SJ_chain;
      }

      if ($is_filtered_chain) {
        next;
      }

      my $start;
      my $end;
      if ($SJ_ID_i == 0) {
        $start = $chain_start;
      } else {
        my $prev_sj_id = $SJ_IDs[$SJ_ID_i - 1];
        my $prev_SJ = $all_SJ_group_sort_ref->[$prev_sj_id];
        my $sj_details = split_sj_string($prev_SJ);
        $start = $sj_details->{'end'};
      }

      if ($SJ_ID_i == $#SJ_IDs) {
        $end = $chain_end;
      } else {
        my $next_sj_id = $SJ_IDs[$SJ_ID_i + 1];
        my $next_SJ = $all_SJ_group_sort_ref->[$next_sj_id];
        my $sj_details = split_sj_string($next_SJ);
        $end = $sj_details->{'start'};
      }

      $SJ_chain_read_exon_end_ref->{$SJ_chain}{'start'}{$SJ_ID} = $start;
      $SJ_chain_read_exon_end_ref->{$SJ_chain}{'end'}{$SJ_ID} = $end;
      if ($should_create_default_handle) {
        print $default_handle "$chr:exon_end: $SJ_chain\t$SJ_ID\t$start\t$end\n";
      }
    }
  }
}

sub set_candidate_compatible_chains {
  my ($SJ_included_chain_ref, $SJ_chain_cat_ref, $SJ_chain_read_ref) = @_;

  # For each chain check all other chains that share at least 1 SJ.
  # A chain won't be marked as compatible to a chain with an 'x' SJ.
  # A chain won't be marked as compatible to a shorter chain.
  # Only a chain that has an 'x' SJ can be marked as compatible to
  # a chain with the same length.
  while (my ($SJ_ID, $SJ_chain_ref) = each %{$SJ_included_chain_ref}) {
    my @chain_length_sort = sort {
      $SJ_chain_read_ref->{$a}{'num_sj'} <=> $SJ_chain_read_ref->{$b}{'num_sj'};
    } @{$SJ_chain_ref};

    my $num_chains = scalar(@chain_length_sort);
    for my $i (0 .. ($num_chains - 2)) {
      my $short_chain = $chain_length_sort[$i];
      my $short_num_sj = $SJ_chain_read_ref->{$short_chain}{'num_sj'};
      my $short_is_filtered = exists(
        $SJ_chain_cat_ref->{'filtered'}{$short_chain});
      for my $j (($i + 1) .. ($num_chains - 1)) {
        my $long_chain = $chain_length_sort[$j];
        my $long_num_sj = $SJ_chain_read_ref->{$long_chain}{'num_sj'};
        my $long_is_filtered = exists(
          $SJ_chain_cat_ref->{'filtered'}{$long_chain});
        if ($short_num_sj == $long_num_sj) {
          if ($long_is_filtered and !$short_is_filtered) {
            $SJ_chain_read_ref->{$long_chain}{'compat_chains'}{$short_chain} ++;
          } elsif ($short_is_filtered and !$long_is_filtered) {
            $SJ_chain_read_ref->{$short_chain}{'compat_chains'}{$long_chain} ++;
          }
        } elsif (!$long_is_filtered) {
          $SJ_chain_read_ref->{$short_chain}{'compat_chains'}{$long_chain} ++;
        }
      }
    }
  }
}

sub check_perfect_reads_for_SJs_and_match_to_isoforms {
  my ($chain, $SJ_reads_ref, $reads_for_chain_ref, $read_info_ref,
      $all_SJ_group_sort_ref, $anno_SJ_ref, $SJ_isoform_SJ_ref, $read_num_cutoff,
      $read_ratio_cutoff, $default_handle) = @_;

  my @SJ_IDs = split('_', $chain);
  if (defined $default_handle) {
    print $default_handle "validation_for_nc_chain: $chain";
  }

  my @SJs = ();
  my %possible_isoforms = ();
  my %SJ_info_nc = ();
  for my $SJ_ID (@SJ_IDs) {
    my $SJ = $all_SJ_group_sort_ref->[$SJ_ID];
    push @SJs, $SJ;

    my $sj_info = split_sj_string($SJ);
    my %sj_details = ();
    $sj_details{'start'} = $sj_info->{'start'};
    $sj_details{'end'} = $sj_info->{'end'};
    $SJ_info_nc{$SJ} = \%sj_details;
    if (exists $anno_SJ_ref->{$SJ}) {
      for my $isoform (@{$anno_SJ_ref->{$SJ}{'isoforms'}}) {
        $possible_isoforms{$isoform} ++;
      }
    }

    if (defined $default_handle) {
      print $default_handle "\t$SJ";
    }
    if (!exists $SJ_reads_ref->{$SJ}) {
      die "0:$chain\t$SJ_ID\t$SJ";
    }
  }
  if (defined $default_handle) {
    print $default_handle "\n";
  }

  my @fsm = ();
  my @ism = ();
  for my $possible_fsm (keys %possible_isoforms) {
    my $comp_result = compare_SJs_to_reference(
      \%SJ_info_nc, $SJ_isoform_SJ_ref->{$possible_fsm}, $default_handle);
    if ($comp_result == 1) {
      push @fsm, $possible_fsm;
    } elsif ($comp_result == 2) {
      push @ism, $possible_fsm;
    }
  }
  %SJ_info_nc = ();

  if (defined $default_handle) {
    if (@fsm >= 1) {
      my $fsms = join_with_trailing(',', \@fsm);
      print $default_handle "unannotated_isoforms $chain: $fsms\n";
    } elsif (@ism >= 1) {
      my $isms = join_with_trailing(',', \@ism);
      print $default_handle "unannotated_isoforms_fragment $chain: $isms\n";
    }
  }

  my $num_reads = scalar(@{$reads_for_chain_ref});
  if ($num_reads == 0) {
    die "no reads loaded for $chain";
  }

  my @SJs_valid = ();
  my @perfect_read_SJs = ();
  for my $SJ (@SJs) {
    my @perfect_reads = grep {
      my $id = $read_info_ref->{${$_}{'id'}}{'ori_ID'};
      $SJ_reads_ref->{$SJ}{$id} >= 1; # is perfect read for this SJ
    } @{$reads_for_chain_ref};

    my $num_perfect_reads = scalar(@perfect_reads);
    push @perfect_read_SJs, $num_perfect_reads;
    if (($num_perfect_reads >= $read_num_cutoff)
        and (($num_perfect_reads / $num_reads) >= $read_ratio_cutoff)) {
      push @SJs_valid, $SJ;
    }
    if (defined $default_handle) {
      printf $default_handle "\t$SJ(%g/%g)", $num_perfect_reads, $num_reads;
    }
  }

  if (defined $default_handle) {
    print $default_handle "\n";
  }
  my @perfect_read_SJs_sort = sort {$a <=> $b} @perfect_read_SJs;
  my $lowest_perfect_count = $perfect_read_SJs_sort[0];
  my %result = ();
  $result{'fsm'} = \@fsm;
  $result{'ism'} = \@ism;
  $result{'lowest_perfect_count'} = $lowest_perfect_count;
  $result{'num_reads'} = $num_reads;
  if (@SJs_valid == @SJs) {
    $result{'is_valid'} = 1;
    if (defined $default_handle) {
      print $default_handle "yes_validated: $chain\n";
    }
    return \%result;
  } else {
    $result{'is_valid'} = 0;
    if (defined $default_handle) {
      print $default_handle "not_validated: $chain\n";
    }
    return \%result;
  }
}

sub find_valid_novel_chains_raw {
  my ($SJ_chain_cat_ref, $all_SJ_reads_ref, $SJ_chain_read_ref, $read_info_ref,
      $all_SJ_group_sort_ref, $anno_SJ_complete_ref, $isoform_SJ_complete_ref,
      $read_num_cutoff, $read_ratio_cutoff, $unannotated_isoforms_ref,
      $longest_nc_valid_ref, $default_handle) = @_;

  for my $cat ('ism', 'nic', 'nnc') {
    if (exists $SJ_chain_cat_ref->{$cat}) {
      for my $chain (keys %{$SJ_chain_cat_ref->{$cat}}) {
        my $result = check_perfect_reads_for_SJs_and_match_to_isoforms(
          $chain, $all_SJ_reads_ref, $SJ_chain_read_ref->{$chain}{'reads'},
          $read_info_ref, $all_SJ_group_sort_ref, $anno_SJ_complete_ref,
          $isoform_SJ_complete_ref, $read_num_cutoff, $read_ratio_cutoff,
          $default_handle);
        $unannotated_isoforms_ref->{$chain} = $result;
        $longest_nc_valid_ref->{$chain} = -1 if $result->{'is_valid'};
      }
    }
  }
}

sub find_valid_novel_chains {
  my ($SJ_chain_ref, $all_SJ_reads_ref, $SJ_chain_read_ref, $read_info_ref,
      $all_SJ_group_sort_ref, $anno_SJ_complete_ref, $isoform_SJ_complete_ref,
      $read_num_cutoff, $read_ratio_cutoff, $unannotated_isoforms_ref,
      $valid_chains_ref, $default_handle) = @_;

  # Sort on num_sj, with compare_SJ_chains_as_numeric_tuples to break ties.
  # The tie break ensures consistent novel isoform IDs.
  my @chain_length_sort = sort {
    my $len_a = $SJ_chain_ref->{$a}{'num_sj'};
    my $len_b = $SJ_chain_ref->{$b}{'num_sj'};
    ($len_a <=> $len_b) or compare_SJ_chains_as_numeric_tuples($a, $b);
  } keys %{$SJ_chain_ref};

  for my $chain(@chain_length_sort) {
    my $result = check_perfect_reads_for_SJs_and_match_to_isoforms(
      $chain, $all_SJ_reads_ref, $SJ_chain_read_ref->{$chain}{'reads'},
      $read_info_ref, $all_SJ_group_sort_ref, $anno_SJ_complete_ref,
      $isoform_SJ_complete_ref, $read_num_cutoff, $read_ratio_cutoff,
      $default_handle);
    if ($result->{'is_valid'}) {
      push @{$valid_chains_ref}, $chain;
      $unannotated_isoforms_ref->{$chain} = $result;
    }
  }
}

sub check_endpoint_compatibility {
  my ($expected_start, $is_expected_start_terminal, $expected_end,
      $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
      $internal_boundary_limit, $allow_longer_terminal_exons,
      $summary_data_ref) = @_;

  my $is_ok = 1;
  my $amount_past_start = $expected_start - $actual_start;
  if ($is_expected_start_terminal) {
    if ((!$allow_longer_terminal_exons) and ($amount_past_start > $SJ_dist)) {
      $summary_data_ref->{'num_terminal_exon_boundary_check_fails'} ++;
      $is_ok = 0;
    }
  } elsif ($amount_past_start > $internal_boundary_limit) {
    $summary_data_ref->{'num_internal_exon_boundary_check_fails'} ++;
    $is_ok = 0;
  }

  my $amount_past_end = $actual_end - $expected_end;
  if ($is_expected_end_terminal) {
    if ((!$allow_longer_terminal_exons) and ($amount_past_end > $SJ_dist)) {
      $summary_data_ref->{'num_terminal_exon_boundary_check_fails'} ++;
      $is_ok = 0;
    }
  } elsif ($amount_past_end > $internal_boundary_limit) {
    $summary_data_ref->{'num_internal_exon_boundary_check_fails'} ++;
    $is_ok = 0
  }

  return $is_ok;
}

sub check_ism_chain_endpoints_and_valid_sjs {
  my ($SJ_dist, $internal_boundary_limit, $allow_longer_terminal_exons,
      $should_create_default_handle, $fsm_chain_ref, $anno_SJ_ref,
      $multi_exon_isoform_end_ref, $all_SJ_reads_ref, $SJ_chain_read_ref,
      $read_info_ref, $all_SJ_group_sort_ref, $anno_SJ_complete_ref,
      $isoform_SJ_complete_ref, $read_num_cutoff, $read_ratio_cutoff,
      $SJ_chain_ref, $unannotated_isoforms_ref, $possible_nc_ism_ref,
      $chain_w_compatible_annotated_fsm_ref, $summary_data_ref,
      $default_handle) = @_;

  # Sort on chain length, but break ties with compare_SJ_chains_as_numeric_tuples.
  # The tie break ensures consistent novel isoform IDs.
  my @chain_length_sort = sort {
    my $len_a = $SJ_chain_ref->{$a}{'num_sj'};
    my $len_b = $SJ_chain_ref->{$b}{'num_sj'};
    ($len_a <=> $len_b) or compare_SJ_chains_as_numeric_tuples($a, $b);
  } keys %{$SJ_chain_ref};

  for my $chain (@chain_length_sort) {
    my @isoforms_fsm_valid = ();
    my @SJs_i = split('_', $chain);
    if ($should_create_default_handle) {
      my $isoforms = join_with_trailing(',',
                                        $SJ_chain_ref->{$chain}{'isoforms'});
      print $default_handle "ism\t$chain\t$isoforms";
    }
    my $compatible_annotated_fsm_num = 0;
    for my $isoform_fsm (@{$SJ_chain_ref->{$chain}{'isoforms'}}) {
      if (!exists $multi_exon_isoform_end_ref->{$isoform_fsm}) {
        next;
      }
      my $end_ref = $multi_exon_isoform_end_ref->{$isoform_fsm};
      my $first_sj = $all_SJ_group_sort_ref->[$SJs_i[0]];
      my $first_sj_start = $anno_SJ_ref->{$first_sj}{'start'};
      my $last_sj = $all_SJ_group_sort_ref->[$SJs_i[-1]];
      my $last_sj_end = $anno_SJ_ref->{$last_sj}{'end'};
      if ((!exists $end_ref->{$first_sj_start})
          or (!exists $end_ref->{$last_sj_end})) {
        next;
      }
      my $expected_first_start = $end_ref->{'0'};
      my $expected_last_end = $end_ref->{'1'};
      my $expected_start = $end_ref->{$first_sj_start};
      my $is_expected_start_terminal = $expected_start == $expected_first_start;
      my $expected_end = $end_ref->{$last_sj_end};
      my $is_expected_end_terminal = $expected_end == $expected_last_end;
      my $actual_start = $SJ_chain_read_ref->{$chain}{'start'};
      my $actual_end = $SJ_chain_read_ref->{$chain}{'end'};
      if (check_endpoint_compatibility(
            $expected_start, $is_expected_start_terminal, $expected_end,
            $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
            $internal_boundary_limit, $allow_longer_terminal_exons,
            $summary_data_ref)) {
        $compatible_annotated_fsm_num ++;
        if (exists $SJ_chain_read_ref->{$fsm_chain_ref->{$isoform_fsm}}) {
          # There is a read for the annotated chain
          push @isoforms_fsm_valid, $isoform_fsm;
        }
      }
    }

    if ($should_create_default_handle) {
      my $isoforms = join_with_trailing(',', \@isoforms_fsm_valid);
      print $default_handle "\t$isoforms\n";
    }
    $SJ_chain_ref->{$chain}{'isoforms'} = \@isoforms_fsm_valid;
    if (scalar(@isoforms_fsm_valid) == 0) {
      # This chain could be a valid novel chain
      if ($compatible_annotated_fsm_num > 0) {
        $chain_w_compatible_annotated_fsm_ref->{$chain} = 1;
      }
      my $result = check_perfect_reads_for_SJs_and_match_to_isoforms(
        $chain, $all_SJ_reads_ref, $SJ_chain_read_ref->{$chain}{'reads'},
        $read_info_ref, $all_SJ_group_sort_ref, $anno_SJ_complete_ref,
        $isoform_SJ_complete_ref, $read_num_cutoff, $read_ratio_cutoff,
        $default_handle);
      if ($result->{'is_valid'}) {
        push @{$possible_nc_ism_ref}, $chain;
        $unannotated_isoforms_ref->{$chain} = $result;
      }
    }
  }
}

sub find_sub_chains {
  my ($SJ_dist, $internal_boundary_limit, $allow_longer_terminal_exons,
      $should_create_default_handle, $chains_ref, $SJ_chain_read_ref,
      $SJ_chain_read_exon_end_ref, $substring_nc_ref, $summary_data_ref,
      $default_handle) = @_;

  # chains are sorted by length (number of SJs)
  my $num_chains = scalar(@{$chains_ref});
  for my $i (0 .. ($num_chains - 2)) {
    my $chain_i = $chains_ref->[$i];
    my @SJs_i = split('_', $chain_i);
    my $first_sj_i = $SJs_i[0];
    my $last_sj_i = $SJs_i[-1];
    my $chain_i_details = $SJ_chain_read_ref->{$chain_i};
    my $chain_i_start = $chain_i_details->{'start'};
    my $chain_i_end = $chain_i_details->{'end'};
    for my $j (($i + 1) .. ($num_chains - 1)) {
      my $chain_j = $chains_ref->[$j];
      if (!exists $chain_i_details->{'compat_chains'}{$chain_j}) {
        next;
      }

      my @SJs_j = split('_', $chain_j);
      my $first_sj_j = $SJs_j[0];
      my $last_sj_j = $SJs_j[-1];
      my $chain_j_exon_details = $SJ_chain_read_exon_end_ref->{$chain_j};
      my $chain_j_exon_start = $chain_j_exon_details->{'start'}{$first_sj_i};
      my $chain_j_exon_end = $chain_j_exon_details->{'end'}{$last_sj_i};

      my ($comp_result, $num_leading_unpaired_SJs) = comp_SJ_chain(
        $chain_i, $chain_j);
      if ($should_create_default_handle) {
        write_tsv_columns($default_handle,
                          [$chain_i, $chain_j, $comp_result, $chain_j_exon_start,
                           $chain_j_exon_end, $chain_i_start, $chain_i_end]);
      }
      if ($comp_result == 2) {
        my $expected_start = $chain_j_exon_start;
        my $is_expected_start_terminal = $first_sj_i == $first_sj_j;
        my $expected_end = $chain_j_exon_end;
        my $is_expected_end_terminal = $last_sj_i == $last_sj_j;
        my $actual_start = $chain_i_start;
        my $actual_end = $chain_i_end;

        my $ends_ok = check_endpoint_compatibility(
          $expected_start, $is_expected_start_terminal, $expected_end,
          $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
          $internal_boundary_limit, $allow_longer_terminal_exons,
          $summary_data_ref);
        if ($ends_ok) {
          push @{$substring_nc_ref->{$chain_i}}, $chain_j;
        }
      }
    }
  }
}

sub add_longest_nnc_valid_chains {
  my ($nnc_valid_all_ref, $substring_nc_ref, $nc_valid_longest_ref,
      $longest_nc_valid_ref) = @_;

  my $num_longest = scalar(@{$nc_valid_longest_ref});
  for my $chain (@{$nnc_valid_all_ref}) {
    if (!exists $substring_nc_ref->{$chain}) {
      push @{$nc_valid_longest_ref}, $chain;
      $longest_nc_valid_ref->{$chain} = $num_longest;
      $num_longest ++;
    }
  }
}

sub add_longest_nic_valid_chains {
  my ($SJ_dist, $internal_boundary_limit, $allow_longer_terminal_exons,
      $nic_valid_all_ref, $SJ_chain_read_ref, $SJ_chain_read_exon_end_ref,
      $substring_nc_ref, $nc_valid_longest_ref, $longest_nc_valid_ref,
      $summary_data_ref) = @_;

  my %longest_nic_valid = ();
  my $num_longest = scalar(@{$nc_valid_longest_ref});
  for my $chain (@{$nic_valid_all_ref}) {
    my @SJs_i = split('_', $chain);
    my $first_sj_i = $SJs_i[0];
    my $last_sj_i = $SJs_i[-1];
    my $chain_i_details = $SJ_chain_read_ref->{$chain};
    my $chain_i_start = $chain_i_details->{'start'};
    my $chain_i_end = $chain_i_details->{'end'};
    for my $nc_chain (keys %{$longest_nc_valid_ref}) {
      if (!exists $chain_i_details->{'compat_chains'}{$nc_chain}) {
        next;
      }

      my ($comp_result, $num_leading_unpaired_SJs) = comp_SJ_chain(
        $chain, $nc_chain);
      if ($comp_result == 2) {
        my @SJs_j = split('_', $nc_chain);
        my $first_sj_j = $SJs_j[0];
        my $last_sj_j = $SJs_j[-1];
        my $chain_j_exon_details = $SJ_chain_read_exon_end_ref->{$nc_chain};
        my $chain_j_exon_start = $chain_j_exon_details->{'start'}{$first_sj_i};
        my $chain_j_exon_end = $chain_j_exon_details->{'end'}{$last_sj_i};

        my $expected_start = $chain_j_exon_start;
        my $is_expected_start_terminal = $first_sj_i == $first_sj_j;
        my $expected_end = $chain_j_exon_end;
        my $is_expected_end_terminal = $last_sj_i == $last_sj_j;
        my $actual_start = $chain_i_start;
        my $actual_end = $chain_i_end;

        my $ends_ok = check_endpoint_compatibility(
          $expected_start, $is_expected_start_terminal, $expected_end,
          $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
          $internal_boundary_limit, $allow_longer_terminal_exons,
          $summary_data_ref);
       if ($ends_ok) {
          push @{$substring_nc_ref->{$chain}}, $nc_chain;
        }
      }
    }

    if (!exists $substring_nc_ref->{$chain}) {
      push @{$nc_valid_longest_ref}, $chain;
      $longest_nic_valid{$chain} = $num_longest;
      $num_longest ++;
    }
  }

  for my $chain (keys %longest_nic_valid) {
    $longest_nc_valid_ref->{$chain} = $longest_nic_valid{$chain};
  }
}

sub add_longest_ism_valid_chains {
  my ($SJ_dist, $internal_boundary_limit, $allow_longer_terminal_exons,
      $possible_nc_ism_ref, $SJ_chain_read_ref, $SJ_chain_read_exon_end_ref,
      $chain_w_compatible_annotated_fsm_ref, $substring_nc_ref,
      $nc_valid_longest_ref, $longest_nc_valid_ref, $SJ_chain_ref,
      $summary_data_ref) = @_;

  my $num_longest = scalar(@{$nc_valid_longest_ref});
  # Process possible_nc_ism chains from longest to shortest in this loop.
  # That way the longer chains can be added as nc_valid_longest before the
  # shorter chains are checked against all nc_valid_longest chains.
  for my $chain (reverse @{$possible_nc_ism_ref}) {
    my @SJs_i = split('_', $chain);
    my $first_sj_i = $SJs_i[0];
    my $last_sj_i = $SJs_i[-1];
    my $chain_i_details = $SJ_chain_read_ref->{$chain};
    my $chain_i_start = $chain_i_details->{'start'};
    my $chain_i_end = $chain_i_details->{'end'};
    for my $nc_chain (@{$nc_valid_longest_ref}) {
      if (!exists $chain_i_details->{'compat_chains'}{$nc_chain}) {
        next;
      }
      # nc_chain could be a possible_nc_ism chain.
      # In that case it may already be added as a substring.
      if (grep {/^$nc_chain$/} @{$substring_nc_ref->{$chain}}) {
        next;
      }

      my ($comp_result, $num_leading_unpaired_SJs) = comp_SJ_chain(
        $chain, $nc_chain);
      if ($comp_result == 2) {
        my @SJs_j = split('_', $nc_chain);
        my $first_sj_j = $SJs_j[0];
        my $last_sj_j = $SJs_j[-1];
        my $chain_j_exon_details = $SJ_chain_read_exon_end_ref->{$nc_chain};
        my $chain_j_exon_start = $chain_j_exon_details->{'start'}{$first_sj_i};
        my $chain_j_exon_end = $chain_j_exon_details->{'end'}{$last_sj_i};

        my $expected_start = $chain_j_exon_start;
        my $is_expected_start_terminal = $first_sj_i == $first_sj_j;
        my $expected_end = $chain_j_exon_end;
        my $is_expected_end_terminal = $last_sj_i == $last_sj_j;
        my $actual_start = $chain_i_start;
        my $actual_end = $chain_i_end;

        my $ends_ok = check_endpoint_compatibility(
          $expected_start, $is_expected_start_terminal, $expected_end,
          $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
          $internal_boundary_limit, $allow_longer_terminal_exons,
          $summary_data_ref);
        if ($ends_ok) {
          push @{$substring_nc_ref->{$chain}}, $nc_chain;
        }
      }
    }

    if (exists $substring_nc_ref->{$chain}) {
      my @nc_isoforms = ();
      for my $nc_chain (@{$substring_nc_ref->{$chain}}) {
        if (exists $longest_nc_valid_ref->{$nc_chain}) {
          push @nc_isoforms, $longest_nc_valid_ref->{$nc_chain};
        }
      }
      $SJ_chain_ref->{$chain}{'isoforms'} = \@nc_isoforms;
    } elsif (exists $chain_w_compatible_annotated_fsm_ref->{$chain}) {
      push @{$nc_valid_longest_ref}, $chain;
      $longest_nc_valid_ref->{$chain} = $num_longest;
      $num_longest ++;
    }
  }
}

sub initialize_isoform_assignment_counts {
  my ($SJ_chain_read_ref, $read_info_ref, $samples_sort_ref,
      $read_count_isoform_ref) = @_;

  while (my ($SJ_chain, $chain_info_ref) = each %{$SJ_chain_read_ref}) {
    my $reads_ref = $chain_info_ref->{'reads'};
    my %read_count_sample = ();
    for my $read_details (@{$reads_ref}) {
      my $sample = $read_info_ref->{$read_details->{'id'}}{'sample'};
      $read_count_sample{$sample} ++;
    }

    for my $sample (@{$samples_sort_ref}) {
      my $count = 0;
      if (exists $read_count_sample{$sample}) {
        $count = $read_count_sample{$sample};
      }

      my %details = ();
      $details{'initial'} = $count;
      $details{'direct'} = $count;
      $details{'EM'} = 0;
      $details{'current'} = 0;
      $details{'previous'} = 0;
      $read_count_isoform_ref->{$SJ_chain}{$sample} = \%details;
    }
  }
}

sub output_chain_details_to_default {
  my ($chr, $should_create_default_handle, $SJ_chain_read_ref,
      $all_SJ_group_sort_ref, $SJ_chain_cat_ref, $longest_nc_valid_ref,
      $substring_nc_ref, $unannotated_isoforms_ref, $samples_sort_ref,
      $read_count_isoform_ref, $default_handle) = @_;
  if (!$should_create_default_handle) {
    return;
  }

  while (my ($SJ_chain, $chain_info_ref) = each %{$SJ_chain_read_ref}) {
    my @SJ_IDs = split('_', $SJ_chain);
    my $cat = $chain_info_ref->{'cat'};
    my $chain_start = $chain_info_ref->{'start'};
    my $chain_end = $chain_info_ref->{'end'};
    my $read_info = $chain_info_ref->{'reads'};
    my $isoform_chain_string = join("\t", ('isoform_chain', $chr, $SJ_chain,
                                           $cat, scalar(@{$read_info}),
                                           $chain_start, $chain_end));
    $isoform_chain_string = $isoform_chain_string . "\t";
    print $default_handle $isoform_chain_string;

    for my $SJ_ID (@SJ_IDs) {
      if ($SJ_ID eq 'x') {
        print $default_handle 'x,';
      } else {
        my $sj = $all_SJ_group_sort_ref->[$SJ_ID];
        print $default_handle "$sj,";
      }
    }
    print $default_handle "\t";

    my $by_chain = $SJ_chain_cat_ref->{$cat};
    if ((exists $by_chain->{$SJ_chain})
        and (defined $by_chain->{$SJ_chain}{'isoforms'})
        and (scalar(@{$by_chain->{$SJ_chain}{'isoforms'}}) > 0)) {
      my $isoforms = join_with_trailing(',', $by_chain->{$SJ_chain}{'isoforms'});
      print $default_handle $isoforms;
    } else {
      print $default_handle 'NA';
    }
    print $default_handle "\t";

    if ($cat eq 'fsm') {
      print $default_handle 'NA';
    } elsif (exists ($longest_nc_valid_ref->{$SJ_chain})) {
      print $default_handle "valid($longest_nc_valid_ref->{$SJ_chain})";
    } else {
      print $default_handle 'non_valid';
    }
    print $default_handle "\t";

    if (exists ($substring_nc_ref->{$SJ_chain})) {
      my $string = join_with_trailing(',', $substring_nc_ref->{$SJ_chain});
      print $default_handle $string;
    } else {
      print $default_handle 'NA';
    }
    print $default_handle "\t";

    if (!exists $unannotated_isoforms_ref->{$SJ_chain}) {
      print $default_handle "NA\tNA\tNA\tNA";
    } else {
      my $details = $unannotated_isoforms_ref->{$SJ_chain};
      if (scalar(@{$details->{'fsm'}}) > 0) {
        my $isoforms = join_with_trailing(',', $details->{'fsm'});
        print $default_handle $isoforms;
      } else {
        print $default_handle 'NA';
      }
      print $default_handle "\t";

      if (scalar(@{$details->{'ism'}}) > 0) {
        my $isoforms = join_with_trailing(',', $details->{'ism'});
        print $default_handle $isoforms;
      } else {
        print $default_handle 'NA';
      }
      print $default_handle "\t";

      print $default_handle "$details->{'lowest_perfect_count'}\t";
      print $default_handle "$details->{'num_reads'}";
    }

    for my $sample (@{$samples_sort_ref}) {
      print $default_handle "\t";
      my $details = $read_count_isoform_ref->{$SJ_chain}{$sample};
      print $default_handle "$details->{'initial'},";
      print $default_handle "$details->{'direct'},";
      print $default_handle "$details->{'EM'},";
      print $default_handle "$details->{'current'},";
      print $default_handle "$details->{'previous'},";
    }
    print $default_handle "\n";
  }
}

sub sort_sj_chains_by_number_of_sjs {
  my ($SJ_chain_cat_ref, $cat_SJ_chains_sort_ref) = @_;

  for my $cat ('ism', 'fsm', 'nic', 'nnc', 'filtered') {
    @{$cat_SJ_chains_sort_ref->{$cat}} = sort {
      my $by_chain = $SJ_chain_cat_ref->{$cat};
      $by_chain->{$b}{'num_sj'} <=> $by_chain->{$a}{'num_sj'};
    } keys %{$SJ_chain_cat_ref->{$cat}};
  }
}

sub assign_reads_to_compatible_nic_and_nnc_chains {
  my ($chr, $n, $should_create_tsv_compt_handle, $should_create_default_handle,
      $cat_SJ_chains_sort_ref, $substring_nc_ref, $longest_nc_valid_ref,
      $nc_valid_longest_ref, $samples_sort_ref, $SJ_chain_read_ref,
      $read_info_ref, $SJ_chain_cat_ref, $read_count_isoform_ref,
      $summary_data_ref, $tsv_compt_handle, $default_handle) = @_;

  for my $cat ('nic', 'nnc') {
    my $compat_class_string = uc($cat);
    for my $chain (@{$cat_SJ_chains_sort_ref->{$cat}}) {
      my $reads_ref = $SJ_chain_read_ref->{$chain}{'reads'};
      my $num_reads = scalar(@{$reads_ref});
      if (exists $substring_nc_ref->{$chain}) {
        my @nc_isoforms = ();
        for my $nc_chain (@{$substring_nc_ref->{$chain}}) {
          if (exists $longest_nc_valid_ref->{$nc_chain}) {
            push @nc_isoforms, $longest_nc_valid_ref->{$nc_chain};
          }
        }

        if (scalar(@nc_isoforms) == 0) {
          if ($should_create_tsv_compt_handle) {
            for my $read_details (@{$reads_ref}) {
              my $read_ID = $read_details->{'id'};
              my $ori_id = $read_info_ref->{$read_ID}{'ori_ID'};
              my $sample = $read_info_ref->{$read_ID}{'sample'};
              write_tsv_columns($tsv_compt_handle, [$ori_id, $sample,
                                                    $compat_class_string, 'NA']);
            }
          }
          $summary_data_ref->{'num_reads_not_assigned'} += $num_reads;
          if ($should_create_default_handle) {
            write_tsv_columns($default_handle, ['read_count_summary', $chr, $cat,
                                                "$n:$chain", 'NA1', $num_reads]);
          }
          next;
        }

        $SJ_chain_cat_ref->{$cat}{$chain}{'isoforms'} = \@nc_isoforms;
        my $type;
        if (scalar(@nc_isoforms) == 1) {
          $type = "direct";
          my $isoform = $nc_valid_longest_ref->[$nc_isoforms[0]];
          my $read_count_total = 0;
          for my $sample (@{$samples_sort_ref}) {
            my $initial = $read_count_isoform_ref->{$chain}{$sample}{'initial'};
            $read_count_isoform_ref->{$isoform}{$sample}{'direct'} += $initial;
            $read_count_total += $initial;
          }
          if ($should_create_default_handle
              and ($read_count_total != $num_reads)) {
            write_tsv_columns($default_handle,
                              ['!!read_count_summary', $chr, $cat, "$n:$chain",
                               $type, "$num_reads!=$read_count_total"]);
          }
        } else {
          $type = "EM";
        }
        if ($should_create_tsv_compt_handle) {
          my @isoform_strings = ();
          for my $isoform_id (@nc_isoforms) {
            push @isoform_strings, "ESPRESSO:$chr:$n:$isoform_id";
          }
          my $isoforms_string = join_with_trailing(',', \@isoform_strings);
          for my $read_details (@{$reads_ref}) {
            my $read_ID = $read_details->{'id'};
            my $ori_id = $read_info_ref->{$read_ID}{'ori_ID'};
            my $sample = $read_info_ref->{$read_ID}{'sample'};
            write_tsv_columns($tsv_compt_handle,
                              [$ori_id, $sample, $compat_class_string,
                               $isoforms_string]);
          }
        }
        $summary_data_ref->{'num_reads_assigned'} += $num_reads;
        if ($should_create_default_handle) {
          write_tsv_columns($default_handle, ['read_count_summary', $chr, $cat,
                                              "$n:$chain", $type, $num_reads]);
        }
        next;
      }

      if (exists $longest_nc_valid_ref->{$chain}) {
        if ($should_create_tsv_compt_handle) {
          my $isoform_id = $longest_nc_valid_ref->{$chain};
          my $isoform_string = "ESPRESSO:$chr:$n:$isoform_id";
          for my $read_details (@{$reads_ref}) {
            my $read_ID = $read_details->{'id'};
            my $ori_id = $read_info_ref->{$read_ID}{'ori_ID'};
            my $sample = $read_info_ref->{$read_ID}{'sample'};
            write_tsv_columns($tsv_compt_handle,
                              [$ori_id, $sample, $compat_class_string,
                               $isoform_string]);
          }
        }
        $summary_data_ref->{'num_reads_assigned'} += $num_reads;
        if ($should_create_default_handle) {
          write_tsv_columns($default_handle,
                            ['read_count_summary', $chr, $cat, "$n:$chain",
                             'itself', $num_reads]);
        }
        next;
      }

      if ($should_create_tsv_compt_handle) {
        for my $read_details (@{$reads_ref}) {
          my $read_ID = $read_details->{'id'};
          my $ori_id = $read_info_ref->{$read_ID}{'ori_ID'};
          my $sample = $read_info_ref->{$read_ID}{'sample'};
          write_tsv_columns($tsv_compt_handle, [$ori_id, $sample,
                                                $compat_class_string, 'NA']);
        }
      }
      $summary_data_ref->{'num_reads_not_assigned'} += $num_reads;
      if ($should_create_default_handle) {
        write_tsv_columns($default_handle, ['read_count_summary', $chr, $cat,
                                            "$n:$chain", 'NA2', $num_reads]);
      }
    }
  }
}

sub assign_reads_to_compatible_ism_chains {
  my ($chr, $n, $group_start, $group_end, $SJ_dist, $internal_boundary_limit,
      $allow_longer_terminal_exons, $should_create_tsv_compt_handle,
      $should_create_default_handle, $cat_SJ_chains_sort_ref, $substring_nc_ref,
      $longest_nc_valid_ref, $nc_valid_longest_ref, $SJ_chain_read_ref,
      $SJ_chain_read_exon_end_ref, $read_info_ref, $multi_exon_isoform_end_ref,
      $anno_SJ_ref, $fsm_chain_ref, $all_SJ_group_sort_ref,
      $ism_mt1_chain_count_ref, $SJ_chain_cat_ref, $read_count_isoform_ref,
      $summary_data_ref, $tsv_compt_handle, $default_handle) = @_;

  my $ism_chain_ref = $SJ_chain_cat_ref->{'ism'};
  for my $ism_SJ_chain (@{$cat_SJ_chains_sort_ref->{'ism'}}) {
    my $chain_read_ref = $SJ_chain_read_ref->{$ism_SJ_chain};
    my $reads_ref = $chain_read_ref->{'reads'};
    my $compat_isoforms_ref = $ism_chain_ref->{$ism_SJ_chain}{'isoforms'};
    my $num_reads = scalar(@{$reads_ref});
    my @ism_SJ_order = split('_', $ism_SJ_chain);
    my $first_sj_i = $ism_SJ_order[0];
    my $last_sj_i = $ism_SJ_order[-1];
    if (exists $longest_nc_valid_ref->{$ism_SJ_chain}) {
      # Reads for all chains have already been assigned to themselves in
      # read_count_isoform (initial and direct).
      # ISM chains have their individual read endpoints checked.
      # The direct assignment count for longest novel ISM chains is updated here.
      my $isoform_id = $longest_nc_valid_ref->{$ism_SJ_chain};
      my $espresso_id = "ESPRESSO:$chr:$n:$isoform_id";
      my $chain_end_ref = $SJ_chain_read_exon_end_ref->{$ism_SJ_chain};
      my $expected_start = $chain_end_ref->{'start'}{$first_sj_i};
      my $is_expected_start_terminal = 1;
      my $expected_end = $chain_end_ref->{'end'}{$last_sj_i};
      my $is_expected_end_terminal = 1;
      for my $read_details (@{$reads_ref}) {
        my $read_ID = $read_details->{'id'};
        my $ori_id = $read_info_ref->{$read_ID}{'ori_ID'};
        my $sample = $read_info_ref->{$read_ID}{'sample'};
        my $read_start = $read_details->{'start'};
        my $read_end = $read_details->{'end'};
        my $actual_start = $read_start;
        my $actual_end = $read_end;
        my $ends_ok = check_endpoint_compatibility(
          $expected_start, $is_expected_start_terminal, $expected_end,
          $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
          $internal_boundary_limit, $allow_longer_terminal_exons,
          $summary_data_ref);
        if ($ends_ok) {
          if ($should_create_tsv_compt_handle) {
            write_tsv_columns($tsv_compt_handle, [$ori_id, $sample, 'NIC',
                                                  "$espresso_id,"]);
          }
          $summary_data_ref->{'num_reads_assigned'} ++;
        } else {
          $read_count_isoform_ref->{$ism_SJ_chain}{$sample}{'direct'} --;
          if ($should_create_tsv_compt_handle) {
            write_tsv_columns($tsv_compt_handle,
                              [$ori_id, $sample, 'ISM', 'NA']);
          }
          $summary_data_ref->{'num_reads_not_assigned'} ++;
        }
      }

      if ($should_create_default_handle) {
        write_tsv_columns($default_handle,
                          ['read_count_summary', $chr, 'ism', "$n:$ism_SJ_chain",
                           'itself', $num_reads]);
      }
      next;
    }

    if (!exists $substring_nc_ref->{$ism_SJ_chain}) {
      while (my ($nc_chain, $nc_chain_ID) = each %{$longest_nc_valid_ref}) {
        if (!exists $chain_read_ref->{'compat_chains'}{$nc_chain}) {
          next;
        }
        my ($comp_result, $num_leading_unpaired_SJs) = comp_SJ_chain(
          $ism_SJ_chain, $nc_chain);
        if ($comp_result == 2) {
          my @SJs_j = split('_', $nc_chain);
          my $first_sj_j = $SJs_j[0];
          my $last_sj_j = $SJs_j[-1];
          my $chain_j_exon_details = $SJ_chain_read_exon_end_ref->{$nc_chain};
          my $chain_j_exon_start = $chain_j_exon_details->{'start'}{$first_sj_i};
          my $chain_j_exon_end = $chain_j_exon_details->{'end'}{$last_sj_i};

          my $expected_start = $chain_j_exon_start;
          my $is_expected_start_terminal = $first_sj_i == $first_sj_j;
          my $expected_end = $chain_j_exon_end;
          my $is_expected_end_terminal = $last_sj_i == $last_sj_j;
          my $actual_start = $chain_read_ref->{'start'};
          my $actual_end = $chain_read_ref->{'end'};
          my $ends_ok = check_endpoint_compatibility(
            $expected_start, $is_expected_start_terminal, $expected_end,
            $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
            $internal_boundary_limit, $allow_longer_terminal_exons,
            $summary_data_ref);
          if ($ends_ok) {
            push @{$compat_isoforms_ref}, $nc_chain_ID;
          }
        }
      }
    }

    my $num_candidate_isoforms = scalar(@{$compat_isoforms_ref});
    if ($num_candidate_isoforms == 0) {
      if ($should_create_tsv_compt_handle) {
        for my $read_details (@{$reads_ref}) {
          my $read_ID = $read_details->{'id'};
          my $ori_id = $read_info_ref->{$read_ID}{'ori_ID'};
          my $sample = $read_info_ref->{$read_ID}{'sample'};
          write_tsv_columns($tsv_compt_handle, [$ori_id, $sample, 'ISM', 'NA']);
        }
      }
      $summary_data_ref->{'num_reads_not_assigned'} += $num_reads;
      if ($should_create_default_handle) {
        write_tsv_columns($default_handle,
                          ['read_count_summary', $chr, 'ism', "$n:$ism_SJ_chain",
                           'NA0', $num_reads]);
      }
      next;
    }

    my $first_sj_string = $all_SJ_group_sort_ref->[$first_sj_i];
    my $last_sj_string = $all_SJ_group_sort_ref->[$last_sj_i];
    my $first_sj_start = $anno_SJ_ref->{$first_sj_string}{'start'};
    my $last_sj_end = $anno_SJ_ref->{$last_sj_string}{'end'};
    my %chains_string_freq = ();
    for my $read_details (@{$reads_ref}) {
      my $read_ID = $read_details->{'id'};
      my $read_start = $read_details->{'start'};
      my $read_end = $read_details->{'end'};
      my $ori_id = $read_info_ref->{$read_ID}{'ori_ID'};
      my $sample = $read_info_ref->{$read_ID}{'sample'};
      my %fsm_SJ_chain_right = ();
      my %read_compatible_isoform = ();
      for my $possible_isoform (@{$compat_isoforms_ref}) {
        if (exists $multi_exon_isoform_end_ref->{$possible_isoform}) {
          my $isoform_end_ref = $multi_exon_isoform_end_ref->{$possible_isoform};
          my $fsm_chain_from_isoform = $fsm_chain_ref->{$possible_isoform};
          # Multiple isoforms could have the same chain.
          if (!exists $fsm_SJ_chain_right{$fsm_chain_from_isoform}) {
            $fsm_SJ_chain_right{$fsm_chain_from_isoform} = 0;
          }
          $read_compatible_isoform{$possible_isoform} = 0;
          if ((!exists $isoform_end_ref->{$first_sj_start})
              or (!exists $isoform_end_ref->{$last_sj_end})) {
            next;
          }

          my $expected_first_start = $isoform_end_ref->{'0'};
          my $expected_last_end = $isoform_end_ref->{'1'};
          my $expected_start = $isoform_end_ref->{$first_sj_start};
          my $is_expected_start_terminal = (
            $expected_start == $expected_first_start);
          my $expected_end = $isoform_end_ref->{$last_sj_end};
          my $is_expected_end_terminal = $expected_end == $expected_last_end;
          my $actual_start = $read_start;
          my $actual_end = $read_end;
          my $ends_ok = check_endpoint_compatibility(
            $expected_start, $is_expected_start_terminal, $expected_end,
            $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
            $internal_boundary_limit, $allow_longer_terminal_exons,
            $summary_data_ref);
          if ($ends_ok) {
            $fsm_SJ_chain_right{$fsm_chain_from_isoform} ++;
            $read_compatible_isoform{$possible_isoform} ++;
          }
        } elsif (($possible_isoform =~ /^\d+$/)
                 and (defined $nc_valid_longest_ref->[$possible_isoform])) {
          my $novel_chain = $nc_valid_longest_ref->[$possible_isoform];
          my $novel_name = "ESPRESSO:$chr:$n:$possible_isoform";
          if (!exists $fsm_SJ_chain_right{$novel_chain}) {
            $fsm_SJ_chain_right{$novel_chain} = 0;
          }
          $read_compatible_isoform{$novel_name} = 0;

          my @SJs_j = split('_', $novel_chain);
          my $first_sj_j = $SJs_j[0];
          my $last_sj_j = $SJs_j[-1];
          my $chain_j_exon_details = $SJ_chain_read_exon_end_ref->{$novel_chain};
          my $chain_j_exon_start = $chain_j_exon_details->{'start'}{$first_sj_i};
          my $chain_j_exon_end = $chain_j_exon_details->{'end'}{$last_sj_i};

          my $expected_start = $chain_j_exon_start;
          my $is_expected_start_terminal = $first_sj_i == $first_sj_j;
          my $expected_end = $chain_j_exon_end;
          my $is_expected_end_terminal = $last_sj_i == $last_sj_j;
          my $actual_start = $read_start;
          my $actual_end = $read_end;
          my $ends_ok = check_endpoint_compatibility(
            $expected_start, $is_expected_start_terminal, $expected_end,
            $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
            $internal_boundary_limit, $allow_longer_terminal_exons,
            $summary_data_ref);
          if ($ends_ok) {
            $fsm_SJ_chain_right{$novel_chain} ++;
            $read_compatible_isoform{$novel_name} ++;
          }
        } else {
          my $die_msg = "don't know what fsm is for $ism_SJ_chain";
          $die_msg = $die_msg . " in group $n($group_start, $group_end)";
          $die_msg = $die_msg . ": $possible_isoform";
          die $die_msg;
        }
      }

      my @valid_SJ_chain = grep {
        $fsm_SJ_chain_right{$_} > 0;
      } keys %fsm_SJ_chain_right;
      if (scalar(@valid_SJ_chain) == 0) {
        if ($should_create_tsv_compt_handle) {
          write_tsv_columns($tsv_compt_handle, [$ori_id, $sample, 'ISM', 'NA']);
        }
        $summary_data_ref->{'num_reads_not_assigned'} ++;
        next;
      }

      # If there was only 1 isoform and it passed the read endpoint check then
      # it can be direct assigned.
      if ($num_candidate_isoforms == 1) {
        my $compat_chain = $valid_SJ_chain[0];
        $read_count_isoform_ref->{$compat_chain}{$sample}{'direct'} ++;
      } else {
        my @sort_chains = sort {$a cmp $b} @valid_SJ_chain;
        my $chains_compt_read_string = join(':', @sort_chains);
        $chains_string_freq{$chains_compt_read_string}{$sample} ++;
      }

      if ($should_create_tsv_compt_handle) {
        my @valid_isoforms = grep {
          $read_compatible_isoform{$_} > 0;
        } keys %read_compatible_isoform;

        my $isoforms_string = join_with_trailing(',', \@valid_isoforms);
        write_tsv_columns($tsv_compt_handle, [$ori_id, $sample, 'ISM',
                                              $isoforms_string]);
      }
      $summary_data_ref->{'num_reads_assigned'} ++;
    }

    my %mt1_for_chain = ();
    my $has_any_mt1 = 0;
    while (my ($chains_string, $sample_ref) = each %chains_string_freq) {
      while (my ($sample, $count) = each %{$sample_ref}) {
        my %details = ();
        $details{'chains'} = $chains_string;
        $details{'count'} = $count;
        push @{$mt1_for_chain{$sample}}, \%details;
        $has_any_mt1 = 1;
      }
    }
    if ($has_any_mt1) {
      $ism_mt1_chain_count_ref->{$ism_SJ_chain} = \%mt1_for_chain;
    }

    if ($should_create_default_handle) {
      my $type = 'EM';
      if ($num_candidate_isoforms == 1) {
        $type = 'direct';
      }
      write_tsv_columns($default_handle,
                        ['read_count_summary', $chr, 'ism', "$n:$ism_SJ_chain",
                         $type, $num_reads]);
    }
  }
}

sub assign_filtered_reads_to_compatible_chains {
  my ($chr, $n, $SJ_dist, $internal_boundary_limit, $allow_longer_terminal_exons,
      $should_create_tsv_compt_handle, $should_create_default_handle,
      $cat_SJ_chains_sort_ref, $SJ_chain_cat_ref, $longest_nc_valid_ref,
      $SJ_chain_read_ref, $SJ_chain_read_exon_end_ref, $read_info_ref,
      $samples_sort_ref, $possible_fsm_filtered_ref, $read_count_isoform_ref,
      $summary_data_ref, $tsv_compt_handle, $default_handle) = @_;

  for my $filtered_SJ_chain (@{$cat_SJ_chains_sort_ref->{'filtered'}}) {
    my @SJs_i = split('_', $filtered_SJ_chain);
    my $filtered_details = $SJ_chain_read_ref->{$filtered_SJ_chain};
    my $filtered_chain_start = $filtered_details->{'start'};
    my $filtered_chain_end = $filtered_details->{'end'};
    my $filtered_compat = $filtered_details->{'compat_chains'};
    my $filtered_reads = $filtered_details->{'reads'};
    my $num_reads = scalar(@{$filtered_reads});
    my @SJ_chains = ();
    my @isoforms_compatible = ();
    for my $fsm_SJ_chain (@{$cat_SJ_chains_sort_ref->{'fsm'}}) {
      if (!exists $filtered_compat->{$fsm_SJ_chain}) {
        next;
      }

      my ($comp_result, $num_leading_unpaired_SJs) = comp_SJ_chain(
        $filtered_SJ_chain, $fsm_SJ_chain);
      if ($comp_result == 3) {
        my @SJs_j = split('_', $fsm_SJ_chain);
        my $first_sj_j = $SJs_j[$num_leading_unpaired_SJs];
        my $last_paired_index = $num_leading_unpaired_SJs + (scalar @SJs_i) - 1;
        my $last_sj_j = $SJs_j[$last_paired_index];
        my $fsm_chain_details = $SJ_chain_read_exon_end_ref->{$fsm_SJ_chain};
        my $expected_start = $fsm_chain_details->{'start'}{$first_sj_j};
        my $is_expected_start_terminal = $num_leading_unpaired_SJs == 0;
        my $expected_end = $fsm_chain_details->{'end'}{$last_sj_j};
        my $is_expected_end_terminal = $last_paired_index == $#SJs_j;
        my $actual_start = $filtered_chain_start;
        my $actual_end = $filtered_chain_end;
        my $ends_ok = check_endpoint_compatibility(
          $expected_start, $is_expected_start_terminal, $expected_end,
          $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
          $internal_boundary_limit, $allow_longer_terminal_exons,
          $summary_data_ref);
        if ($ends_ok) {
          push @SJ_chains, $fsm_SJ_chain;
          my $fsm_details = $SJ_chain_cat_ref->{'fsm'}{$fsm_SJ_chain};
          push @isoforms_compatible, @{$fsm_details->{'isoforms'}};
        }
      }
    }

    while (my ($nc_chain, $nc_chain_ID) = each %{$longest_nc_valid_ref}) {
      if (!exists $filtered_compat->{$nc_chain}) {
        next;
      }

      my ($comp_result, $num_leading_unpaired_SJs) = comp_SJ_chain(
        $filtered_SJ_chain, $nc_chain);
      if ($comp_result == 3) {
        my @SJs_j = split('_', $nc_chain);
        my $first_sj_j = $SJs_j[$num_leading_unpaired_SJs];
        my $last_paired_index = $num_leading_unpaired_SJs + (scalar @SJs_i) - 1;
        my $last_sj_j = $SJs_j[$last_paired_index];
        my $nc_chain_details = $SJ_chain_read_exon_end_ref->{$nc_chain};
        my $expected_start = $nc_chain_details->{'start'}{$first_sj_j};
        my $is_expected_start_terminal = $num_leading_unpaired_SJs == 0;
        my $expected_end = $nc_chain_details->{'end'}{$last_sj_j};
        my $is_expected_end_terminal = $last_paired_index == $#SJs_j;
        my $actual_start = $filtered_chain_start;
        my $actual_end = $filtered_chain_end;
        my $ends_ok = check_endpoint_compatibility(
          $expected_start, $is_expected_start_terminal, $expected_end,
          $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
          $internal_boundary_limit, $allow_longer_terminal_exons,
          $summary_data_ref);
        if ($ends_ok) {
          push @SJ_chains, $nc_chain;
          push @isoforms_compatible, "ESPRESSO:$chr:$n:$nc_chain_ID";
        }
      }
    }

    my $num_matched_chains = scalar(@SJ_chains);
    my $isoforms_string = 'NA';
    my $type = 'NA';
    if ($num_matched_chains >= 1) {
      $isoforms_string = join_with_trailing(',', \@isoforms_compatible);
      $type = 'EM';
      $possible_fsm_filtered_ref->{$num_matched_chains}{$filtered_SJ_chain} = (
        \@SJ_chains);
      $summary_data_ref->{'num_reads_assigned'} += $num_reads;
    } else {
      $summary_data_ref->{'num_reads_not_assigned'} += $num_reads;
    }

    if ($should_create_tsv_compt_handle) {
      for my $read_details (@{$filtered_reads}) {
        my $read_ID = $read_details->{'id'};
        my $ori_id = $read_info_ref->{$read_ID}{'ori_ID'};
        my $sample = $read_info_ref->{$read_ID}{'sample'};
        write_tsv_columns($tsv_compt_handle, [$ori_id, $sample, 'NCD',
                                              $isoforms_string]);
      }
    }
    # ($num_matched_chains == 1) are handled later
    if ($should_create_default_handle and ($num_matched_chains != 1)) {
      write_tsv_columns($default_handle,
                        ['read_count_summary', $chr, 'filtered',
                         "$n:$filtered_SJ_chain", $type, $num_reads]);
    }
  }

  my $one_chain_ref = $possible_fsm_filtered_ref->{'1'};
  while (my ($filtered_SJ_chain, $chain_ref) = each %{$one_chain_ref}) {
    my $filtered_counts_ref = $read_count_isoform_ref->{$filtered_SJ_chain};
    my $compat_chain = $chain_ref->[0];
    for my $sample (@{$samples_sort_ref}) {
      my $initial = $filtered_counts_ref->{$sample}{'initial'};
      $read_count_isoform_ref->{$compat_chain}{$sample}{'direct'} += $initial;
    }

    if ($should_create_default_handle) {
      my $filtered_reads = $SJ_chain_read_ref->{$filtered_SJ_chain}{'reads'};
      my $num_reads = scalar(@{$filtered_reads});
      write_tsv_columns($default_handle,
                        ['read_count_summary', $chr, 'filtered',
                         "$n:$filtered_SJ_chain", 'direct', $num_reads]);
    }
  }
}

sub initialize_em_from_direct_counts {
  my ($read_count_isoform_ref) = @_;

  for my $sample_ref (values %{$read_count_isoform_ref}) {
    for my $count_ref (values %{$sample_ref}) {
      $count_ref->{'EM'} = $count_ref->{'direct'};
    }
  }
}

# The EM algorithm loops until the largest_diff drops below a threshold.
# The outcome of the largest_diff check can change based on the order of
# floating point operations. For example the largest_diff could be calculated as
# 0.09999 instead of 0.1 which would lead to an extra loop iteration.
# Sort the data structures to ensure a consistent result.
sub sort_chains_for_em {
  my ($ism_mt1_chain_count_ref, $possible_fsm_filtered_ref,
      $cat_SJ_chains_sort_ref, $SJ_chain_cat_ref, $sorted_ism_mt1_chains_ref,
      $sorted_samples_by_ism_mt1_chain_ref, $sorted_possible_fsm_num_ref,
      $sorted_chains_by_fsm_num_ref) = @_;

  while (my ($chain, $sample_ref) = each %{$ism_mt1_chain_count_ref}) {
    push @{$sorted_ism_mt1_chains_ref}, $chain;
    my @sorted_samples_for_chain = ();
    while (my ($sample, $chains_count_ref) = each %{$sample_ref}) {
      push @sorted_samples_for_chain, $sample;
      @{$chains_count_ref} = sort {
        $a->{'chains'} cmp $b->{'chains'}
      } @{$chains_count_ref};
    }
    @sorted_samples_for_chain = sort {$a cmp $b} @sorted_samples_for_chain;
    $sorted_samples_by_ism_mt1_chain_ref->{$chain} = \@sorted_samples_for_chain;
  }
  @{$sorted_ism_mt1_chains_ref} = sort {
    $a cmp $b
  } @{$sorted_ism_mt1_chains_ref};

  for my $cat ('nic', 'nnc') {
    my $by_chain = $SJ_chain_cat_ref->{$cat};
    for my $chain (@{$cat_SJ_chains_sort_ref->{$cat}}) {
      my $chain_details = $by_chain->{$chain};
      if ((defined $chain_details->{'isoforms'})
          and (scalar(@{$chain_details->{'isoforms'}}) > 1)) {
        @{$chain_details->{'isoforms'}} = sort {
          $a cmp $b
        } @{$chain_details->{'isoforms'}};
      }
    }
  }

  while (my ($fsm_num, $by_chain_ref) = each %{$possible_fsm_filtered_ref}) {
    push @{$sorted_possible_fsm_num_ref}, $fsm_num;
    my @sorted_chains_for_fsm_num = ();
    while (my ($filtered_chain, $fsm_chain_ref) = each %{$by_chain_ref}) {
      push @sorted_chains_for_fsm_num, $filtered_chain;
      @{$fsm_chain_ref} = sort {$a cmp $b} @{$fsm_chain_ref};
    }
    @sorted_chains_for_fsm_num = sort {$a cmp $b} @sorted_chains_for_fsm_num;
    $sorted_chains_by_fsm_num_ref->{$fsm_num} = (
      \@sorted_chains_for_fsm_num);
  }
  @{$sorted_possible_fsm_num_ref} = sort {
    $a <=> $b
  } @{$sorted_possible_fsm_num_ref};
}

sub round_to_decimal_places {
  my ($value, $places) = @_;

  return sprintf('%.' . "$places" . 'f', $value);
}

sub expectation_maximization {
  my ($chr, $n, $max_iterate, $should_create_default_handle, $samples_sort_ref,
      $cat_SJ_chains_sort_ref, $SJ_chain_cat_ref, $nc_valid_longest_ref,
      $ism_mt1_chain_count_ref, $sorted_ism_mt1_chains_ref,
      $sorted_samples_by_ism_mt1_chain_ref, $possible_fsm_filtered_ref,
      $sorted_possible_fsm_num_ref, $sorted_chains_by_fsm_num_ref,
      $read_count_isoform_ref, $default_handle) = @_;

  my $iteration_threshold = 0.1;
  my $iter_i = 0;
  while (1) {
    $iter_i ++;

    for my $ism_SJ_chain (@{$sorted_ism_mt1_chains_ref}) {
      my $sample_ref = $ism_mt1_chain_count_ref->{$ism_SJ_chain};
      my $total_count = 0;
      for my $sample (@{$sorted_samples_by_ism_mt1_chain_ref->{$ism_SJ_chain}}) {
        my $chains_count_ref = $sample_ref->{$sample};
        for my $chains_string_freq_ref (@{$chains_count_ref}) {
          my $chains_string = $chains_string_freq_ref->{'chains'};
          my $count = $chains_string_freq_ref->{'count'};
          $total_count += $count;
          my @chains = split(':', $chains_string);
          my $num_chains = scalar(@chains);
          my $total_chain_em = 0;
          for my $chain (@chains) {
            my $em_count = $read_count_isoform_ref->{$chain}{$sample}{'EM'};
            $total_chain_em += $em_count;
          }
          if ($total_chain_em > 0) {
            for my $chain (@chains) {
              my $sample_counts = $read_count_isoform_ref->{$chain}{$sample};
              my $ratio = $sample_counts->{'EM'} / $total_chain_em;
              my $abun = $ratio * $count;
              $sample_counts->{'current'} += round_to_decimal_places($abun, 2);
            }
          } else {
            my $equal_split = round_to_decimal_places($count / $num_chains, 2);
            for my $chain (@chains) {
              $read_count_isoform_ref->{$chain}{$sample}{'current'} += (
                $equal_split);
            }
          }
          if ($should_create_default_handle) {
            my $chains_str = '';
            for my $chain (@chains) {
              $chains_str = $chains_str . "$chain($count)"
            }
            write_tsv_columns($default_handle,
                              ["special_ISM_read_assigned:$ism_SJ_chain",
                               $sample, $num_chains, $chains_str]);
          }
        }
      }
      if ($should_create_default_handle) {
        write_tsv_columns($default_handle,
                          ['read_count_summary_in_EM', $chr, 'ISM',
                           "$n:$ism_SJ_chain", 'EM', $total_count, $iter_i]);
      }
    }

    for my $cat ('nic', 'nnc') {
      my $by_chain = $SJ_chain_cat_ref->{$cat};
      for my $chain (@{$cat_SJ_chains_sort_ref->{$cat}}) {
        my $chain_details = $by_chain->{$chain};
        if ((!defined $chain_details->{'isoforms'})
            or (scalar(@{$chain_details->{'isoforms'}}) <= 1)) {
          next;
        }

        my @longest_chains = map {
          $nc_valid_longest_ref->[$_]
        } @{$chain_details->{'isoforms'}};
        my $num_chains = scalar(@longest_chains);
        my $counts_by_sample = $read_count_isoform_ref->{$chain};
        my $total_count = 0;
        for my $sample (@{$samples_sort_ref}) {
          my $initial = $counts_by_sample->{$sample}{'initial'};
          if ($initial == 0) {
            next;
          }

          $total_count += $initial;
          my $total_chain_em = 0;
          for my $longest_chain (@longest_chains) {
            $total_chain_em += (
              $read_count_isoform_ref->{$longest_chain}{$sample}{'EM'});
          }
          if ($total_chain_em > 0) {
            for my $longest_chain (@longest_chains) {
              my $count_details = (
                $read_count_isoform_ref->{$longest_chain}{$sample});
              my $ratio = $count_details->{'EM'} / $total_chain_em;
              my $abun = $ratio * $initial;
              $count_details->{'current'} += round_to_decimal_places($abun, 2);
            }
          } else {
            my $equal_split = round_to_decimal_places($initial / $num_chains, 2);
            for my $longest_chain (@longest_chains) {
              $read_count_isoform_ref->{$longest_chain}{$sample}{'current'} += (
                $equal_split);
            }
          }
        }
        if ($should_create_default_handle) {
          write_tsv_columns($default_handle,
                            ['read_count_summary_in_EM', $chr, $cat, "$n:$chain",
                             'EM', $total_count, $iter_i]);
        }
      }
    }

    for my $fsm_num (@{$sorted_possible_fsm_num_ref}) {
      if ($fsm_num <= 1) {
        next;
      }

      my $isoforms_by_chain = $possible_fsm_filtered_ref->{$fsm_num};
      for my $filtered_chain (@{$sorted_chains_by_fsm_num_ref->{$fsm_num}}) {
        my $fsm_chains = $isoforms_by_chain->{$filtered_chain};
        my $num_chains = scalar(@{$fsm_chains});
        my $counts_by_sample = $read_count_isoform_ref->{$filtered_chain};
        my $total_count = 0;
        for my $sample (@{$samples_sort_ref}) {
          my $sample_counts = $counts_by_sample->{$sample};
          my $initial = $sample_counts->{'initial'};
          if ($initial == 0) {
            next;
          }
          $total_count += $initial;
          my $total_chain_em = 0;
          for my $chain (@{$fsm_chains}) {
            $total_chain_em += $read_count_isoform_ref->{$chain}{$sample}{'EM'};
          }
          if ($total_chain_em > 0) {
            for my $chain (@{$fsm_chains}) {
              my $count_details = $read_count_isoform_ref->{$chain}{$sample};
              my $ratio = $count_details->{'EM'} / $total_chain_em;
              my $abun = $ratio * $initial;
              $count_details->{'current'} += round_to_decimal_places($abun, 2);
            }
          } else {
            my $equal_split = round_to_decimal_places($initial / $num_chains, 2);
            for my $chain (@{$fsm_chains}) {
              $read_count_isoform_ref->{$chain}{$sample}{'current'} += (
                $equal_split);
            }
          }
        }
        if ($should_create_default_handle) {
          write_tsv_columns($default_handle,
                            ['read_count_summary_in_EM', $chr, 'filtered',
                             "$n:$filtered_chain", 'EM', $total_count, $iter_i]);
        }
      }
    }

    my $largest_diff = 0;
    for my $chain (keys %{$read_count_isoform_ref}) {
      my $counts_by_sample = $read_count_isoform_ref->{$chain};
      for my $sample (@{$samples_sort_ref}) {
        my $sample_counts = $counts_by_sample->{$sample};
        $sample_counts->{'EM'} = ($sample_counts->{'direct'}
                                  + $sample_counts->{'current'});
        my $abs_diff = abs($sample_counts->{'current'}
                           - $sample_counts->{'previous'});
        if ($largest_diff < $abs_diff) {
          $largest_diff = $abs_diff;
        }
        $sample_counts->{'previous'} = $sample_counts->{'current'};
        $sample_counts->{'current'} = 0;
      }
    }

    if ($largest_diff < $iteration_threshold) {
      if ($should_create_default_handle) {
        write_tsv_columns($default_handle, ["converge: $n", $largest_diff,
                                            $iter_i]);
      }
      last;
    } elsif ($iter_i >= $max_iterate) {
      if ($should_create_default_handle) {
        write_tsv_columns($default_handle, ["not_converge: $n", $largest_diff,
                                            $iter_i]);
      }
      last;
    }
  }
}

sub output_fsm_abundance {
  my ($chr, $n, $should_create_tsv_compt_handle, $should_create_default_handle,
      $strand_symbol_ref, $SJ_chain_cat_fsm_ref, $SJ_chain_read_ref,
      $all_SJ_group_sort_ref, $read_count_isoform_ref, $read_info_ref,
      $isoform_perfect_match_ref, $multi_exon_isoform_end_ref, $isoform_info_ref,
      $samples_sort_ref, $summary_data_ref, $abu_handle, $gtf_handle,
      $tsv_compt_handle, $default_handle) = @_;

  while (my ($fsm_chain, $chain_info_ref) = each %{$SJ_chain_cat_fsm_ref}) {
    my $chain_details = $SJ_chain_read_ref->{$fsm_chain};
    my $chain_isoforms = $chain_info_ref->{'isoforms'};
    my $num_isoforms = scalar(@{$chain_isoforms});
    my $isSameStrand = $chain_details->{'strand_isoform'};
    my @SJ_IDs = split('_', $fsm_chain);
    my @SJs = map {$all_SJ_group_sort_ref->[$_]} @SJ_IDs;
    my $counts_by_sample = $read_count_isoform_ref->{$fsm_chain};
    if ($should_create_default_handle) {
      my $initial_string = '';
      my $direct_string = '';
      my $em_string = '';
      # trailing commas are used here
      for my $sample (@{$samples_sort_ref}) {
        $initial_string = ($initial_string
                           . $counts_by_sample->{$sample}{'initial'} . ',');
        my $rounded_direct = round_to_decimal_places(
          $counts_by_sample->{$sample}{'direct'}, 2);
        $direct_string = $direct_string . $rounded_direct . ',';
        my $rounded_em = round_to_decimal_places(
          $counts_by_sample->{$sample}{'EM'}, 2);
        $em_string = $em_string . $rounded_em . ',';
      }
      my $isoforms_string = join_with_trailing(',', $chain_isoforms);
      my $coord_string = "$chain_details->{'start'},";
      $coord_string = $coord_string . join(',', @SJs);
      $coord_string = $coord_string . ",$chain_details->{'end'}";
      write_tsv_columns($default_handle,
                        ['assign_reads', $fsm_chain, $initial_string,
                         $direct_string, $em_string, $isoforms_string,
                         $coord_string]);
    }

    my %total_perfect_match = ();
    for my $isoform (@{$chain_isoforms}) {
      for my $sample (@{$samples_sort_ref}) {
        if ((exists $isoform_perfect_match_ref->{$isoform})
            and (exists $isoform_perfect_match_ref->{$isoform}{$sample})) {
          my $num_reads = scalar(
            @{$isoform_perfect_match_ref->{$isoform}{$sample}});
          $total_perfect_match{$sample} += $num_reads;
        }
      }
    }

    for my $isoform (@{$chain_isoforms}) {
      my $isoform_info = $isoform_info_ref->{$isoform};
      my $gene_ID = $isoform_info->{'gene_id'};
      my $isoform_name = $isoform_info->{'isoform_name'};
      my $has_abundance = 0;
      my @abu_columns = ($isoform, $isoform_name, $gene_ID);
      for my $sample (@{$samples_sort_ref}) {
        my $em_counts = $counts_by_sample->{$sample}{'EM'};
        if ((exists $total_perfect_match{$sample})
            and ($total_perfect_match{$sample} > 0)) {
          if ((exists $isoform_perfect_match_ref->{$isoform})
              and (exists $isoform_perfect_match_ref->{$isoform}{$sample})) {
            my $num_perfect = scalar(
              @{$isoform_perfect_match_ref->{$isoform}{$sample}});
            my $ratio = $num_perfect / $total_perfect_match{$sample};
            my $abun = $ratio * $em_counts;
            $abun = round_to_decimal_places($abun, 2);
            push @abu_columns, $abun;
            $summary_data_ref->{'total_fsm_abundance'} += $abun;
            $has_abundance = 1;
          } else {
            push @abu_columns, 0;
          }
        } else {
          if ($em_counts > 0) {
            my $equal_split = $em_counts / $num_isoforms;
            my $abu_amount = round_to_decimal_places($equal_split, 2);
            push @abu_columns, $abu_amount;
            $summary_data_ref->{'total_fsm_abundance'} += $abu_amount;
            $has_abundance = 1;
          } else {
            push @abu_columns, 0;
          }
        }
      }
      if ($has_abundance == 1) {
        $summary_data_ref->{'num_fsm_isoforms_detected'} ++;
        write_tsv_columns($abu_handle, \@abu_columns);
        my $start = $multi_exon_isoform_end_ref->{$isoform}{'0'};
        my $end = $multi_exon_isoform_end_ref->{$isoform}{'1'};
        gtf_output($gtf_handle, $chr, $isoform, $gene_ID, [$start, $end], \@SJs,
                   $strand_symbol_ref->{$isSameStrand}, 'annotated_isoform');
      }
    }

    if ($should_create_tsv_compt_handle) {
      for my $read_details (@{$chain_details->{'reads'}}) {
        my $read_ID = $read_details->{'id'};
        my $ori_id = $read_info_ref->{$read_ID}{'ori_ID'};
        my $sample = $read_info_ref->{$read_ID}{'sample'};
        my $isoforms_string = join_with_trailing(',', $chain_isoforms);
        write_tsv_columns($tsv_compt_handle, [$ori_id, $sample, 'FSM',
                                              $isoforms_string]);
      }
    }
    my $num_reads = scalar(@{$chain_details->{'reads'}});
    $summary_data_ref->{'num_reads_assigned'} += $num_reads;
    if ($should_create_default_handle) {
      write_tsv_columns($default_handle,
                        ['read_count_summary', $chr, 'fsm', "$n:$fsm_chain",
                         'itself', $num_reads]);
    }
  }
}

sub output_novel_abundance {
  my ($chr, $n, $should_create_default_handle, $strand_symbol_ref,
      $longest_nc_valid_ref, $SJ_chain_cat_ref, $SJ_chain_read_ref,
      $all_SJ_group_sort_ref, $read_count_isoform_ref, $unannotated_isoforms_ref,
      $isoform_info_ref, $anno_SJ_ref, $anno_SJ_complete_ref, $samples_sort_ref,
      $summary_data_ref, $abu_handle, $gtf_handle, $default_handle) = @_;

  while (my ($nc_chain, $nc_chain_ID) = each %{$longest_nc_valid_ref}) {
    # longest novel should be one of ism, nic, nnc
    my $chain_cat = '';
    if (exists $SJ_chain_cat_ref->{'ism'}{$nc_chain}) {
      $chain_cat = 'novel_ism';
    } elsif (exists $SJ_chain_cat_ref->{'nic'}{$nc_chain}) {
      $chain_cat = 'nic';
    } elsif (exists $SJ_chain_cat_ref->{'nnc'}{$nc_chain}) {
      $chain_cat = 'nnc';
    }
    my $total_abun_key = "total_$chain_cat"."_abundance";
    my $num_isoforms_key = "num_$chain_cat"."_isoforms_detected";

    my $chain_details = $SJ_chain_read_ref->{$nc_chain};
    my $counts_by_sample = $read_count_isoform_ref->{$nc_chain};
    my $isSameStrand = $chain_details->{'strand_isoform'};
    my @SJ_IDs = split('_', $nc_chain);
    my @SJs = map {$all_SJ_group_sort_ref->[$_]} @SJ_IDs;
    if ($should_create_default_handle) {
      my $initial_string = '';
      my $direct_string = '';
      my $em_string = '';
      # trailing commas are used here
      for my $sample (@{$samples_sort_ref}) {
        $initial_string = ($initial_string
                           . $counts_by_sample->{$sample}{'initial'} . ',');
        my $rounded_direct = round_to_decimal_places(
          $counts_by_sample->{$sample}{'direct'}, 2);
        $direct_string = $direct_string . $rounded_direct . ',';
        my $rounded_em = round_to_decimal_places(
          $counts_by_sample->{$sample}{'EM'}, 2);
        $em_string = $em_string . $rounded_em . ',';
      }

      my $id_str = "$n:$nc_chain_ID";
      if ((exists $unannotated_isoforms_ref->{$nc_chain})
          and (scalar(@{$unannotated_isoforms_ref->{$nc_chain}{'fsm'}}) > 0)) {
        $id_str = join_with_trailing(
          ',', $unannotated_isoforms_ref->{$nc_chain}{'fsm'});
      }

      my $perfect_count = 0;
      my $num_reads = 0;
      if (exists $unannotated_isoforms_ref->{$nc_chain}) {
        my $unanno_details = $unannotated_isoforms_ref->{$nc_chain};
        $perfect_count = $unanno_details->{'lowest_perfect_count'};
        $num_reads = $unanno_details->{'num_reads'};
      }
      my $coord_string = "$chain_details->{'start'},";
      $coord_string = $coord_string . join(',', @SJs);
      $coord_string = $coord_string . ",$chain_details->{'end'}";
      write_tsv_columns($default_handle,
                        ['assign_reads', $nc_chain, $initial_string,
                         $direct_string, $em_string, $id_str, $perfect_count,
                         $num_reads, $coord_string]);
    }

    my $type;
    my $isoform_ID;
    my %possible_genes = ();
    my %possible_isoform_names = ();
    if ((exists $unannotated_isoforms_ref->{$nc_chain})
        and (@{$unannotated_isoforms_ref->{$nc_chain}{'fsm'}} > 0)) {
      $type = 'annotated_isoform';
      my $isoforms_ref = $unannotated_isoforms_ref->{$nc_chain}{'fsm'};
      $isoform_ID = join('/', @{$isoforms_ref});
      for my $isoform (@{$isoforms_ref}) {
        my $isoform_details = $isoform_info_ref->{$isoform};
        $possible_genes{$isoform_details->{'gene_id'}} ++;
        $possible_isoform_names{$isoform_details->{'isoform_name'}} ++;
      }
    } else {
      $type = 'novel_isoform';
      $isoform_ID = "ESPRESSO:$chr:$n:$nc_chain_ID";
      for my $SJ (@SJs) {
        if (exists $anno_SJ_ref->{$SJ}) {
          for my $isoform (@{$anno_SJ_ref->{$SJ}{'isoforms'}}) {
            $possible_genes{$isoform_info_ref->{$isoform}{'gene_id'}} ++;
          }
        } elsif (exists $anno_SJ_complete_ref->{$SJ}) {
          for my $isoform (@{$anno_SJ_complete_ref->{$SJ}{'isoforms'}}) {
            $possible_genes{$isoform_info_ref->{$isoform}{'gene_id'}} ++;
          }
        }
      }
    }

    my @possible_genes_sort = sort {
      $possible_genes{$b} <=> $possible_genes{$a}
    } keys %possible_genes;
    my @possible_isoform_names_sort = sort {
      $possible_isoform_names{$b} <=> $possible_isoform_names{$a}
    } keys %possible_isoform_names;

    my $isoform_name = 'NA';
    if (scalar(@possible_isoform_names_sort) > 0) {
      $isoform_name = join(',', @possible_isoform_names_sort);
    }
    my $possible_genes_ID = 'NA';
    if (scalar(@possible_genes_sort) > 0) {
      $possible_genes_ID = join(',', @possible_genes_sort);
    }

    my @abu_columns = ($isoform_ID, $isoform_name, $possible_genes_ID);
    for my $sample (@{$samples_sort_ref}) {
      my $abun = $counts_by_sample->{$sample}{'EM'};
      $abun = round_to_decimal_places($abun, 2);
      push @abu_columns, $abun;
      $summary_data_ref->{$total_abun_key} += $abun;
    }
    write_tsv_columns($abu_handle, \@abu_columns);
    $summary_data_ref->{$num_isoforms_key} ++;

    # gtf_output expects a 0-based start coordinate
    my $zero_based_start = $chain_details->{'start'} - 1;
    gtf_output($gtf_handle, $chr, $isoform_ID, $possible_genes_ID,
               [$zero_based_start, $chain_details->{'end'}], \@SJs,
               $strand_symbol_ref->{$isSameStrand}, $type);
  }
}

sub sort_single_exon_isoforms_and_reads {
  my ($n, $single_exon_isoforms_by_group_ref, $single_exon_reads_ref,
      $single_exon_isoform_end_ref, $read_info_ref,
      $single_exon_isoforms_group_sort_ref, $single_exon_reads_sort_ref) = @_;

  if (exists $single_exon_isoforms_by_group_ref->{$n}) {
    @{$single_exon_isoforms_group_sort_ref} = sort {
      my $a_info = $single_exon_isoform_end_ref->{$a};
      my $a_start = $a_info->{'0'};
      my $a_end = $a_info->{'1'};
      my $b_info = $single_exon_isoform_end_ref->{$b};
      my $b_start = $b_info->{'0'};
      my $b_end = $b_info->{'1'};
      ($a_start <=> $b_start) or ($a_end <=> $b_end);
    } @{$single_exon_isoforms_by_group_ref->{$n}};
  }

  @{$single_exon_reads_sort_ref} = sort {
    my $a_info = $read_info_ref->{$a};
    my $a_start = $a_info->{'start'};
    my $a_end = $a_info->{'end'};
    my $b_info = $read_info_ref->{$b};
    my $b_start = $b_info->{'start'};
    my $b_end = $b_info->{'end'};
    ($a_start <=> $b_start) or ($a_end <=> $b_end);
  } @{$single_exon_reads_ref};
}

sub initialize_single_exon_assignment_counts {
  my ($single_exon_isoforms_group_sort_ref, $samples_sort_ref,
      $single_exon_isoform_read_count_ref) = @_;
  for my $isoform (@{$single_exon_isoforms_group_sort_ref}) {
    for my $sample (@{$samples_sort_ref}) {
      my %details = ();
      $details{'direct'} = 0;
      $details{'current'} = 0;
      $details{'previous'} = 0;
      $single_exon_isoform_read_count_ref->{$isoform}{$sample} = \%details;
    }
  }
}

sub assign_reads_to_single_exon_isoforms {
  my ($chr, $n, $should_create_tsv_compt_handle, $should_create_default_handle,
      $SJ_dist, $internal_boundary_limit, $allow_longer_terminal_exons,
      $single_exon_reads_sort_ref, $single_exon_isoforms_group_sort_ref,
      $single_exon_isoform_end_ref, $single_exon_reads_ref, $read_info_ref,
      $possible_single_exon_isoform_ref, $single_exon_isoform_read_count_ref,
      $summary_data_ref, $tsv_compt_handle, $default_handle) = @_;

  my $total_assigned = 0;
  my $total_not_assigned = 0;
  # All @single_exon_reads_sort have their endpoints checked against
  # @single_exon_isoforms_group_sort. Since both lists are sorted,
  # keep track of which isoform to start from for each read.
  my $num_isoforms = scalar(@{$single_exon_isoforms_group_sort_ref});
  my $min_isoform_i = 0;
  for my $read_ID (@{$single_exon_reads_sort_ref}) {
    my $read_details = $read_info_ref->{$read_ID};
    my $ori_id = $read_details->{'ori_ID'};
    my $sample = $read_details->{'sample'};
    my $read_start = $read_details->{'start'};
    my $read_end = $read_details->{'end'};
    while ($min_isoform_i < $num_isoforms) {
      my $isoform_id = $single_exon_isoforms_group_sort_ref->[$min_isoform_i];
      my $isoform_info = $single_exon_isoform_end_ref->{$isoform_id};
      my $isoform_end = $isoform_info->{'1'};
      if ($read_start > ($isoform_end + $SJ_dist)) {
        # All remaining reads have a start past this isoform end.
        $min_isoform_i += 1;
        next;
      }
      last;
    }

    for my $isoform_i ($min_isoform_i .. ($num_isoforms - 1)) {
      my $isoform = $single_exon_isoforms_group_sort_ref->[$isoform_i];
      my $expected_start = $single_exon_isoform_end_ref->{$isoform}{'0'};
      my $is_expected_start_terminal = 1;
      my $expected_end = $single_exon_isoform_end_ref->{$isoform}{'1'};
      my $is_expected_end_terminal = 1;
      my $actual_start = $read_start;
      my $actual_end = $read_end;
      my $ends_ok = check_endpoint_compatibility(
        $expected_start, $is_expected_start_terminal, $expected_end,
        $is_expected_end_terminal, $actual_start, $actual_end, $SJ_dist,
        $internal_boundary_limit, $allow_longer_terminal_exons,
        $summary_data_ref);
      if ($ends_ok) {
        push @{$possible_single_exon_isoform_ref->{$read_ID}}, $isoform;
      }
      if ($read_end < ($expected_start - $SJ_dist)) {
        # All remaining isoforms have a start past this read end.
        last;
      }
    }

    if (exists $possible_single_exon_isoform_ref->{$read_ID}) {
      $total_assigned ++;
      my $isoforms_ref = $possible_single_exon_isoform_ref->{$read_ID};
      my $num_compat_isoforms = scalar(@{$isoforms_ref});
      if ($should_create_tsv_compt_handle) {
        my $isoforms_string = join_with_trailing(',', $isoforms_ref);
        write_tsv_columns($tsv_compt_handle, [$ori_id, $sample, 'single-exon',
                                              $isoforms_string]);
      }

      for my $isoform (@{$isoforms_ref}) {
        my $equal_split = 1 / $num_compat_isoforms;
        $single_exon_isoform_read_count_ref->{$isoform}{$sample}{'direct'} += (
          round_to_decimal_places($equal_split, 2));
      }

      $summary_data_ref->{'num_reads_assigned'} ++;
    } else {
      $total_not_assigned ++;
      if ($should_create_tsv_compt_handle) {
        write_tsv_columns($tsv_compt_handle, [$ori_id, $sample, 'single-exon',
                                              'NA']);
      }
      $summary_data_ref->{'num_reads_not_assigned'} ++;
    }
  }

  if ($should_create_default_handle) {
    write_tsv_columns($default_handle,
                      ['read_count_summary', $chr, 'one_exon', "$n:total", 'EM',
                       $total_assigned]);
    write_tsv_columns($default_handle,
                      ['read_count_summary', $chr, 'one_exon', "$n:total", 'NA',
                       $total_not_assigned]);
    my $num_reads = scalar(@{$single_exon_reads_ref});
    write_tsv_columns($default_handle, ['read_count_summary_one_exon', $chr, $n,
                                        $num_reads]);
  }
}

sub initialize_em_from_single_exon_counts {
  my ($single_exon_isoforms_group_sort_ref, $samples_sort_ref,
      $single_exon_isoform_read_count_ref) = @_;

  for my $isoform (@{$single_exon_isoforms_group_sort_ref}) {
    for my $sample (@{$samples_sort_ref}) {
      my $counts = $single_exon_isoform_read_count_ref->{$isoform}{$sample};
      $counts->{'previous'} = $counts->{'direct'};
    }
  }
}

sub sort_single_exon_isoforms_for_em {
  my ($possible_single_exon_isoform_ref,
      $sorted_possible_single_exon_reads_ref) = @_;

  # Sort to ensure consistent result of EM algorithm
  while (my ($read_ID, $isoforms) = each %{$possible_single_exon_isoform_ref}) {
    push @{$sorted_possible_single_exon_reads_ref}, $read_ID;
    @{$isoforms} = sort {$a cmp $b} @{$isoforms};
  }
  @{$sorted_possible_single_exon_reads_ref} = sort {
    $a cmp $b;
  } @{$sorted_possible_single_exon_reads_ref};
}

sub single_exon_expectation_maximization {
  my ($chr, $n, $max_iterate, $should_create_default_handle,
      $sorted_possible_single_exon_reads_ref,
      $single_exon_isoforms_group_sort_ref, $possible_single_exon_isoform_ref,
      $read_info_ref, $samples_sort_ref, $single_exon_isoform_read_count_ref,
      $default_handle) = @_;

  my $iteration_threshold = 0.1;
  my $iter_i = 0;
  while (1) {
    $iter_i ++;

    for my $read_ID (@{$sorted_possible_single_exon_reads_ref}) {
      my $isoforms_ref = $possible_single_exon_isoform_ref->{$read_ID};
      my $sample = $read_info_ref->{$read_ID}{'sample'};
      my $total_em = 0;
      for my $isoform (@{$isoforms_ref}) {
        $total_em += (
          $single_exon_isoform_read_count_ref->{$isoform}{$sample}{'previous'});
      }
      if ($total_em > 0) {
        for my $isoform (@{$isoforms_ref}) {
          my $sample_counts = (
            $single_exon_isoform_read_count_ref->{$isoform}{$sample});
          my $ratio = $sample_counts->{'previous'} / $total_em;
          $sample_counts->{'current'} += round_to_decimal_places($ratio, 2);
        }
      } else {
        my $num_isoforms = scalar(@{$isoforms_ref});
        for my $isoform (@{$isoforms_ref}) {
          my $sample_counts = (
            $single_exon_isoform_read_count_ref->{$isoform}{$sample});
          my $equal_split = 1 / $num_isoforms;
          $sample_counts->{'current'} += round_to_decimal_places(
            $equal_split, 2);
        }
      }
    }

    if ($should_create_default_handle) {
      my $num_reads = scalar(keys %{$possible_single_exon_isoform_ref});
      write_tsv_columns($default_handle,
                        ['read_count_summary_in_EM', $chr, 'one_exon',
                         "$n:total", 'EM', $num_reads, $iter_i]);
    }

    my $largest_diff = 0;
    for my $isoform (@{$single_exon_isoforms_group_sort_ref}) {
      my $counts_by_sample = $single_exon_isoform_read_count_ref->{$isoform};
      for my $sample (@{$samples_sort_ref}) {
        my $sample_counts = $counts_by_sample->{$sample};
        my $abs_diff = abs($sample_counts->{'current'}
                           - $sample_counts->{'previous'});
        if ($largest_diff < $abs_diff) {
          $largest_diff = $abs_diff ;
        }
        $sample_counts->{'previous'} = $sample_counts->{'current'};
        $sample_counts->{'current'} = 0;
      }
    }

    if ($largest_diff < $iteration_threshold) {
      if ($should_create_default_handle) {
        write_tsv_columns($default_handle, ["single_converge: $n", $largest_diff,
                                            $iter_i]);
      }
      last;
    } elsif ($iter_i >= $max_iterate) {
      if ($should_create_default_handle) {
        write_tsv_columns($default_handle, ["single_not_converge: $n",
                                            $largest_diff, $iter_i]);
      }
      last;
    }
  }
}

sub update_single_exon_totals {
  my ($should_create_default_handle, $single_exon_isoform_read_count_ref,
      $single_exon_isoform_end_ref, $samples_sort_ref,
      $single_exon_isoform_total_ref, $default_handle) = @_;

  while (my ($isoform, $counts) = each %{$single_exon_isoform_read_count_ref}) {
    my $total_em = 0;
    for my $sample (@{$samples_sort_ref}) {
      $total_em += $counts->{$sample}{'previous'};
    }
    if ($total_em > 0) {
      if ($should_create_default_handle) {
        my $direct_string = '';
        my $em_string = '';
        # trailing commas are used here
        for my $sample (@{$samples_sort_ref}) {
          $direct_string = $direct_string . $counts->{$sample}{'direct'} . ',';
          my $rounded_em = round_to_decimal_places(
            $counts->{$sample}{'previous'}, 2);
          $em_string = $em_string . $rounded_em . ',';
        }
        my $isoform_string = "$isoform,";
        my $start = $single_exon_isoform_end_ref->{$isoform}{'0'};
        my $end = $single_exon_isoform_end_ref->{$isoform}{'1'};
        my $coord_string = "$start,$end";
        write_tsv_columns($default_handle,
                          ['assign_reads', 'NA', 'NA', $direct_string,
                           $em_string, $isoform_string, $coord_string]);
      }
      if (!exists $single_exon_isoform_total_ref->{$isoform}) {
        my %total_by_sample = ();
        $single_exon_isoform_total_ref->{$isoform} = \%total_by_sample;
      }
      my $isoform_total = $single_exon_isoform_total_ref->{$isoform};
      for my $sample (@{$samples_sort_ref}) {
        $isoform_total->{$sample} += $counts->{$sample}{'previous'};
      }
    }
  }
}

sub output_single_exon_abundance {
  my ($chr, $single_exon_isoform_total_ref, $isoform_info_ref,
      $single_exon_isoform_end_ref, $samples_sort_ref, $summary_data_ref,
      $abu_handle, $gtf_handle) = @_;

  while (my ($isoform, $sample_ref) = each %{$single_exon_isoform_total_ref}) {
    my $info_ref = $isoform_info_ref->{$isoform};
    my $gene_ID = $info_ref->{'gene_id'};
    my $isoform_name = $info_ref->{'isoform_name'};
    my @abu_columns = ($isoform, $isoform_name, $gene_ID);
    for my $sample (@{$samples_sort_ref}) {
      my $abun = round_to_decimal_places($sample_ref->{$sample}, 2);
      push @abu_columns, $abun;
      $summary_data_ref->{'total_single_exon_abundance'} += $abun;
    }
    write_tsv_columns($abu_handle, \@abu_columns);
    $summary_data_ref->{'num_single_exon_isoforms_detected'} ++;

    my $end_ref = $single_exon_isoform_end_ref->{$isoform};
    my $start = $end_ref->{'0'};
    my $end = $end_ref->{'1'};
    my $strand = $end_ref->{'strand'};
    gtf_output($gtf_handle, $chr, $isoform, $gene_ID, [$start, $end], [],
               $strand, 'annotated_isoform');
  }
}

sub output_totals_to_default {
  my ($should_create_default_handle, $chr_type_ref, $default_handle) = @_;

  if (!$should_create_default_handle) {
    return;
  }

  my @types = ('fsm', 'ism', 'nic', 'nnc');
  my %total_type = ();
  while (my ($chr, $type_ref) = each %{$chr_type_ref}) {
    print $default_handle "$chr";
    for my $type (@types) {
      if (exists $type_ref->{$type}) {
        my $count = $type_ref->{$type};
        print $default_handle " $count";
        $total_type{$type} += $count;
      } else {
        print $default_handle " NA";
      }
    }
    print $default_handle "\n";
  }

  for my $type (@types) {
    if (exists $total_type{$type}) {
      print $default_handle " $total_type{$type}";
    } else {
      print $default_handle " NA";
    }
  }
  print $default_handle "\n";
}

sub main {
  my ($stored_shared_arguments, $stored_chr_arguments) = @_;
  my $shared_arguments_ref = retrieve($stored_shared_arguments);
  my $chr_arguments_ref = retrieve($stored_chr_arguments);

  my $out_dir = $shared_arguments_ref->{'out_dir'};
  my $samples_sort_ref = $shared_arguments_ref->{'samples_sort'};
  my $should_create_default_handle = (
    $shared_arguments_ref->{'should_create_default_handle'});
  my $should_create_tsv_compt_handle = (
    $shared_arguments_ref->{'should_create_tsv_compt_handle'});
  my $keep_tmp = $shared_arguments_ref->{'keep_tmp'};
  my $target_col_index = $shared_arguments_ref->{'target_col_index'};
  my $SJ_dist = $shared_arguments_ref->{'SJ_dist'};
  my $internal_boundary_limit = (
    $shared_arguments_ref->{'internal_boundary_limit'});
  my $allow_longer_terminal_exons = (
    $shared_arguments_ref->{'allow_longer_terminal_exons'});
  my $raw = $shared_arguments_ref->{'raw'};
  my $read_num_cutoff = $shared_arguments_ref->{'read_num_cutoff'};
  my $read_ratio_cutoff = $shared_arguments_ref->{'read_ratio_cutoff'};
  my $max_iterate = $shared_arguments_ref->{'max_iterate'};
  my $sort_buffer_size = $shared_arguments_ref->{'sort_buffer_size'};

  my $chr = $chr_arguments_ref->{'chr'};
  my $anno_SJ_ref = $chr_arguments_ref->{'anno_SJ'};
  my $has_tmp_sj = $chr_arguments_ref->{'has_tmp_sj'};
  my $single_exon_isoform_end_ref = (
    $chr_arguments_ref->{'single_exon_isoform_end'});
  my $isoform_SJ_ref = $chr_arguments_ref->{'isoform_SJ'};
  my $anno_SJ_complete_ref = $chr_arguments_ref->{'anno_SJ_complete'};
  my $isoform_SJ_complete_ref = $chr_arguments_ref->{'isoform_SJ_complete'};
  my $isoform_info_ref = $chr_arguments_ref->{'isoform_info'};
  my $anno_SS_ref = $chr_arguments_ref->{'anno_SS'};
  my $multi_exon_isoform_end_ref = (
    $chr_arguments_ref->{'multi_exon_isoform_end'});
  my $read_final_paths_ref = $chr_arguments_ref->{'read_final_paths'};

  process_results_for_chr(
    $chr, $out_dir, $anno_SJ_ref, $has_tmp_sj, $single_exon_isoform_end_ref,
    $anno_SS_ref, $isoform_SJ_ref, $multi_exon_isoform_end_ref,
    $anno_SJ_complete_ref, $isoform_SJ_complete_ref, $samples_sort_ref,
    $isoform_info_ref, $read_final_paths_ref, $should_create_default_handle,
    $should_create_tsv_compt_handle, $keep_tmp, $target_col_index, $SJ_dist,
    $internal_boundary_limit, $allow_longer_terminal_exons, $raw,
    $read_num_cutoff, $read_ratio_cutoff, $max_iterate, $sort_buffer_size);
}

sub process_results_for_chr {
  my ($chr, $out_dir, $anno_SJ_ref, $has_tmp_sj, $single_exon_isoform_end_ref,
      $anno_SS_ref, $isoform_SJ_ref, $multi_exon_isoform_end_ref,
      $anno_SJ_complete_ref, $isoform_SJ_complete_ref, $samples_sort_ref,
      $isoform_info_ref, $read_final_paths_ref, $should_create_default_handle,
      $should_create_tsv_compt_handle, $keep_tmp, $target_col_index, $SJ_dist,
      $internal_boundary_limit, $allow_longer_terminal_exons, $raw,
      $read_num_cutoff, $read_ratio_cutoff, $max_iterate,
      $sort_buffer_size) = @_;

  my $default_handle;
  if ($should_create_default_handle) {
    my $default_path = default_path_for_chr($chr, $out_dir);
    open($default_handle, '>', $default_path)
      or die "cannot write $default_path: $!";
  }

  my $tsv_compt_handle;
  if ($should_create_tsv_compt_handle) {
    my $tsv_compt_path = tsv_compt_path_for_chr($chr, $out_dir);
    open($tsv_compt_handle, '>', $tsv_compt_path)
      or die "cannot write $tsv_compt_path: $!";
  }

  my $abu_path = abu_path_for_chr($chr, $out_dir);
  open(my $abu_handle, '>', $abu_path) or die "cannot write $abu_path: $!";
  my $gtf_path = gtf_path_for_chr($chr, $out_dir);
  open(my $gtf_handle, '>', $gtf_path) or die "cannot write $gtf_path: $!";

  my %unanno_SJ = ();
  my %group_ends = ();
  my @group_sort = ();
  my %sj_to_groups = ();
  my %summary_data = ();
  my $group_info_path = "$out_dir/$chr.group_info.tmp";
  my $group_info_path_unsorted = $group_info_path . '.unsorted';
  load_read_final($out_dir, $chr, $keep_tmp, $target_col_index,
                  $sort_buffer_size, $group_info_path_unsorted,
                  $read_final_paths_ref, $anno_SJ_ref, \%unanno_SJ, \%group_ends,
                  \@group_sort, \%sj_to_groups, \%summary_data);
  add_perfect_reads_to_group_info($out_dir, $chr, $keep_tmp, $has_tmp_sj,
                                  \%sj_to_groups, $group_info_path_unsorted);
  %sj_to_groups = ();
  sort_by_numeric_first_column($group_info_path_unsorted, $group_info_path,
                               $sort_buffer_size);
  unlink($group_info_path_unsorted) if !$keep_tmp;

  my %single_exon_isoforms_by_group = ();
  sort_single_exon_isoforms($SJ_dist, $single_exon_isoform_end_ref, \@group_sort,
                            \%group_ends, \%single_exon_isoforms_by_group);

  my @annotated_SJ_chr = ();
  my @unannotated_SJ_chr = ();
  sort_splice_junctions($anno_SJ_ref, \%unanno_SJ, \@annotated_SJ_chr,
                        \@unannotated_SJ_chr);

  my $min_anno_sj_i = 0;
  my $min_unanno_sj_i = 0;
  my %single_exon_isoform_total = ();
  my %chr_type = ();
  my %group_info_offsets = ();

  my %strand_symbol = ('0' => '+', '1' => '-', 'unknown' => 'unknown');
  for my $n (@group_sort) {
    my %SJ_read_improved = ();
    my %SJ_read_order = ();
    my %all_SJ_reads = ();
    my %read_info = ();
    my %read_filtered = ();
    my @group_reads = ();
    load_group_info($n, $group_info_path, \%group_info_offsets,
                    \%SJ_read_improved, \%SJ_read_order, \%all_SJ_reads,
                    \%read_info, \%read_filtered, \@group_reads);

    my %fsm = ();
    my %ism = ();
    my %nic = ();
    my %nnc = ();
    my %isoform_perfect_match = ();
    my @single_exon_reads = ();
    classify_reads_and_match_to_annotated_isoforms(
      $chr, $n, $SJ_dist, $should_create_default_handle, $anno_SJ_ref,
      $anno_SS_ref, $isoform_SJ_ref, $multi_exon_isoform_end_ref, \@group_reads,
      \%read_filtered, \%SJ_read_improved, \%read_info, \%chr_type, \%fsm, \%ism,
      \%nic, \%nnc, \%isoform_perfect_match, \@single_exon_reads, \%summary_data,
      $default_handle);
    %SJ_read_improved = ();
    %read_filtered = ();

    my %all_SJ_group = ();
    my ($group_start, $group_end) = @{$group_ends{$n}};
    if ($should_create_default_handle) {
      print $default_handle "group $n: ($group_start, $group_end)\n";
    }
    $min_anno_sj_i = get_sjs_for_group($min_anno_sj_i, $group_start, $group_end,
                                       1, \@annotated_SJ_chr, $anno_SJ_ref,
                                       \%all_SJ_group);
    $min_unanno_sj_i = get_sjs_for_group($min_unanno_sj_i, $group_start,
                                         $group_end, 0, \@unannotated_SJ_chr,
                                         \%unanno_SJ, \%all_SJ_group);

    my %annotated_isoforms_included = ();
    add_other_sjs_for_annotated_isoforms($anno_SJ_ref, $isoform_SJ_ref,
                                         \%annotated_isoforms_included,
                                         \%all_SJ_group);
    my @all_SJ_group_sort = ();
    my %all_SJ_group_ID = ();
    create_sorted_sj_mapping(\%all_SJ_group, \@all_SJ_group_sort,
                             \%all_SJ_group_ID);
    %all_SJ_group = ();

    my %fsm_chain = ();
    create_sj_id_chains($should_create_default_handle, $isoform_SJ_ref,
                        \%annotated_isoforms_included, \%all_SJ_group_ID,
                        \%fsm_chain, $default_handle);
    %annotated_isoforms_included = ();

    my %SJ_chain_read = ();
    my %SJ_chain_cat = ();
    create_sj_chain_mappings($n, \@group_reads, \%SJ_read_order,
                             \%all_SJ_group_ID, \%read_info, \%fsm, \%ism, \%nic,
                             \%nnc, \%SJ_chain_read, \%SJ_chain_cat,
                             \%summary_data);
    %SJ_read_order = ();
    @group_reads = ();
    %fsm = ();
    %ism = ();
    %nic = ();
    %nnc = ();
    %all_SJ_group_ID = ();

    my %SJ_chain_read_exon_end = ();
    my %SJ_included_chain = ();
    create_endpoint_mapping_by_chain($chr, $raw, $should_create_default_handle,
                                     \@all_SJ_group_sort, \%SJ_chain_cat,
                                     \%SJ_chain_read, \%SJ_chain_read_exon_end,
                                     \%SJ_included_chain, $default_handle);
    set_candidate_compatible_chains(\%SJ_included_chain, \%SJ_chain_cat,
                                    \%SJ_chain_read);
    %SJ_included_chain = ();

    my %unannotated_isoforms = ();
    my %longest_nc_valid = ();
    my @nc_valid_longest = ();
    my %substring_nc = ();
    if (defined $raw) {
      find_valid_novel_chains_raw(
        \%SJ_chain_cat, \%all_SJ_reads, \%SJ_chain_read, \%read_info,
        \@all_SJ_group_sort, $anno_SJ_complete_ref, $isoform_SJ_complete_ref,
        $read_num_cutoff, $read_ratio_cutoff, \%unannotated_isoforms,
        \%longest_nc_valid, $default_handle);
    } else {
      my @nnc_valid_all = ();
      if (exists $SJ_chain_cat{'nnc'}) {
        find_valid_novel_chains(
          $SJ_chain_cat{'nnc'}, \%all_SJ_reads, \%SJ_chain_read, \%read_info,
          \@all_SJ_group_sort, $anno_SJ_complete_ref, $isoform_SJ_complete_ref,
          $read_num_cutoff, $read_ratio_cutoff, \%unannotated_isoforms,
          \@nnc_valid_all, $default_handle);

        find_sub_chains($SJ_dist, $internal_boundary_limit,
                        $allow_longer_terminal_exons,
                        $should_create_default_handle, \@nnc_valid_all,
                        \%SJ_chain_read, \%SJ_chain_read_exon_end,
                        \%substring_nc, \%summary_data, $default_handle);

        add_longest_nnc_valid_chains(\@nnc_valid_all, \%substring_nc,
                                     \@nc_valid_longest, \%longest_nc_valid);
      }
      $summary_data{'num_validated_nnc_chains'} += scalar(@nnc_valid_all);
      @nnc_valid_all = ();

      my @nic_valid_all = ();
      if (exists $SJ_chain_cat{'nic'}) {
        find_valid_novel_chains(
          $SJ_chain_cat{'nic'}, \%all_SJ_reads, \%SJ_chain_read, \%read_info,
          \@all_SJ_group_sort, $anno_SJ_complete_ref, $isoform_SJ_complete_ref,
          $read_num_cutoff, $read_ratio_cutoff, \%unannotated_isoforms,
          \@nic_valid_all, $default_handle);

        find_sub_chains($SJ_dist, $internal_boundary_limit,
                        $allow_longer_terminal_exons,
                        $should_create_default_handle, \@nic_valid_all,
                        \%SJ_chain_read, \%SJ_chain_read_exon_end,
                        \%substring_nc, \%summary_data, $default_handle);

        add_longest_nic_valid_chains(
          $SJ_dist, $internal_boundary_limit, $allow_longer_terminal_exons,
          \@nic_valid_all, \%SJ_chain_read, \%SJ_chain_read_exon_end,
          \%substring_nc, \@nc_valid_longest, \%longest_nc_valid,
          \%summary_data);
      }
      $summary_data{'num_validated_nic_chains'} += scalar @nic_valid_all;
      @nic_valid_all = ();

      if (exists $SJ_chain_cat{'ism'}) {
        my @possible_nc_ism = ();
        my %chain_w_compatible_annotated_fsm = ();
        check_ism_chain_endpoints_and_valid_sjs(
          $SJ_dist, $internal_boundary_limit, $allow_longer_terminal_exons,
          $should_create_default_handle, \%fsm_chain, $anno_SJ_ref,
          $multi_exon_isoform_end_ref, \%all_SJ_reads, \%SJ_chain_read,
          \%read_info, \@all_SJ_group_sort, $anno_SJ_complete_ref,
          $isoform_SJ_complete_ref, $read_num_cutoff, $read_ratio_cutoff,
          $SJ_chain_cat{'ism'}, \%unannotated_isoforms, \@possible_nc_ism,
          \%chain_w_compatible_annotated_fsm, \%summary_data, $default_handle);

        find_sub_chains($SJ_dist, $internal_boundary_limit,
                        $allow_longer_terminal_exons,
                        $should_create_default_handle, \@possible_nc_ism,
                        \%SJ_chain_read, \%SJ_chain_read_exon_end,
                        \%substring_nc, \%summary_data, $default_handle);

        add_longest_ism_valid_chains(
          $SJ_dist, $internal_boundary_limit, $allow_longer_terminal_exons,
          \@possible_nc_ism, \%SJ_chain_read, \%SJ_chain_read_exon_end,
          \%chain_w_compatible_annotated_fsm, \%substring_nc, \@nc_valid_longest,
          \%longest_nc_valid, $SJ_chain_cat{'ism'}, \%summary_data);
      }
    }
    %all_SJ_reads = ();

    my %read_count_isoform = ();
    initialize_isoform_assignment_counts(
      \%SJ_chain_read, \%read_info, $samples_sort_ref, \%read_count_isoform);

    if (defined $raw) {
      output_chain_details_to_default(
        $chr, $should_create_default_handle, \%SJ_chain_read,
        \@all_SJ_group_sort, \%SJ_chain_cat, \%longest_nc_valid, \%substring_nc,
        \%unannotated_isoforms, $samples_sort_ref, \%read_count_isoform,
        $default_handle);
      next;
    }

    my %cat_SJ_chains_sort = ();
    sort_sj_chains_by_number_of_sjs(\%SJ_chain_cat, \%cat_SJ_chains_sort);

    assign_reads_to_compatible_nic_and_nnc_chains(
      $chr, $n, $should_create_tsv_compt_handle, $should_create_default_handle,
      \%cat_SJ_chains_sort, \%substring_nc, \%longest_nc_valid,
      \@nc_valid_longest, $samples_sort_ref, \%SJ_chain_read, \%read_info,
      \%SJ_chain_cat, \%read_count_isoform, \%summary_data, $tsv_compt_handle,
      $default_handle);

    # mt1 = 'more than one'
    my %ism_mt1_chain_count = ();
    assign_reads_to_compatible_ism_chains(
      $chr, $n, $group_start, $group_end, $SJ_dist, $internal_boundary_limit,
      $allow_longer_terminal_exons, $should_create_tsv_compt_handle,
      $should_create_default_handle, \%cat_SJ_chains_sort, \%substring_nc,
      \%longest_nc_valid, \@nc_valid_longest, \%SJ_chain_read,
      \%SJ_chain_read_exon_end, \%read_info, $multi_exon_isoform_end_ref,
      $anno_SJ_ref, \%fsm_chain, \@all_SJ_group_sort, \%ism_mt1_chain_count,
      \%SJ_chain_cat, \%read_count_isoform, \%summary_data, $tsv_compt_handle,
      $default_handle);
    %fsm_chain = ();

    my %possible_fsm_filtered = ();
    assign_filtered_reads_to_compatible_chains(
      $chr, $n, $SJ_dist, $internal_boundary_limit, $allow_longer_terminal_exons,
      $should_create_tsv_compt_handle, $should_create_default_handle,
      \%cat_SJ_chains_sort, \%SJ_chain_cat, \%longest_nc_valid, \%SJ_chain_read,
      \%SJ_chain_read_exon_end, \%read_info, $samples_sort_ref,
      \%possible_fsm_filtered, \%read_count_isoform, \%summary_data,
      $tsv_compt_handle, $default_handle);
    %SJ_chain_read_exon_end = ();

    initialize_em_from_direct_counts(\%read_count_isoform);

    output_chain_details_to_default(
      $chr, $should_create_default_handle, \%SJ_chain_read, \@all_SJ_group_sort,
      \%SJ_chain_cat, \%longest_nc_valid, \%substring_nc, \%unannotated_isoforms,
      $samples_sort_ref, \%read_count_isoform, $default_handle);
    %substring_nc = ();

    my @sorted_ism_mt1_chains = ();
    my %sorted_samples_by_ism_mt1_chain = ();
    my @sorted_possible_fsm_num = ();
    my %sorted_chains_by_fsm_num = ();
    sort_chains_for_em(
      \%ism_mt1_chain_count, \%possible_fsm_filtered, \%cat_SJ_chains_sort,
      \%SJ_chain_cat, \@sorted_ism_mt1_chains, \%sorted_samples_by_ism_mt1_chain,
      \@sorted_possible_fsm_num, \%sorted_chains_by_fsm_num);

    expectation_maximization(
      $chr, $n, $max_iterate, $should_create_default_handle, $samples_sort_ref,
      \%cat_SJ_chains_sort, \%SJ_chain_cat, \@nc_valid_longest,
      \%ism_mt1_chain_count, \@sorted_ism_mt1_chains,
      \%sorted_samples_by_ism_mt1_chain, \%possible_fsm_filtered,
      \@sorted_possible_fsm_num, \%sorted_chains_by_fsm_num,
      \%read_count_isoform, $default_handle);

    @sorted_ism_mt1_chains = ();
    %sorted_samples_by_ism_mt1_chain = ();
    @sorted_possible_fsm_num = ();
    %sorted_chains_by_fsm_num = ();
    @nc_valid_longest = ();
    %cat_SJ_chains_sort = ();
    %ism_mt1_chain_count = ();
    %possible_fsm_filtered = ();

    output_fsm_abundance(
      $chr, $n, $should_create_tsv_compt_handle, $should_create_default_handle,
      \%strand_symbol, $SJ_chain_cat{'fsm'}, \%SJ_chain_read,
      \@all_SJ_group_sort, \%read_count_isoform, \%read_info,
      \%isoform_perfect_match, $multi_exon_isoform_end_ref, $isoform_info_ref,
      $samples_sort_ref, \%summary_data, $abu_handle, $gtf_handle,
      $tsv_compt_handle, $default_handle);
    %isoform_perfect_match = ();

    output_novel_abundance(
      $chr, $n, $should_create_default_handle, \%strand_symbol,
      \%longest_nc_valid, \%SJ_chain_cat, \%SJ_chain_read, \@all_SJ_group_sort,
      \%read_count_isoform, \%unannotated_isoforms, $isoform_info_ref,
      $anno_SJ_ref, $anno_SJ_complete_ref, $samples_sort_ref, \%summary_data,
      $abu_handle, $gtf_handle, $default_handle);
    @all_SJ_group_sort = ();
    %SJ_chain_read = ();
    %SJ_chain_cat = ();
    %unannotated_isoforms = ();
    %longest_nc_valid = ();
    %read_count_isoform = ();

    my @single_exon_isoforms_group_sort = ();
    my @single_exon_reads_sort = ();
    sort_single_exon_isoforms_and_reads(
      $n, \%single_exon_isoforms_by_group, \@single_exon_reads,
      $single_exon_isoform_end_ref, \%read_info,
      \@single_exon_isoforms_group_sort, \@single_exon_reads_sort);

    my %single_exon_isoform_read_count = ();
    initialize_single_exon_assignment_counts(
      \@single_exon_isoforms_group_sort, $samples_sort_ref,
      \%single_exon_isoform_read_count);

    my %possible_single_exon_isoform = ();
    assign_reads_to_single_exon_isoforms(
      $chr, $n, $should_create_tsv_compt_handle, $should_create_default_handle,
      $SJ_dist, $internal_boundary_limit, $allow_longer_terminal_exons,
      \@single_exon_reads_sort, \@single_exon_isoforms_group_sort,
      $single_exon_isoform_end_ref, \@single_exon_reads, \%read_info,
      \%possible_single_exon_isoform, \%single_exon_isoform_read_count,
      \%summary_data, $tsv_compt_handle, $default_handle);
    @single_exon_reads = ();
    @single_exon_reads_sort = ();

    initialize_em_from_single_exon_counts(
      \@single_exon_isoforms_group_sort, $samples_sort_ref,
      \%single_exon_isoform_read_count);

    my @sorted_possible_single_exon_reads = ();
    sort_single_exon_isoforms_for_em(\%possible_single_exon_isoform,
                                     \@sorted_possible_single_exon_reads);

    single_exon_expectation_maximization(
      $chr, $n, $max_iterate, $should_create_default_handle,
      \@sorted_possible_single_exon_reads, \@single_exon_isoforms_group_sort,
      \%possible_single_exon_isoform, \%read_info, $samples_sort_ref,
      \%single_exon_isoform_read_count, $default_handle);
    @sorted_possible_single_exon_reads = ();
    %read_info = ();
    @single_exon_isoforms_group_sort = ();
    %possible_single_exon_isoform = ();

    update_single_exon_totals(
      $should_create_default_handle, \%single_exon_isoform_read_count,
      $single_exon_isoform_end_ref, $samples_sort_ref,
      \%single_exon_isoform_total, $default_handle);
  }
  unlink($group_info_path) if !$keep_tmp;
  %unanno_SJ = ();
  %group_ends = ();
  @group_sort = ();
  %single_exon_isoforms_by_group = ();
  @annotated_SJ_chr = ();
  @unannotated_SJ_chr = ();
  %group_info_offsets = ();

  output_single_exon_abundance(
    $chr, \%single_exon_isoform_total, $isoform_info_ref,
    $single_exon_isoform_end_ref, $samples_sort_ref, \%summary_data, $abu_handle,
    $gtf_handle);
  %single_exon_isoform_total = ();

  output_totals_to_default($should_create_default_handle, \%chr_type,
                           $default_handle);
  %chr_type = ();

  close $default_handle if $should_create_default_handle;
  close $tsv_compt_handle if $should_create_tsv_compt_handle;
  close $abu_handle;
  close $gtf_handle;

  my $summary_path = summary_path_for_chr($chr, $out_dir);
  store \%summary_data, $summary_path;
}

sub gtf_output {
  my ($gtf_handle, $chr, $isoform_ID, $gene_ID, $ends_ref, $SJs_ref, $strand,
      $type) = @_;
  my $start = $ends_ref->[0] + 1; # convert to 1-based
  my $end = $ends_ref->[1];
  my $attributes = "transcript_id \"$isoform_ID\"; gene_id \"$gene_ID\";";
  write_tsv_columns($gtf_handle, [$chr, $type, 'transcript', $start, $end, '.',
                                  $strand, '.', $attributes]);
  my @exons = ();
  if (scalar(@{$SJs_ref}) == 0) {
    @exons = ([$start, $end]);
  } else {
    @exons = ([$start]);
    for my $SJ (@{$SJs_ref}) {
      my $SJ_info = split_sj_string($SJ);
      push @{$exons[-1]}, $SJ_info->{'start'};
      push @exons, [$SJ_info->{'end'} + 1];
    }
    push @{$exons[-1]}, $end;
  }

  for my $exon_ID (0 .. $#exons) {
    my $exon_i = $exon_ID + 1;
    my $exon_start = $exons[$exon_ID][0];
    my $exon_end = $exons[$exon_ID][1];
    my $exon_attributes = "$attributes exon_number \"$exon_i\";";
    write_tsv_columns($gtf_handle, [$chr, $type, 'exon', $exon_start, $exon_end,
                                    '.', $strand, '.', $exon_attributes]);
  }
}

sub compare_SJ_chains_as_numeric_tuples {
  my ($chain_1, $chain_2) = @_;
  my @vals_1 = split('_', $chain_1);
  my @vals_2 = split('_', $chain_2);
  my $len_1 = scalar(@vals_1);
  my $len_2 = scalar(@vals_2);
  my $shorter_len = $len_1 < $len_2 ? $len_1 : $len_2;
  for my $i (0 .. ($shorter_len - 1)) {
    my $cmp_result = $vals_1[$i] cmp $vals_2[$i];
    if ($cmp_result == 0) {
      next;
    }
    my $val_1_is_numeric = $vals_1[$i] =~ /^[0-9]+$/;
    my $val_2_is_numeric = $vals_2[$i] =~ /^[0-9]+$/;
    if ($val_1_is_numeric) {
      if ($val_2_is_numeric) {
        return $vals_1[$i] < $vals_2[$i] ? -1 : 1;
      }
      return -1;
    }
    if ($val_2_is_numeric) {
      return 1;
    }
    return $vals_1[$i] lt $vals_2[$i] ? -1 : 1;
  }
  if ($len_1 == $len_2) {
    return 0;
  }
  return $len_1 < $len_2 ? -1 : 1;
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

sub summary_path_for_chr {
  my ($chr, $out_dir) = @_;
  return $out_dir."/$chr.espresso_q_summary.tmp";
}

# package must return true
1;

use strict;
use warnings;

use threads;
use threads::shared;
use Thread::Queue;
use Getopt::Long;
use Storable qw(freeze thaw);

use File::Basename qw(dirname);
use lib dirname(__FILE__);
use ESPRESSO_Version;

sub parse_args {
  my %args = ();
  $args{'arguments_before_parsing'} = "@ARGV";

  $args{'SJFS_dist'} = undef;
  $args{'SJFS_dist_add'} = undef;
  $args{'SJ_bed'} = undef;
  $args{'anno'} = undef;
  $args{'chrM'} = undef;
  $args{'cont_del_max'} = undef;
  $args{'fa'} = undef;
  $args{'group_extra_range'} = undef;
  $args{'alignment_read_groups'} = undef;
  $args{'help'} = undef;
  $args{'inserted_cont_cutoff'} = undef;
  $args{'keep_tmp'} = undef;
  $args{'list_samples'} = undef;
  $args{'mapq_cutoff'} = undef;
  $args{'num_thread'} = undef;
  $args{'out'} = undef;
  $args{'read_num_cutoff'} = undef;
  $args{'read_ratio_cutoff'} = undef;
  $args{'sort_buffer_size'} = undef;

  Getopt::Long::GetOptions(
    'SJFS_dist=i' => \$args{'SJFS_dist'},
    'SJFS_dist_add=i' => \$args{'SJFS_dist_add'},
    'SJ_bed|B=s' => \$args{'SJ_bed'},
    'anno|A=s' => \$args{'anno'},
    'chrM|M=s' => \$args{'chrM'},
    'cont_del_max|C=i' => \$args{'cont_del_max'},
    'fa|F=s' => \$args{'fa'},
    'group_extra_range|E=i' => \$args{'group_extra_range'},
    'alignment_read_groups' => \$args{'alignment_read_groups'},
    'help|H!' => \$args{'help'},
    'inserted_cont_cutoff=i' => \$args{'inserted_cont_cutoff'},
    'keep_tmp|K!' => \$args{'keep_tmp'},
    'list_samples|L=s' => \$args{'list_samples'},
    'mapq_cutoff|Q=i' => \$args{'mapq_cutoff'},
    'num_thread|T=i' => \$args{'num_thread'},
    'out|O=s' => \$args{'out'},
    'read_num_cutoff|N=i' => \$args{'read_num_cutoff'},
    'read_ratio_cutoff|R=f' => \$args{'read_ratio_cutoff'},
    'sort_buffer_size=s' => \$args{'sort_buffer_size'}
  );

  if (!defined $args{'out'}) {
    $args{'out'} = '.';
  }
  if (!defined $args{'read_num_cutoff'}) {
    $args{'read_num_cutoff'} = 2;
  }
  if (!defined $args{'read_ratio_cutoff'}) {
    $args{'read_ratio_cutoff'} = 0;
  }
  if (!defined $args{'group_extra_range'}) {
    $args{'group_extra_range'} = 0;
  }
  if (!defined $args{'mapq_cutoff'}) {
    $args{'mapq_cutoff'} = 1;
  }
  if (!defined $args{'inserted_cont_cutoff'}) {
    $args{'inserted_cont_cutoff'} = 20;
  }
  if (!defined $args{'SJFS_dist'}) {
    $args{'SJFS_dist'} = 10;
  }
  if (!defined $args{'SJFS_dist_add'}) {
    $args{'SJFS_dist_add'} = 15;
  }
  if (!defined $args{'num_thread'}) {
    $args{'num_thread'} = 5;
  }
  if (!defined $args{'chrM'}) {
    $args{'chrM'} = 'chrM';
  }
  if (!defined $args{'cont_del_max'}) {
    $args{'cont_del_max'} = 50;
  }
  if (defined $args{'sort_buffer_size'}) {
    $args{'sort_buffer_size'} = "--buffer-size=$args{'sort_buffer_size'}";
  } else {
    $args{'sort_buffer_size'} = '--buffer-size=2G';
  }

  return \%args;
}

sub show_help_message {
  my $version_number = ESPRESSO_Version::get_version_number();
  my $version_string = "S_$version_number";

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
}

sub starts_with {
  my ($string, $value) = @_;

  return index($string, $value) == 0;
}

sub check_bit_from_base_10 {
  my ($num, $bit_i) = @_;

  my $binary_string = sprintf('%b', $num);
  my $len_string = length($binary_string);
  if ($bit_i >= $len_string) {
    return 0;
  }
  my $string_i = $len_string - ($bit_i + 1);
  return substr($binary_string, $string_i, 1);
}

sub split_sj_string {
  my ($sep_SJ, $sj_string) = @_;

  my @parts = split($sep_SJ, $sj_string);
  my %result;
  $result{'chr'} = $parts[0];
  $result{'start'} = $parts[1];
  $result{'end'} = $parts[2];
  return \%result
}

sub write_tsv_columns {
  my ($handle, $columns_ref) = @_;

  my $with_tabs = join("\t", @{$columns_ref});
  print $handle "$with_tabs\n";
}

sub write_group_and_sj_list_lines {
  my ($group_ID, $group_reads_ref, $group_info_ref, $SJ_read_ref, $group_handle,
      $sj_list_handle) = @_;

  my $comma_reads = join(',', @{$group_reads_ref});
  $comma_reads .= ',';  # TODO trailing comma
  write_tsv_columns(
    $group_handle,
    [$group_ID, $group_info_ref->{'chr'}, $group_info_ref->{'start'},
     $group_info_ref->{'end'}, $comma_reads]);
  if (scalar(keys %{$SJ_read_ref}) > 0) {
    while (my ($SJ, $SJ_info_ref) = each %{$SJ_read_ref}) {
      my $num_reads = 0;
      my @reads;
      my @perfect_reads;
      while (my ($read_id, $is_perfect)
             = each %{$SJ_info_ref->{'perfect_by_read'}}) {
        $num_reads ++;
        push @reads, $read_id;
        if ($is_perfect == 1) {
          push @perfect_reads, $read_id;
        }
      }
      my $comma_perfect = 'NA';
      my $num_perfect = scalar(@perfect_reads);
      if ($num_perfect > 0) {
        $comma_perfect = join(',', @perfect_reads);
        $comma_perfect .= ',';  # TODO trailing comma
      }
      my $comma_reads = join(',', @reads);
      $comma_reads .= ',';  # TODO trailing comma
      write_tsv_columns(
        $sj_list_handle,
        [$group_ID, $SJ, $SJ_info_ref->{'chr'}, $SJ_info_ref->{'start'},
         $SJ_info_ref->{'end'}, $num_perfect, $num_reads, $comma_perfect,
         $comma_reads]);
    }
  }
}

sub min {
  my ($x, $y) = @_;

  if ($x < $y) {
    return $x;
  }
  return $y;
}

sub max {
  my ($x, $y) = @_;

  if ($x > $y) {
    return $x;
  }
  return $y;
}

# Assumes that length($x) == length($y)
sub hamming_distance {
  my ($x, $y) = @_;

  # Any matching chars will XOR to null (\0).
  my $xor_bytes = $x ^ $y;
  # tr counts the nulls
  my $num_matches = $xor_bytes =~ tr[\0][\0];
  return length($x) - $num_matches;
}

sub comp_rev {
  my ($orig) = @_;

  my $seq = reverse($orig);
  $seq =~ tr/ATCG/TAGC/;
  return $seq;
}

sub process_cigar_string {
  my ($cigar, $continuous_deletion_max) = @_;

  my $read_length = 0;
  my $mapped_nt_genome = 0;
  my $mapped_nt_genomeN = 0;

  # TODO handle hard clipping differently
  $cigar =~ s/H/S/g;
  my @counts = split(/[MSIDN=X]/, $cigar);
  my @operators = split /\d+/, $cigar;
  shift @operators;

  my @exonIntron = (0);
  my @SJ_dist_read = (0);
  my @ID_from_current_SJ = ([]);
  my @M_from_current_SJ = ([]);
  my $dist_from_prev_SJ = 0;
  my $dist_from_read_start = 0;
  my $max_I = 0;
  my @clip_ends = (0, 0);
  for my $i (0 .. $#operators) {
    my $count = $counts[$i];
    my $op = $operators[$i];
    if (($op eq 'N') or (($op eq 'D') and ($count > $continuous_deletion_max))) {
      push @exonIntron, $count;
      push @exonIntron, 0;
      push @SJ_dist_read, $SJ_dist_read[-1];
      $dist_from_prev_SJ = 0;
      push @ID_from_current_SJ, [];
      push @M_from_current_SJ, [];
      $mapped_nt_genomeN += $count;
    } elsif (($op eq 'M') or ($op eq '=') or ($op eq 'X')) {
      $exonIntron[-1] += $count;
      $SJ_dist_read[-1] += $count;
      $read_length += $count;
      $mapped_nt_genome += $count;
      $mapped_nt_genomeN += $count;
      for my $id_op_ref (@{$ID_from_current_SJ[-1]}) {
        $id_op_ref->{'ref_nt_to_next_sj'} += $count;
      }
      my $m_ops_ref = $M_from_current_SJ[-1];
      for my $m_op_ref (@{$m_ops_ref}) {
        $m_op_ref->{'ref_nt_to_next_sj'} += $count;
      }

      my %details;
      $details{'count'} = $count;
      $details{'ref_nt_to_next_sj'} = 0;
      $details{'ref_nt_from_sj'} = $dist_from_prev_SJ;
      $details{'dist_from_read_start'} = $dist_from_read_start;
      push @{$m_ops_ref}, \%details;
      $dist_from_prev_SJ += $count;
      $dist_from_read_start += $count;
    } elsif ($op eq 'D') {
      $exonIntron[-1] += $count;
      $mapped_nt_genome += $count;
      $mapped_nt_genomeN += $count;
      for my $m_op_ref (@{$M_from_current_SJ[-1]}) {
        $m_op_ref->{'ref_nt_to_next_sj'} += $count;
      }
      my $id_ops_ref = $ID_from_current_SJ[-1];
      for my $id_op_ref (@{$id_ops_ref}) {
        $id_op_ref->{'ref_nt_to_next_sj'} += $count;
      }

      my %details;
      $details{'count'} = 0 - $count;
      $details{'ref_nt_to_next_sj'} = 0;
      $details{'ref_nt_from_sj'} = $dist_from_prev_SJ;
      push @{$id_ops_ref}, \%details;
      $dist_from_prev_SJ += $count;
    } elsif ($op eq 'I') {
      $SJ_dist_read[-1] += $count;
      $read_length += $count;
      if ($count > $max_I) {
        $max_I = $count;
      }

      my %details;
      $details{'count'} = $count;
      $details{'ref_nt_to_next_sj'} = 0;
      $details{'ref_nt_from_sj'} = $dist_from_prev_SJ;
      push @{$ID_from_current_SJ[-1]}, \%details;
      $dist_from_read_start += $count;
    } elsif ($op eq 'S') {
      $SJ_dist_read[-1] += $count;
      $read_length += $count;
      $dist_from_read_start += $count;
      if ($mapped_nt_genome > 0) {
        $clip_ends[1] += $count;
      } else {
        $clip_ends[0] += $count;
      }
    } else {
      die "unexpected cigar operator $op: $cigar";
    }
  }
  pop @SJ_dist_read;

  my %result;
  $result{'exonIntronRef'} = \@exonIntron;
  $result{'read_length'} = $read_length;
  $result{'mappedGenome'} = $mapped_nt_genome;
  $result{'mappedGenomeN'} = $mapped_nt_genomeN;
  $result{'ID_from_current_SJ'} = \@ID_from_current_SJ;
  $result{'M_from_current_SJ'} = \@M_from_current_SJ;
  $result{'SJ_dist_read'} = \@SJ_dist_read;
  $result{'clip_ends'} = \@clip_ends;
  $result{'max_I'} = $max_I;
  return \%result;
}

sub maybe_start_new_alignment_group {
  my ($chr, $start, $end, $file, $group_extra_range, $chr_seq_len_ref, $group_ID,
      $group_info_ref, $processed_chrs_ref, $group_reads_ref, $SJ_read_ref,
      $group_handle, $sj_list_handle) = @_;

  my $is_chr_switch = 0;
  my $is_switch_to_old_chr = 0;
  my $is_lower_start_coord = 0;
  if (exists $group_info_ref->{'chr'}) {
    $is_chr_switch = $chr ne $group_info_ref->{'chr'};
    $is_switch_to_old_chr = ($is_chr_switch
                             and exists $processed_chrs_ref->{$chr});
    $is_lower_start_coord = (!$is_chr_switch
                             and ($start < $group_info_ref->{'start'}));
  }

  if ($is_switch_to_old_chr or $is_lower_start_coord) {
    die "$file is not sorted";
  }
  $processed_chrs_ref->{$chr} = 1;

  if (!exists $group_info_ref->{'chr'} or $is_chr_switch
      or ($start > ($group_info_ref->{'end'} + $group_extra_range))) {
    if (exists $group_info_ref->{'chr'}) {
      write_group_and_sj_list_lines(
        $group_ID, $group_reads_ref, $group_info_ref, $SJ_read_ref,
        $group_handle, $sj_list_handle);
      $group_ID ++;
    }
    %{$SJ_read_ref} = ();
    @{$group_reads_ref} = ();

    $group_info_ref->{'chr'} = $chr;
    $group_info_ref->{'start'} = max($start - $group_extra_range, 0);
    $group_info_ref->{'end'} = min($end + $group_extra_range,
                                   $chr_seq_len_ref->{$chr});
  } else {
    my $new_end = min($end + $group_extra_range, $chr_seq_len_ref->{$chr});
    if ($group_info_ref->{'end'} < $new_end) {
      $group_info_ref->{'end'} = $new_end;
    }
  }

  return $group_ID;
}

sub check_indels_and_substitutions_by_sj {
  my ($read_id, $chr, $start, $seq_len, $sequence, $SJFS_dist, $sep_SJ,
      $chr_seq_ref, $align_info, $IDS_SJ_ref, $SJ_cor_seq_ref,
      $SJ_read_ref) = @_;

  my $prev_exon_start = $start - 1;
  my $insert_num_bg = 0;
  my $delete_num_bg = 0;
  my $substi_num_bg = 0;
  my $total_bg = 0;
  my $num_sjs = $#{$align_info->{'ID_from_current_SJ'}};
  # Loop over the SJs. When $i == $num_sjs, just check the final exon.
  # TODO handle the final exon in the ($num_sjs - 1) iteration.
  # TODO only process each ID_from_current_SJ and M_from_current_SJ once
  #      (instead of as a prev_exon and as a next_exon)
  for my $i (0 .. $num_sjs) {
    my $prev_exon_i = $i * 2;
    my $intron_i = $prev_exon_i + 1;
    my $prev_exon_end = ($prev_exon_start
                         + $align_info->{'exonIntronRef'}[$prev_exon_i]);

    # check region not within SJFS_dist of SJ (background)
    my $prev_exon_ID_ref = $align_info->{'ID_from_current_SJ'}[$i];
    my $prev_exon_M_ref = $align_info->{'M_from_current_SJ'}[$i];
    # TODO should $i == 0 check ref_nt_from_sj since this is before the 1st SJ?
    # TODO should $i == $num_sjs check ref_nt_to_next_sj?
    # TODO should deletion consider overlap with the SJFS_dist cutoff?
    for my $ID_info (@{$prev_exon_ID_ref}) {
      if (($ID_info->{'ref_nt_to_next_sj'} > $SJFS_dist)
          and ($ID_info->{'ref_nt_from_sj'} > $SJFS_dist)) {
        if ($ID_info->{'count'} > 0) {
          $insert_num_bg += $ID_info->{'count'};
        } else {
          $delete_num_bg -= $ID_info->{'count'};
        }
        $total_bg += abs($ID_info->{'count'});
      }
    }

    # TODO should $i == 0 check ref_nt_from_sj?
    # TODO should $i == $num_sjs check ref_nt_to_next_sj?
    for my $M_info (@{$prev_exon_M_ref}) {
      my $count = $M_info->{'count'};
      if (($M_info->{'ref_nt_to_next_sj'} + $count <= $SJFS_dist)
          or ($M_info->{'ref_nt_from_sj'} + $count <= $SJFS_dist)) {
        next;
      }

      my $substr_length = $M_info->{'count'};
      if ($M_info->{'ref_nt_to_next_sj'} < $SJFS_dist) {
        my $overlap = $SJFS_dist - $M_info->{'ref_nt_to_next_sj'};
        $substr_length -= $overlap;
      }

      my $ref_start = $prev_exon_start + $M_info->{'ref_nt_from_sj'};
      my $read_start = $M_info->{'dist_from_read_start'};
      if ($M_info->{'ref_nt_from_sj'} < $SJFS_dist) {
        my $overlap = $SJFS_dist - $M_info->{'ref_nt_from_sj'};
        $substr_length -= $overlap;
        $ref_start += $overlap;
        $read_start += $overlap;
      }

      if ($substr_length <= 0) {
        next;
      }

      my $seq_ref = substr($chr_seq_ref->{$chr}, $ref_start, $substr_length);
      if (($read_start + $substr_length) > $seq_len) {
        next; # TODO due to supplementary alignment? Or treating H like S?
      }
      my $seq_read = substr($sequence, $read_start, $substr_length);
      $substi_num_bg += hamming_distance("\U$seq_read", "\U$seq_ref");
      $total_bg += $substr_length;
    }

    if ($i == $num_sjs) {
      last;
    }

    my $next_exon_start = ($prev_exon_end
                           + $align_info->{'exonIntronRef'}[$intron_i]);
    my $SJ = $chr.$sep_SJ.$prev_exon_end.$sep_SJ.$next_exon_start;

    my $next_exon_ID_ref = $align_info->{'ID_from_current_SJ'}[$i + 1];
    my $next_exon_M_ref = $align_info->{'M_from_current_SJ'}[$i + 1];

    my $insert_num = 0;
    my $deletion_num = 0;
    my $subst_num = 0;
    # check region within SJFS_dist before SJ
    for my $ID_info (@{$prev_exon_ID_ref}) {
      if ($ID_info->{'ref_nt_to_next_sj'} < $SJFS_dist) {
        if ($ID_info->{'count'} > 0) {
          $insert_num += $ID_info->{'count'};
        } else {
          my $del_count = -$ID_info->{'count'};
          my $possible_overlap = $SJFS_dist - $ID_info->{'ref_nt_to_next_sj'};
          $deletion_num += min($del_count, $possible_overlap);
        }
      }
    }

    for my $M_info (@{$prev_exon_M_ref}) {
      if ($M_info->{'ref_nt_to_next_sj'} >= $SJFS_dist) {
        next;
      }

      my $possible_overlap = $SJFS_dist - $M_info->{'ref_nt_to_next_sj'};
      my $substr_length;
      my $ref_start = $prev_exon_end - $M_info->{'ref_nt_to_next_sj'};
      my $read_start = $M_info->{'dist_from_read_start'};
      if ($M_info->{'count'} < $possible_overlap) {
        $substr_length = $M_info->{'count'};
        $ref_start -= $M_info->{'count'};
      } else {
        $substr_length = $possible_overlap;
        $ref_start -= $possible_overlap;
        $read_start += $M_info->{'count'} - $possible_overlap;
      }

      my $seq_ref = substr($chr_seq_ref->{$chr}, $ref_start, $substr_length);
      my $required_len = $read_start + $substr_length;
      if ($seq_len < $required_len) {
        next; # TODO due to supplementary alignment?
      }
      my $seq_read = substr($sequence, $read_start, $substr_length);
      $subst_num += hamming_distance("\U$seq_read", "\U$seq_ref");
    }

    # check region within SJFS_dist after SJ
    for my $ID_info (@{$next_exon_ID_ref}) {
      if ($ID_info->{'ref_nt_from_sj'} < $SJFS_dist) {
        if ($ID_info->{'count'} > 0) {
          $insert_num += $ID_info->{'count'};
        } else {
          my $del_count = -$ID_info->{'count'};
          my $possible_overlap = $SJFS_dist - $ID_info->{'ref_nt_from_sj'};
          $deletion_num += min($del_count, $possible_overlap);
        }
      }
    }

    for my $M_info (@{$next_exon_M_ref}) {
      if ($M_info->{'ref_nt_from_sj'} >= $SJFS_dist) {
        next;
      }

      my $possible_overlap = $SJFS_dist - $M_info->{'ref_nt_from_sj'};
      my $ref_start = $next_exon_start + $M_info->{'ref_nt_from_sj'};
      my $read_start = $M_info->{'dist_from_read_start'};
      my $substr_length = min($M_info->{'count'}, $possible_overlap);

      my $seq_ref = substr($chr_seq_ref->{$chr}, $ref_start, $substr_length);
      my $required_len = $read_start + $substr_length;
      if ($seq_len < $required_len) {
        next; # TODO due to supplementary alignment?
      }
      my $seq_read = substr($sequence, $read_start, $substr_length);
      $subst_num += hamming_distance("\U$seq_read", "\U$seq_ref");
    }

    my %sj_cor_seq_details;
    $sj_cor_seq_details{'read_offset'} = $align_info->{'SJ_dist_read'}[$i];
    $sj_cor_seq_details{'sj'} = $SJ;
    push @{$SJ_cor_seq_ref}, \%sj_cor_seq_details;
    my %ids_details;
    $ids_details{'insertion'} = $insert_num;
    $ids_details{'deletion'} = $deletion_num;
    $ids_details{'substitution'} = $subst_num;
    push @{$IDS_SJ_ref}, \%ids_details;

    my $isPerfect = 0;
    if (($insert_num + $deletion_num + $subst_num) == 0) {
      $isPerfect = 1;
    }
    if (!exists $SJ_read_ref->{$SJ}) {
      my $sj_info = split_sj_string($sep_SJ, $SJ);
      my %details;
      $details{'chr'} = $sj_info->{'chr'};
      $details{'start'} = $sj_info->{'start'};
      $details{'end'} = $sj_info->{'end'};
      my %perfect_by_read;
      $perfect_by_read{$read_id} = $isPerfect;
      $details{'perfect_by_read'} = \%perfect_by_read;
      $SJ_read_ref->{$SJ} = \%details;
    } elsif ($isPerfect == 1
             or !exists $SJ_read_ref->{$SJ}{'perfect_by_read'}{$read_id}) {
      $SJ_read_ref->{$SJ}{'perfect_by_read'}{$read_id} = $isPerfect;
    }

    $prev_exon_start = $next_exon_start;
  }

  my %result;
  $result{'insert_num'} = $insert_num_bg;
  $result{'deletion_num'} = $delete_num_bg;
  $result{'substitution_num'} = $substi_num_bg;
  $result{'total_num'} = $total_bg;
  return \%result;
}

sub write_sam_list_line {
  my ($group_ID, $line_num, $read_id, $sample, $flag, $chr, $start, $seq_len,
      $sequence, $mapq, $end, $notSameStrand, $attribute_string, $align_info,
      $IDS_SJ_ref, $background_counts, $SJ_cor_seq_ref, $sam_list_handle) = @_;

  my @out_columns;
  push @out_columns, $group_ID;
  push @out_columns, $line_num;
  push @out_columns, $read_id;
  push @out_columns, $sample;
  push @out_columns, $align_info->{'read_length'};
  push @out_columns, $flag;
  push @out_columns, $chr;
  push @out_columns, $start;
  push @out_columns, $mapq;
  push @out_columns, $end;
  push @out_columns, $notSameStrand;
  push @out_columns, $align_info->{'mappedGenome'};
  my $comma_ends = join(',', @{$align_info->{'clip_ends'}});
  $comma_ends .= ',';  # TODO
  push @out_columns, $comma_ends;
  my $comma_exons = join(',', @{$align_info->{'exonIntronRef'}});
  $comma_exons .= ',';  # TODO
  push @out_columns, $comma_exons;

  if (@{$IDS_SJ_ref} > 0) {
    my @sj_strings;
    for my $details (@{$IDS_SJ_ref}) {
      my $string = join(';', ($details->{'insertion'}, $details->{'deletion'},
                              $details->{'substitution'}));
      push @sj_strings, $string;
    }
    my $comma_sjs = join(',', @sj_strings);
    $comma_sjs .= ',';  # TODO
    push @out_columns, $comma_sjs;
  } else {
    push @out_columns, 'NA';
  }

  my $NM_num = '';
  if (defined $attribute_string) {
    $attribute_string =~ /NM:i:(\d+)/;
    $NM_num = $1;
  }
  push @out_columns, $NM_num;
  push @out_columns, $background_counts->{'insert_num'};
  push @out_columns, $background_counts->{'deletion_num'};
  push @out_columns, $background_counts->{'substitution_num'};
  push @out_columns, $background_counts->{'total_num'};
  if (@{$SJ_cor_seq_ref} > 0) {
    my @sj_strings;
    for my $details (@{$SJ_cor_seq_ref}) {
      my $string = join(';', ($details->{'read_offset'}, $details->{'sj'},
                              '0'));
      push @sj_strings, $string;
    }
    my $comma_sjs = join(',', @sj_strings);
    $comma_sjs .= ',';  # TODO
    push @out_columns, $comma_sjs;
  } else {
    push @out_columns, 'NA';
  }

  if ($align_info->{'read_length'} == $seq_len) {
    if ($notSameStrand == 1) {
      push @out_columns, comp_rev($sequence);
    } else {
      push @out_columns, $sequence;
    }
  } elsif ($align_info->{'read_length'} > $seq_len) {
    push @out_columns, 'NA';
  } else {
    my $message = 'read length from CIGAR was shorter than read sequence:';
    $message .= " $align_info->{'read_length'}, $seq_len";
    die $message;
  }

  write_tsv_columns($sam_list_handle, \@out_columns);
}

sub parallel_scan1 {
  my ($file, $is_sam, $input_ref, $summary_data_ref) = @_;

  my $file_ID_ref = $input_ref->{'file_ID'};
  my $out_dir = $input_ref->{'out'};
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

  my $num = $file_ID_ref->{$file};
  my $sample = $sam_sample_ref->{$file};

  my $sam_handle;
  if ($is_sam) {
    open($sam_handle, '<', $file) or die "cannot open $file: $!";
  } else {
    open($sam_handle, "samtools view -h $file |")
      or die "cannot open $file using samtools: $!";
  }

  my $num_dir = "$out_dir/$num";
  if (!-d $num_dir) {
    mkdir($num_dir) or die "cannot mkdir $num_dir: $!";
  }

  my $sam_list_path = "$num_dir/sam.list";
  my $group_list_path = "$num_dir/group.list";
  my $sj_list_path = "$num_dir/sj.list";
  open(my $sam_list_handle, '>', $sam_list_path)
    or die "cannot write tmp $sam_list_path: $!";
  open(my $group_handle, '>', $group_list_path)
    or die "cannot write tmp $group_list_path: $!";
  open(my $sj_list_handle, '>', $sj_list_path)
    or die "cannot write tmp $sj_list_path: $!";

  my %SJ_read;
  my %group_info;
  my @group_reads;
  my %processed_chrs;
  my $group_ID = 1;
  my $line_num = 0;
  my $past_headers = 0;
  while (<$sam_handle>) {
    chomp;
    $line_num++;
    my @line = split /\t/;
    if (!$past_headers and starts_with($line[0], '@')) {
      next;
    }
    $past_headers = 1;

    my $read_id = $line[0];
    my $flag = $line[1];
    my $chr = $line[2];
    my $start = $line[3];  # 1-based
    my $mapq = $line[4];
    my $cigar = $line[5];
    my $sequence = $line[9];
    my $attribute_string = $line[11];
    my $seq_len = length($sequence);
    my $notSameStrand = check_bit_from_base_10($flag, 4);
    my $is_secondary = check_bit_from_base_10($flag, 8);
    if ($chr eq $chrM) {
      $summary_data_ref->{'num_alignments_filtered_for_chrM'} ++;
      next;
    }
    if ($is_secondary) {
      $summary_data_ref->{'num_alignments_filtered_for_secondary'} ++;
      next;
    }
    if ($mapq < $mapq_cutoff) {
      $summary_data_ref->{'num_alignments_filtered_for_mapping_quality'} ++;
      next;
    }
    if (!exists $chr_seq_len_ref->{$chr}) {
      die "chr name ($chr) not recognized for $read_id in $file";
    }

    my $align_info = process_cigar_string($cigar, $cont_del_max);
    if ($align_info->{'max_I'} >= $inserted_cont_cutoff) {
      $summary_data_ref->{'num_alignments_filtered_for_max_insertion'} ++;
      next;
    }
    my $end = ($start - 1) + $align_info->{'mappedGenomeN'};
    if ($end > $chr_seq_len_ref->{$chr}) {
      $summary_data_ref->{'num_alignments_filtered_for_past_chr_end'} ++;
      next;
    }

    $group_ID = maybe_start_new_alignment_group(
      $chr, $start, $end, $file, $group_extra_range, $chr_seq_len_ref, $group_ID,
      \%group_info, \%processed_chrs, \@group_reads, \%SJ_read, $group_handle,
      $sj_list_handle);

    push @group_reads, $read_id;

    my @IDS_SJ;
    my @SJ_cor_seq;
    my $background_counts = check_indels_and_substitutions_by_sj(
      $read_id, $chr, $start, $seq_len, $sequence, $SJFS_dist, $sep_SJ,
      $chr_seq_ref, $align_info, \@IDS_SJ, \@SJ_cor_seq, \%SJ_read);

    write_sam_list_line($group_ID, $line_num, $read_id, $sample, $flag, $chr,
                        $start, $seq_len, $sequence, $mapq, $end, $notSameStrand,
                        $attribute_string, $align_info, \@IDS_SJ,
                        $background_counts, \@SJ_cor_seq, $sam_list_handle);
  }

  close $sam_handle;
  close $sam_list_handle;

  if (exists $group_info{'chr'}) {
    write_group_and_sj_list_lines($group_ID, \@group_reads, \%group_info,
                                  \%SJ_read, $group_handle, $sj_list_handle);
  }

  close $group_handle;
  close $sj_list_handle;
}

sub select_best_alignment_and_get_sequence {
  my ($read_id, $right_length_index, $read_alignments_ref, $read_info_ref) = @_;

  my $is_short = 0;
  my @sort_alignment = sort {
    $b->{'mappedGenome'} <=> $a->{'mappedGenome'}
      or $b->{'mapq'} <=> $a->{'mapq'};
  } @{$read_alignments_ref};
  my $best = $sort_alignment[0];

  my %details;
  $details{'read_length'} = $best->{'read_length'};
  $details{'line_num'} = $best->{'line_num'};
  $details{'mappedGenome'} = $best->{'mappedGenome'};
  if (scalar(@sort_alignment) > 1) {
    $read_info_ref->{$read_id} = \%details;
  }

  if (length($best->{'readSeq'}) != $best->{'read_length'}) {
    my $seq;
    if ($right_length_index >= 0) {
      $seq = $read_alignments_ref->[$right_length_index]{'readSeq'};
    } else {
      $seq = 'short';
      $is_short = 1;
    }
    $details{'seq'} = $seq;
    $read_info_ref->{$read_id} = \%details;
  }

  return $is_short;
}

sub load_info_for_reads_with_multiple_alignments_or_missing_sequence {
  my ($sam_list2_path, $output_titles_ID_ref, $read_info_ref) = @_;

  my $any_short = 0;
  my $read_id;
  my @read_alignments;
  my $right_length_index = -1;
  open(my $sam_list2_handle, '<', $sam_list2_path)
    or die "cannot read tmp $sam_list2_path: $!";
  while (<$sam_list2_handle>) {
    chomp;
    my @line = split /\t/;
    my $line_read_id = $line[$output_titles_ID_ref->{'readID'}];
    if ($line_read_id eq 'readID') {
      next;
    }

    if (defined $read_id and $read_id ne $line_read_id) {
      my $is_short = select_best_alignment_and_get_sequence(
        $read_id, $right_length_index, \@read_alignments, $read_info_ref);
      if ($is_short) {
        $any_short = 1;
      }
      @read_alignments = ();
      $right_length_index = -1;
    }

    $read_id = $line_read_id;
    my %details;
    $details{'read_length'} = $line[$output_titles_ID_ref->{'read_length'}];
    $details{'mappedGenome'} = $line[$output_titles_ID_ref->{'mappedGenome'}];
    $details{'mapq'} = $line[$output_titles_ID_ref->{'mapq'}];
    $details{'line_num'} = $line[$output_titles_ID_ref->{'line_num'}];
    $details{'readSeq'} = $line[$output_titles_ID_ref->{'readSeq'}];
    push @read_alignments, \%details;
    if ((length($details{'readSeq'}) == $details{'read_length'})
        and ($details{'readSeq'} ne 'NA')) {
      $right_length_index = $#read_alignments;
    }
  }
  close $sam_list2_handle;

  if (defined $read_id) {
    my $is_short = select_best_alignment_and_get_sequence(
      $read_id, $right_length_index, \@read_alignments, $read_info_ref);
    if ($is_short) {
      $any_short = 1;
    }
    @read_alignments = ();
    $right_length_index = -1;
  }

  return $any_short;
}

sub check_sam_for_full_sequence {
  my ($is_sam, $file, $read_info_ref) = @_;

  my $sam_handle;
  if ($is_sam) {
    open($sam_handle, '<', $file) or die "cannot open $file: $!";
  } else {
    open($sam_handle, "samtools view -h $file |")
      or die "cannot open $file using samtools: $!";
  }
  while (<$sam_handle>) {
    chomp;
    my @line = split /\t/;
    my $read_id = $line[0];
    my $flag = $line[1];
    my $seq = $line[9];
    if (exists $read_info_ref->{$read_id}
        and exists $read_info_ref->{$read_id}{'seq'}
        and ($read_info_ref->{$read_id}{'seq'} eq 'short')
        and (length($seq) == $read_info_ref->{$read_id}{'read_length'})) {
      my $notSameStrand = check_bit_from_base_10($flag, 4);
      if ($notSameStrand == 1) {
        $read_info_ref->{$read_id}{'seq'} = comp_rev($seq);
      } else {
        $read_info_ref->{$read_id}{'seq'} = $seq;
      }
    }
  }
  close $sam_handle;
}

sub print_read3_lines {
  my ($pending_lines_ref, $sam_list3_handle, $summary_data_ref) = @_;

  for my $line (@{$pending_lines_ref}) {
    print $sam_list3_handle $line;
    $summary_data_ref->{'number_of_reads_output'} ++;
  }
  @{$pending_lines_ref} = ();
}

sub print_read3_lines_for_previous_group {
  my ($group, $pending_group, $pending_lines_ref, $sam_list3_handle,
      $summary_data_ref) = @_;

  if ((defined $pending_group) and ($pending_group < $group)) {
    print_read3_lines($pending_lines_ref, $sam_list3_handle, $summary_data_ref);
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
  my ($group_id_header_i, $chr, $read_start, $chr_group_ref, $group_all_info_ref,
      $pending_group, $possible_groups_ref, $pending_lines_ref, $columns_ref,
      $sam_list3_handle, $summary_data_ref) = @_;

  my $group = $possible_groups_ref->[0];
  $columns_ref->[$group_id_header_i] = $group;
  my $out_line = join("\t", @{$columns_ref});
  $out_line .= "\n";

  my $num_possible = scalar(@{$possible_groups_ref});
  if ($num_possible > 2) {
    my $message = ("a read is assigned to more than 2 read groups"
                   . " ($num_possible) $out_line");
    die $message;
  }

  if ($num_possible == 2) {
    my $group_2 = $possible_groups_ref->[1];
    $pending_group = print_read3_lines_for_previous_group(
      $group_2, $pending_group, $pending_lines_ref, $sam_list3_handle,
      $summary_data_ref);

    print $sam_list3_handle $out_line;
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
    $pending_group = print_read3_lines_for_previous_group(
      $group + 1, $pending_group, $pending_lines_ref, $sam_list3_handle,
      $summary_data_ref);
    print $sam_list3_handle $out_line;
    $summary_data_ref->{'number_of_reads_output'} ++;
  } else {
    print_read3_lines_for_previous_group(
      $group, $pending_group, $pending_lines_ref, $sam_list3_handle,
      $summary_data_ref);
    push @{$pending_lines_ref}, $out_line;
    $pending_group = $group;
  }

  return $pending_group;
}

sub load_sj_simplified_for_chr {
  my ($chr, $out, $SJ_updated_all_ref) = @_;

  %{$SJ_updated_all_ref} = ();
  my $sj_simplified_path;
  if (length($chr) <= 5) {
    $sj_simplified_path = "$out/${chr}_SJ_simplified.list";
  } else {
    $sj_simplified_path = "$out/other_SJ_simplified.list";
  }
  open(my $sj_simplified_handle, '<', $sj_simplified_path)
    or die "cannot open tmp $sj_simplified_path: $!";
  while (<$sj_simplified_handle>) {
    chomp;
    my @line = split /\t/;
    if ($line[0] eq 'SJ_cluster') {
      next;
    }
    my $sj = $line[1];
    my $is_high_confidence = $line[-2];
    my $by_start_index = $line[-1];
    my %details;
    $details{'is_high_confidence'} = $is_high_confidence;
    $details{'by_start_index'} = $by_start_index;
    $SJ_updated_all_ref->{$sj} = \%details;
  }

  close $sj_simplified_handle;
}

sub update_sj_cor_seq {
  my ($sep_SJ, $output_titles_ID_ref, $line_ref, $SJ_updated_all_ref) = @_;

  if ($line_ref->[$output_titles_ID_ref->{'SJcorSeqRef'}] eq 'NA') {
    return;
  }

  my @strings = split(',', $line_ref->[$output_titles_ID_ref->{'SJcorSeqRef'}]);
  for my $i (0 .. $#strings) {
    my $string = $strings[$i];
    my @info = split(';', $string);
    my $sj = $info[1];
    my $sj_with_strand;
    if (exists $SJ_updated_all_ref->{$sj.$sep_SJ.'0'}) {
      $sj_with_strand = $sj.$sep_SJ.'0';
    } elsif (exists $SJ_updated_all_ref->{$sj.$sep_SJ.'1'}) {
      $sj_with_strand = $sj.$sep_SJ.'1';
    } elsif (exists $SJ_updated_all_ref->{$sj.$sep_SJ.'x'}) {
      $sj_with_strand = $sj.$sep_SJ.'x';
    } else {
      die "did not find $sj in $_";
    }

    $info[1] = $sj_with_strand;
    $info[2] = $SJ_updated_all_ref->{$sj_with_strand}{'is_high_confidence'};
    push @info, $SJ_updated_all_ref->{$sj_with_strand}{'by_start_index'};
    $strings[$i] = join(';', @info);
  }

  $line_ref->[$output_titles_ID_ref->{'SJcorSeqRef'}] = join(',', @strings);
}

sub output_read_lines_with_updated_sj_info_and_sequence {
  my ($sam_list_path, $sam_list3_path, $gtf_read_groups, $out, $num, $sep_SJ,
      $group_sep_in_all_ref, $group_all_info_ref, $chr_group_ref,
      $output_titles_ID_ref, $read_info_ref, $summary_data_ref) = @_;

  open(my $sam_list_handle, '<', $sam_list_path)
    or die "cannot open tmp $sam_list_path: $!";
  open(my $sam_list3_handle, '>', $sam_list3_path)
    or die "cannot write tmp $sam_list3_path: $!";
  my $prev_chr;
  my %SJ_updated_all = ();
  my @pending_read3_lines = ();
  my $pending_read3_lines_group = undef;
  while (<$sam_list_handle>) {
    chomp;
    my @line = split /\t/;
    my $read_id = $line[$output_titles_ID_ref->{'readID'}];
    my $line_num = $line[$output_titles_ID_ref->{'line_num'}];
    if ($read_id eq 'readID') {
      next;
    }
    if (exists $read_info_ref->{$read_id}
        and ($read_info_ref->{$read_id}{'line_num'} != $line_num)) {
      next;
    }

    my $chr = $line[$output_titles_ID_ref->{'chr'}];
    my $orig_group_id = $line[$output_titles_ID_ref->{'group_ID'}];
    my $read_start = $line[$output_titles_ID_ref->{'start'}];
    my $read_end = $line[$output_titles_ID_ref->{'end'}];
    my @possible_groups = lookup_final_group(
      $gtf_read_groups, $num, $orig_group_id, $chr, $read_start, $read_end,
      $group_sep_in_all_ref, $group_all_info_ref, $chr_group_ref);
    if (scalar(@possible_groups) == 0) {
      $summary_data_ref->{'num_reads_filtered_no_assigned_read_group'} ++;
      next;
    }

    if (!defined($prev_chr) or ($chr ne $prev_chr)) {
      print_read3_lines(\@pending_read3_lines,
                        $sam_list3_handle, $summary_data_ref);
      $pending_read3_lines_group = undef;

      load_sj_simplified_for_chr($chr, $out, \%SJ_updated_all);
    }

    $prev_chr = $chr;
    update_sj_cor_seq($sep_SJ, $output_titles_ID_ref, \@line, \%SJ_updated_all);

    if ($line[$output_titles_ID_ref->{'readSeq'}] eq 'NA') {
      if (exists $read_info_ref->{$read_id}
          and exists $read_info_ref->{$read_id}{'seq'}) {
        my $read_seq = $read_info_ref->{$read_id}{'seq'};
        if ($read_seq eq 'short') {
          $summary_data_ref->{'num_reads_filtered_missing_full_sequence'} ++;
          next;
        }
        $line[$output_titles_ID_ref->{'readSeq'}] = $read_seq;
      } else {
        $line[$output_titles_ID_ref->{'readSeq'}] = 'read_not_recorded';
      }
    }

    $pending_read3_lines_group = maybe_print_read3_lines(
      $output_titles_ID_ref->{'group_ID'}, $chr, $read_start, $chr_group_ref,
      $group_all_info_ref, $pending_read3_lines_group, \@possible_groups,
      \@pending_read3_lines, \@line, $sam_list3_handle, $summary_data_ref);
  }
  close $sam_list_handle;

  print_read3_lines(\@pending_read3_lines, $sam_list3_handle, $summary_data_ref);
  close $sam_list3_handle;
  $pending_read3_lines_group = undef;
}

sub parallel_scan2 {
  my ($file, $is_sam, $input_ref, $summary_data_ref) = @{$_[0]};

  my $file_ID_ref = $input_ref->{'file_ID'};
  my $out = $input_ref->{'out'};
  my $output_titles_ID_ref = $input_ref->{'output_titles_ID'};
  my $group_sep_in_all_ref = $input_ref->{'group_sep_in_all'};
  my $group_all_info_ref = $input_ref->{'group_all_info'};
  my $chr_group_ref = $input_ref->{'chr_group'};
  my $gtf_read_groups = $input_ref->{'gtf_read_groups'};
  my $sep_SJ = $input_ref->{'sep_SJ'};
  my $keep_tmp = $input_ref->{'keep_tmp'};
  my $sort_buffer_size = $input_ref->{'sort_buffer_size'};

  my $num = $file_ID_ref->{$file};
  print "$file\t$num\n";

  my $sam_list_path = "$out/$num/sam.list";
  my $sam_list2_path = "$out/$num/sam.list2";
  my $sam_list3_path = "$out/$num/sam.list3";
  # sort by readID
  my $sort_command = ("sort $sort_buffer_size -k 3,3 $sam_list_path"
                      . " > $sam_list2_path");
  my $exit_sort = system($sort_command);
  if ($exit_sort != 0) {
    die "Failed to $sort_command. Exit code is $exit_sort";
  }

  my %read_info;
  my $any_short = (
    load_info_for_reads_with_multiple_alignments_or_missing_sequence(
      $sam_list2_path, $output_titles_ID_ref, \%read_info));

  if ($any_short) {
    check_sam_for_full_sequence($is_sam, $file, \%read_info);
  }

  output_read_lines_with_updated_sj_info_and_sequence(
    $sam_list_path, $sam_list3_path, $gtf_read_groups, $out, $num, $sep_SJ,
    $group_sep_in_all_ref, $group_all_info_ref, $chr_group_ref,
    $output_titles_ID_ref, \%read_info, $summary_data_ref);

  if (!defined $keep_tmp) {
    unlink $sam_list_path;
    unlink $sam_list2_path;
  }
}

sub parallel_scan_prep {
  my ($main_input_ref, $summary_data_ref) = @_;

  my $files_ref = thaw($main_input_ref->{'files_for_thread'});
  my $all_input_ref = thaw($main_input_ref->{'inputs_for_all_threads'});
  my $scan_num = $all_input_ref->{'scan_number'};

  for my $file (@{$files_ref}) {
    my $is_sam = ($file =~ /sam$/i);
    if ($scan_num == 1) {
      parallel_scan1($file, $is_sam, $all_input_ref, $summary_data_ref);
    } else {
      parallel_scan2([$file, $is_sam, $all_input_ref, $summary_data_ref]);
    }
  }
}

sub thread_work_loop {
  my ($input_queue, $output_queue) = @{$_[0]};

  # Allow the main thread to stop other threads with kill('TERM')
  $SIG{'TERM'} = sub { threads->exit(); };

  while (defined(my $work_details = $input_queue->dequeue())) {
    my %summary_data = ();
    parallel_scan_prep($work_details, \%summary_data);
    my $serialized_summary :shared;
    $serialized_summary = freeze(\%summary_data);
    # signal that work is done
    $output_queue->enqueue($serialized_summary);
  }
  $output_queue->end();
}

sub start_threads {
  my ($num_thread) = @_;

  my @worker_threads = ();
  my @thread_input_queues = ();
  my @thread_output_queues = ();

  for (1 .. $num_thread) {
    push @thread_input_queues, Thread::Queue->new();
    push @thread_output_queues, Thread::Queue->new();
    my $thread = threads->new(
      {'context' => 'void'}, \&thread_work_loop,
      [$thread_input_queues[-1], $thread_output_queues[-1]]);
    push @worker_threads, $thread;
  }

  my %thread_details = ();
  $thread_details{'worker_threads'} = \@worker_threads;
  $thread_details{'thread_input_queues'} = \@thread_input_queues;
  $thread_details{'thread_output_queues'} = \@thread_output_queues;
  return \%thread_details;
}

sub maybe_create_out_dir {
  my ($out, $warn_ref) = @_;

  if (!-d $out) {
    push @{$warn_ref}, "Work directory $out does not exist. Make it by myself\n";
    mkdir($out) or die "cannot mkdir $out: $!";
  }
}

sub initialize_summary_data {
  my ($summary_data_ref) = @_;

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
  $summary_data_ref->{'num_alignments_filtered_for_chrM'} = 0;
  $summary_data_ref->{'num_alignments_filtered_for_secondary'} = 0;
  $summary_data_ref->{'num_alignments_filtered_for_mapping_quality'} = 0;
  $summary_data_ref->{'num_alignments_filtered_for_max_insertion'} = 0;
  $summary_data_ref->{'num_alignments_filtered_for_past_chr_end'} = 0;
  $summary_data_ref->{'num_reads_filtered_missing_full_sequence'} = 0;
  $summary_data_ref->{'num_reads_filtered_no_assigned_read_group'} = 0;
}

sub check_samtools_version {
  my ($die_ref) = @_;

  my $samtools_version_error = ESPRESSO_Version::check_samtools_version();
  if ($samtools_version_error ne '') {
    push @{$die_ref}, "$samtools_version_error\n";
  }
}

sub check_storable_version {
  my ($warn_ref) = @_;

  my $storable_version_warning = ESPRESSO_Version::check_storable_version();
  if ($storable_version_warning ne '') {
    push @{$warn_ref}, "$storable_version_warning\n";
  }
}

sub trim_whitespace {
  my ($str) = @_;

  $str =~ s/^\s+|\s+$//g;
  return $str;
}

sub read_list_samples {
  my ($list_samples, $sam_size_ref, $file_ID_ref, $sam_sample_ref,
      $die_ref) = @_;

  my $file_i = 0;
  open(my $list_handle, '<', $list_samples)
    or die "cannot open $list_samples: $!";
  while(<$list_handle>) {
    my $line_string = trim_whitespace($_);
    if ($line_string eq '') {
      next;
    }
    my @line = split(/\t/, $line_string);
    my $sam_path = $line[0];
    my $sample_name = $line[1];

    if ((exists $sam_sample_ref->{$sam_path})
        and ($sam_sample_ref->{$sam_path} ne $sample_name)) {
      my $message = "$sam_path in $list_samples corresponds to multiple samples";
      $message .= ": $sam_sample_ref->{$sam_path}, $sample_name!\n";
      push @{$die_ref}, $message;
    } elsif (-f $sam_path and -r $sam_path) {
      $sam_size_ref->{$sam_path} = -s $sam_path;
      $sam_sample_ref->{$sam_path} = $sample_name;

      $file_ID_ref->{$sam_path} = $file_i;
      $file_i ++;
    } else {
      my $message = "$sam_path in $list_samples does not exist";
      $message .= " or is not readable!\n";
      push @{$die_ref}, $message;
    }
    if ($sam_path !~ /[sb]am$/i) {
      push @{$die_ref}, "Don't know what format $sam_path is!\n";
    }
  }
  close $list_handle;
}

sub print_with_timestamp {
  my ($message) = @_;

  my $time_value = localtime;
  print "[$time_value] $message\n";
}

sub assign_sam_files_to_threads {
  my ($sam_size_ref, $num_thread) = @_;

  my $num_sam_files = scalar(keys %{$sam_size_ref});
  if ($num_sam_files > 1) {
    my $message = "Calculating how to assign $num_sam_files files into";
    $message .= "$num_thread threads";
    print_with_timestamp($message);
  }

  $num_thread = min($num_sam_files, $num_thread);
  my @ordered_file = sort {
    $sam_size_ref->{$b} <=> $sam_size_ref->{$a};
  } keys %{$sam_size_ref};
  my @file_distributed;
  my @thread_sum;
  my $last_thread_i = $num_thread - 1;
  for my $i (0 .. $last_thread_i) {
    $thread_sum[$i] = 0;
  }

  for my $file (@ordered_file) {
    my ($min_thread, undef) = sort {
      $thread_sum[$a] <=> $thread_sum[$b];
    } (0 .. $last_thread_i);
    unshift @{$file_distributed[$min_thread]}, $file;
    $thread_sum[$min_thread] += $sam_size_ref->{$file};
  }

  my @sort_thread = sort {
    $thread_sum[$b] <=> $thread_sum[$a]
  } (0 .. $last_thread_i);
  my @final_distribution = map {$file_distributed[$_]} (@sort_thread);
  return \@final_distribution;
}

sub file_name_from_path {
  my ($path) = @_;

  my $last_slash_i = rindex($path, '/');
  return substr($path, $last_slash_i + 1);
}

sub write_samples_updated_file {
  my ($orig_path, $out_dir, $file_ID_ref, $sam_sample_ref) = @_;

  my $samples_file_name = file_name_from_path($orig_path);
  my $samples_updated_path = "$out_dir/${samples_file_name}.updated";
  open(my $samples_updated_handle, '>', $samples_updated_path)
    or die "cannot write $samples_updated_path: $!";
  while (my ($file, $ID) = each %{$file_ID_ref}) {
    write_tsv_columns($samples_updated_handle,
                      [$file, $sam_sample_ref->{$file}, $ID]);
  }
  close $samples_updated_handle;
}

sub remove_old_output_files {
  my ($out_dir) = @_;

  opendir(my $dir_handle, $out_dir) or die "cannot opendir $out_dir: $!";
  for my $file (readdir($dir_handle)) {
    my $path = "$out_dir/$file";
    if (-d $path and $file =~ /^\d+$/) {
      opendir(my $dir2_handle, $path) or die "cannot opendir $path: $!";
      for my $file2 (readdir($dir2_handle)) {
        my $path2 = "$path/$file2";
        if (-f $path2 and (starts_with($file, 'sam.list')
                           or ($file eq 'sj.list')
                           or ($file eq 'group.list'))) {
          unlink $path2;
        }
      }
      closedir($dir2_handle);
    } elsif (-f $path and (($file eq 'SJ_group_all.fa')
                           or ($file =~ /SJ_simplified\.list$/))) {
      unlink $path;
    }
  }
  closedir($dir_handle);
}

sub load_ref_fasta {
  my ($path, $chr_seq_ref, $chr_seq_len_ref) = @_;

  my $current_chr;
  open(my $handle, '<', $path) or die "cannot open $path: $!";
  while (my $line = <$handle>) {
    $line =~ s/\r\n//;
    chomp($line);
    if ($line =~ /^>/) {
      my $info = substr($line, 1);
      ($current_chr, undef) = split(/\s+/, $info);
    } elsif (defined $current_chr) {
      $chr_seq_ref->{$current_chr} .= uc($line);
      $chr_seq_len_ref->{$current_chr} += length($line);
    }
  }
  close($handle);
}

sub load_annotation_file {
  my ($anno, $fa, $chr_seq_ref, $strand_digit_ref, $sep_SJ, $gtf_read_groups,
      $anno_coords_by_chr_by_gene_ref, $anno_SJ_ref, $summary_data_ref) = @_;

  my (%anno_exons, %isoform_info, %anno_chr_names);
  open(my $anno_handle, '<', $anno ) or die "cannot open $anno: $!";
  while(<$anno_handle>) {
    chomp;
    my @line = split /\t/;
    next if $line[2] ne 'exon';
    my ($current_isoform, $current_gene);

    if ($line[8] =~ /transcript_id \"(\S+)\"/) {
      $current_isoform = $1;
    } else {
      die "no transcript_id found in $_";
    }
    if ($line[8] =~ /gene_id \"(\S+)\"/) {
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
        my $gene_coords_ref = (
          $anno_coords_by_chr_by_gene_ref->{$chr}{$current_gene});
        $gene_coords_ref->[0] = min($gene_coords_ref->[0], $start);
        $gene_coords_ref->[1] = max($gene_coords_ref->[1], $end);
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
  my $num_chrs_in_anno = scalar(keys %anno_chr_names);
  $num_chrs_only_in_anno = $num_chrs_in_anno - $num_chrs_in_anno_and_fa;

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
    for my $i (1 .. $#exon_start_sort) {
      my $prev_exon_end = $exon_start_ref->{$exon_start_sort[$i-1]};
      my $exon_start = $exon_start_sort[$i];
      my $SJ2_no_strand = $chr.$sep_SJ.$prev_exon_end.$sep_SJ.$exon_start;
      my $SJ2 = $SJ2_no_strand.$sep_SJ.$strand_digit_ref->{$strand};
      if (!exists $anno_SJ_ref->{$chr}{$SJ2}) {
        my %details;
        $details{'strand'} = $strand_digit_ref->{$strand};
        $details{'start'} = $prev_exon_end;
        $details{'end'} = $exon_start;
        $anno_SJ_ref->{$chr}{$SJ2} = \%details;
      }
    }
  }
}

sub load_sj_bed_file {
  my ($SJ_bed, $strand_digit_ref, $chr_seq_len_ref, $sep_SJ, $anno_SJ_ref) = @_;

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
      my $message = 'Downstream splice site coordinate should be less than';
      $message .= " chromosome length($chr_seq_len_ref->{$chr} <= $sj_end): $_";
      die $message;
    }
    if ($sj_start <= 0) {
      my $message = 'Upstream splice site coordinate should be positive';
      $message .= "($sj_start <= 0): $_";
      die $message;
    }

    my $SJ2 = $chr.$sep_SJ.$sj_start.$sep_SJ.$sj_end.$sep_SJ.$strand;
    if (!exists $anno_SJ_ref->{$chr}{$SJ2}) {
      my %details;
      $details{'strand'} = $strand_digit_ref->{$line[5]};
      $details{'start'} = $sj_start;
      $details{'end'} = $sj_end;
      $anno_SJ_ref->{$chr}{$SJ2} = \%details;
    }
  }
  close $bed_handle;
}

sub load_reference_files {
  my ($fa_path, $anno_path, $bed_path, $sep_SJ, $gtf_read_groups, $chr_seq_ref,
      $chr_seq_len_ref, $anno_SJ_ref, $anno_coords_by_chr_by_gene_ref,
      $chr_sort_ref, $summary_data_ref) = @_;

  my %strand_digit = ('+' => 0, '-' => 1);

  load_ref_fasta($fa_path, $chr_seq_ref, $chr_seq_len_ref);
  my $num_chr_seqs = scalar(keys %{$chr_seq_ref});
  if ($num_chr_seqs == 0) {
    die "No sequence information found in $fa_path";
  }

  # sort chr names to get consistent output
  @{$chr_sort_ref} = sort {$a cmp $b} keys %{$chr_seq_ref};

  # $anno_SJ{$chr}{$SJ} = {'strand': strand, 'start': start, 'end': end}
  # where $SJ = chr:start:end:strand
  #       $start is the 0-based 1st position of the junction
  #       $end is the 1-based last position of the junction
  if (defined $anno_path) {
    print_with_timestamp('Loading annotation');
    load_annotation_file(
      $anno_path, $fa_path, $chr_seq_ref, \%strand_digit, $sep_SJ,
      $gtf_read_groups, $anno_coords_by_chr_by_gene_ref, $anno_SJ_ref,
      $summary_data_ref);
  } else {
    $summary_data_ref->{'num_chrs_only_in_fa'} = $num_chr_seqs;
  }

  if (defined $bed_path) {
    print_with_timestamp('Loading custom splice junction');
    load_sj_bed_file($bed_path, \%strand_digit, $chr_seq_len_ref, $sep_SJ,
                     $anno_SJ_ref);
  }

  for my $chr (keys %{$anno_SJ_ref}) {
    my $num_sjs = scalar(keys %{$anno_SJ_ref->{$chr}});
    $summary_data_ref->{'num_annotated_splice_junctions'} += $num_sjs;
  }
}

sub get_output_titles {
  my ($output_titles_ref, $output_titles_ID_ref) = @_;

  @{$output_titles_ref} = (
    'group_ID', 'line_num', 'readID', 'sample', 'read_length', 'flag', 'chr',
    'start', 'mapq', 'end', 'notSameStrand', 'mappedGenome', 'clip_ends',
    'exonIntronRef', 'IDS_SJ_ref', 'NM_num', 'insNumBg', 'delNumBg',
    'substNumBg', 'totalNumBg', 'SJcorSeqRef', 'readSeq');
  my $num_titles = scalar(@{$output_titles_ref});
  for my $i (0 .. ($num_titles - 1)) {
    $output_titles_ID_ref->{$output_titles_ref->[$i]} = $i;
  }
}

sub start_threads_for_parallel_scan1 {
  my ($sep_SJ, $out_dir, $mapq_cutoff, $chrM, $cont_del_max,
      $inserted_cont_cutoff, $group_extra_range, $SJFS_dist, $file_ID_ref,
      $chr_seq_ref, $chr_seq_len_ref, $sam_sample_ref, $sam_distr_ref,
      $thread_details_ref) = @_;

  my %inputs_for_all_threads;
  $inputs_for_all_threads{'scan_number'} = 1;
  $inputs_for_all_threads{'file_ID'} = $file_ID_ref;
  $inputs_for_all_threads{'out'} = $out_dir;
  $inputs_for_all_threads{'mapq_cutoff'} = $mapq_cutoff;
  $inputs_for_all_threads{'chrM'} = $chrM;
  $inputs_for_all_threads{'cont_del_max'} = $cont_del_max;
  $inputs_for_all_threads{'chr_seq_len'} = $chr_seq_len_ref;
  $inputs_for_all_threads{'inserted_cont_cutoff'} = $inserted_cont_cutoff;
  $inputs_for_all_threads{'group_extra_range'} = $group_extra_range;
  $inputs_for_all_threads{'SJFS_dist'} = $SJFS_dist;
  $inputs_for_all_threads{'chr_seq'} = $chr_seq_ref;
  $inputs_for_all_threads{'sep_SJ'} = $sep_SJ;
  $inputs_for_all_threads{'sam_sample'} = $sam_sample_ref;
  start_threads_with_inputs(\%inputs_for_all_threads, $sam_distr_ref,
                            $thread_details_ref);
}

sub start_threads_for_parallel_scan2 {
  my ($sep_SJ, $out_dir, $gtf_read_groups, $keep_tmp, $sort_buffer_size,
      $file_ID_ref, $output_titles_ID_ref, $group_sep_in_all_ref,
      $group_all_info_ref, $chr_group_ref, $sam_distr_ref,
      $thread_details_ref) = @_;

  my %inputs_for_all_threads;
  $inputs_for_all_threads{'scan_number'} = 2;
  $inputs_for_all_threads{'file_ID'} = $file_ID_ref;
  $inputs_for_all_threads{'out'} = $out_dir;
  $inputs_for_all_threads{'output_titles_ID'} = $output_titles_ID_ref;
  $inputs_for_all_threads{'group_sep_in_all'} = $group_sep_in_all_ref;
  $inputs_for_all_threads{'group_all_info'} = $group_all_info_ref;
  $inputs_for_all_threads{'chr_group'} = $chr_group_ref;
  $inputs_for_all_threads{'gtf_read_groups'} = $gtf_read_groups;
  $inputs_for_all_threads{'sep_SJ'} = $sep_SJ;
  $inputs_for_all_threads{'keep_tmp'} = $keep_tmp;
  $inputs_for_all_threads{'sort_buffer_size'} = $sort_buffer_size;
  start_threads_with_inputs(\%inputs_for_all_threads, $sam_distr_ref,
                            $thread_details_ref);
}

sub start_threads_with_inputs {
  my ($inputs_for_all_threads_ref, $sam_distr_ref, $thread_details_ref) = @_;

  my $serialized_inputs_for_all_threads :shared;
  $serialized_inputs_for_all_threads = freeze($inputs_for_all_threads_ref);

  for my $thread_i (0 .. $#{$sam_distr_ref}) {
    my $files_for_thread = $sam_distr_ref->[$thread_i];
    my $serialized_files_for_thread :shared;
    $serialized_files_for_thread = freeze($files_for_thread);
    my %thread_input :shared;
    $thread_input{'files_for_thread'} = $serialized_files_for_thread;
    $thread_input{'inputs_for_all_threads'} = $serialized_inputs_for_all_threads;
    $thread_details_ref->{'thread_input_queues'}[$thread_i]->enqueue(
      \%thread_input);
    print "Worker $thread_i begins to scan: \n @{$files_for_thread}\n";
  }
}

sub merge_summary_results {
  my ($summary_data_ref, $thread_result) = @_;

  my $thread_summary_ref = thaw($thread_result);
  while (my ($key, $value) = each %{$thread_summary_ref}) {
    $summary_data_ref->{$key} += $value;
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
    die 'Exiting due to error in worker thread';
  }
}

sub wait_for_threads_to_complete_work {
  my ($num_threads, $thread_input_queues_ref, $thread_output_queues_ref,
      $worker_threads_ref, $summary_data_ref) = @_;

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
        merge_summary_results($summary_data_ref, $thread_result);
        print "Worker $thread_i finished reporting.\n";
      } elsif ($thread_is_joinable) {
        # The thread should not be joinable until the work queue is ended.
        # This thread must have had an error.
        cleanup_threads_and_exit_if_error($thread_input_queues_ref,
                                          $worker_threads_ref);
      } else {
        push @still_running_threads, $thread_i;
      }
    }

    @running_threads = @still_running_threads;
    @still_running_threads = ();
  }
}

sub sort_numeric_string {
  # for this function, disable the warning:
  #   "isn't numeric in numeric comparison"
  no warnings 'numeric';
  # prefer to sort numerically, but fall back to string comparison
  return (($a <=> $b) or ($a cmp $b));
}

sub load_read_groups_from_sams {
  my ($out_dir, $keep_tmp, $sam_sort_ref, $file_ID_ref, $group_sep_ref) = @_;

  for my $sam_file (@{$sam_sort_ref}) {
    my $num = $file_ID_ref->{$sam_file};
    my $group_path = "$out_dir/$num/group.list";
    open(my $group_handle, '<', $group_path)
      or die "cannot open $group_path: $!";
    while (<$group_handle>) {
      chomp;
      my @line = split /\t/;
      my $orig_group_id = $line[0];
      my $chr = $line[1];
      my $start = $line[2];
      my $end = $line[3];
      $group_sep_ref->{$chr}{$num.'_'.$orig_group_id} = [$start, $end];
    }
    close $group_handle;
    if (!defined $keep_tmp) {
      unlink $group_path;
    }
  }
}

sub create_coordinate_groups_for_genes {
  my ($gene_group_i_offset, $anno_coords_by_gene_ref,
      $anno_gene_groups_ref) = @_;

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
      next;
    }

    my $old_group_coords_ref = $anno_gene_groups_ref->{$gene_group_i};
    if ($gene_start > $old_group_coords_ref->[1]) {
      $gene_group_i ++;
      $anno_gene_groups_ref->{$gene_group_i} = [$gene_start, $gene_end];
    } elsif ($gene_end > $old_group_coords_ref->[1]) {
      $anno_gene_groups_ref->{$gene_group_i}->[1] = $gene_end;
    }
  }
}

sub finalize_read_groups {
  my ($gtf_read_groups, $group_sep_ref, $chr_sort_ref,
      $anno_coords_by_chr_by_gene_ref, $group_all_info_ref, $chr_group_ref,
      $group_sep_in_all_ref) = @_;

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
      create_coordinate_groups_for_genes($gene_group_i_offset,
                                         $anno_coords_by_gene_ref,
                                         \%anno_gene_groups);
      my $highest_coord_used_for_chr = (
        $orig_group_ref->{$orig_group_sort[-1]}[1]);
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
          if (($group_start <= $orig_group_start)
              and ($orig_group_start <= $group_end)) {
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
      my $first_group_ref = $orig_group_ref->{$first_group_id};
      push @{$group_all_info_ref}, [$first_group_ref->[0],
                                    $first_group_ref->[1]];
      $group_sep_in_all_ref->{$first_group_id} = $gene_group_i;
      push @{$chr_group_ref->{$chr}}, $gene_group_i;
      for my $group_id (@orig_group_sort) {
        my $group_start = $orig_group_ref->{$group_id}[0];
        my $group_end = $orig_group_ref->{$group_id}[1];
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
  my ($gtf_read_groups, $file_num, $orig_group_id, $chr, $start, $end,
      $group_sep_in_all_ref, $group_all_info_ref, $chr_group_ref) = @_;

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

sub load_sjs_from_sams {
  my ($out_dir, $gtf_read_groups, $sam_sort_ref, $file_ID_ref,
      $group_sep_in_all_ref, $group_all_info_ref, $chr_group_ref,
      $SJ_all_info_ref, $SJ_group_ref) = @_;

  for my $sam_file (@{$sam_sort_ref}) {
    my $num = $file_ID_ref->{$sam_file};
    my $sj_path = "$out_dir/$num/sj.list";
    open(my $sj_handle, '<', $sj_path) or die "cannot open tmp $sj_path: $!";
    while (<$sj_handle>) {
      chomp;
      my @line = split /\t/;
      my $group_id = $line[0];
      my $sj = $line[1];
      my $chr = $line[2];
      my $start = $line[3];
      my $end = $line[4];
      my $num_perfect = $line[5];
      my $num_all = $line[6];

      if (exists $SJ_all_info_ref->{$sj}) {
        $SJ_all_info_ref->{$sj}{'perfect_read_count'} += $num_perfect;
        $SJ_all_info_ref->{$sj}{'read_count'} += $num_all;
      } else {
        my @possible_groups = lookup_final_group(
          $gtf_read_groups, $num, $group_id, $chr, $start, $end,
          $group_sep_in_all_ref, $group_all_info_ref, $chr_group_ref);
        if (scalar(@possible_groups) > 0) {
          my %info;
          $info{'chr'} = $chr;
          $info{'start'} = $start;
          $info{'end'} = $end;
          $info{'perfect_read_count'} = $num_perfect;
          $info{'read_count'} = $num_all;
          $SJ_all_info_ref->{$sj} = \%info;
        }
        for my $n (@possible_groups) {
          push @{$SJ_group_ref->{$n}}, $sj;
        }
      }
    }
    close $sj_handle;
  }
}

sub determine_sj_details {
  my ($chr, $n, $sep_SJ, $splicing_signals_strand_ref, $chr_seq_ref,
      $chr_sj_ref, $sorted_sjs_group_ref, $SJ_group_ref, $SJ_all_info_ref,
      $SJ_updated_ref) = @_;

  my %recorded_SJ;
  for my $SJ (@{$SJ_group_ref->{$n}}) {
    my $sj_info_ref = $SJ_all_info_ref->{$SJ};
    my $chr2 = $sj_info_ref->{'chr'};
    my $SJ_start = $sj_info_ref->{'start'};
    my $SJ_end = $sj_info_ref->{'end'};
    my $perfect_count = $sj_info_ref->{'perfect_read_count'};
    my $all_count = $sj_info_ref->{'read_count'};
    if ($chr ne $chr2) {
      die "$chr ne $chr2: $SJ";
    }
    my $upstream_2nt = substr($chr_seq_ref->{$chr}, $SJ_start, 2);
    my $downstream_2nt = substr($chr_seq_ref->{$chr}, $SJ_end - 2, 2);
    my $tag = 1;
    my %notSameStrand_SJ_freq;
    if (exists $splicing_signals_strand_ref->{$upstream_2nt.$downstream_2nt}) {
      my $strand = $splicing_signals_strand_ref->{$upstream_2nt.$downstream_2nt};
      $notSameStrand_SJ_freq{$strand} ++;
    }

    my $sj_0_strand = $SJ . $sep_SJ . '0';
    my $sj_1_strand = $SJ . $sep_SJ . '1';
    if (exists $chr_sj_ref->{$sj_0_strand}) {
      $tag = 2;
      $notSameStrand_SJ_freq{'0'} ++;
      $recorded_SJ{$sj_0_strand} = 1;
    }
    if (exists $chr_sj_ref->{$sj_1_strand}) {
      $tag = 2;
      $notSameStrand_SJ_freq{'1'} ++;
      $recorded_SJ{$sj_1_strand} = 1;
    }

    my $notSameStrand_SJ;
    if (exists $notSameStrand_SJ_freq{'1'}
        and exists $notSameStrand_SJ_freq{'0'}) {
      $notSameStrand_SJ = 'x';
    } elsif (exists $notSameStrand_SJ_freq{'1'}) {
      $notSameStrand_SJ = '1';
    } elsif (exists $notSameStrand_SJ_freq{'0'}) {
      $notSameStrand_SJ = '0';
    } else {
      $notSameStrand_SJ = 'x';
      $tag = 0;
    }

    my $SJ2 = $SJ . $sep_SJ . $notSameStrand_SJ;
    my %sj_details;
    $sj_details{'start'} = $SJ_start;
    $sj_details{'end'} = $SJ_end;
    $sj_details{'strand'} = $notSameStrand_SJ;
    $sj_details{'perfect_read_count'} = $perfect_count;
    $sj_details{'read_count'} = $all_count;
    $sj_details{'up_2nt'} = $upstream_2nt;
    $sj_details{'down_2nt'} = $downstream_2nt;
    $sj_details{'tag'} = $tag;
    $sj_details{'has_read'} = 'yes';
    if (exists $recorded_SJ{$sj_0_strand} or $recorded_SJ{$sj_1_strand}) {
      $sj_details{'is_annotated'} = 'yes';
    } else {
      $sj_details{'is_annotated'} = 'no';
    }
    $sj_details{'is_high_confidence'} = 0;
    $SJ_updated_ref->{$SJ2} = \%sj_details;
  }

  for my $SJ2 (@{$sorted_sjs_group_ref}) {
    if (!exists $recorded_SJ{$SJ2}) {
      my $sj2_details = $chr_sj_ref->{$SJ2};
      my $notSameStrand_SJ = $sj2_details->{'strand'};
      my $SJ_start = $sj2_details->{'start'};
      my $SJ_end = $sj2_details->{'end'};

      my %sj_details;
      $sj_details{'start'} = $SJ_start;
      $sj_details{'end'} = $SJ_end;
      $sj_details{'strand'} = $notSameStrand_SJ;
      $sj_details{'perfect_read_count'} = 0;
      $sj_details{'read_count'} = 0;
      $sj_details{'up_2nt'} = 'TBD';
      $sj_details{'down_2nt'} = 'TBD';
      $sj_details{'tag'} = 2;
      $sj_details{'has_read'} = 'no';
      $sj_details{'is_annotated'} = 'yes';
      $sj_details{'is_high_confidence'} = 0;
      $SJ_updated_ref->{$SJ2} = \%sj_details;
    }
  }
}

sub create_sj_clusters {
  my ($two_sjfs_dist_plus_add, $SJ_updated_ref, $sorted_SJ_cluster_ref,
      $sorted_SJ_cluster_ends_ref) = @_;

  my @sorted_SJ = sort {
    $SJ_updated_ref->{$a}{'start'} <=> $SJ_updated_ref->{$b}{'start'};
  } keys %{$SJ_updated_ref};
  my @sorted_SJ_group_f = ([$sorted_SJ[0]]);
  for my $i (1 .. $#sorted_SJ) {
    my $sj = $sorted_SJ[$i];
    my $sj_details = $SJ_updated_ref->{$sj};
    my $prev_sj = $sorted_SJ[$i - 1];
    my $prev_sj_details = $SJ_updated_ref->{$prev_sj};
    my $start_diff = $sj_details->{'start'} - $prev_sj_details->{'start'};
    if ($start_diff > $two_sjfs_dist_plus_add) {
      push @sorted_SJ_group_f, [$sj];
    } else {
      push @{$sorted_SJ_group_f[-1]}, $sj;
    }
  }

  for my $SJ_group_ref (@sorted_SJ_group_f) {
    my @sorted_SJ_group_r = sort {
      $SJ_updated_ref->{$a}{'end'} <=> $SJ_updated_ref->{$b}{'end'};
    } @{$SJ_group_ref};
    my $first_sj = $sorted_SJ_group_r[0];
    my $first_sj_details = $SJ_updated_ref->{$first_sj};
    push @{$sorted_SJ_cluster_ref}, [$first_sj];
    push @{$sorted_SJ_cluster_ends_ref}, [$first_sj_details->{'start'},
                                          $first_sj_details->{'end'}];
    for my $i (1 .. $#sorted_SJ_group_r) {
      my $sj = $sorted_SJ_group_r[$i];
      my $prev_sj = $sorted_SJ_group_r[$i-1];
      my $sj_details = $SJ_updated_ref->{$sj};
      my $prev_sj_details = $SJ_updated_ref->{$prev_sj};
      my $end_diff = $sj_details->{'end'} - $prev_sj_details->{'end'};
      if ($end_diff > $two_sjfs_dist_plus_add) {
        push @{$sorted_SJ_cluster_ref}, [$sj];
        push @{$sorted_SJ_cluster_ends_ref}, [$sj_details->{'start'},
                                              $sj_details->{'end'}];
      } else {
        push @{$sorted_SJ_cluster_ref->[-1]}, $sj;
        $sorted_SJ_cluster_ends_ref->[-1][1] = $sj_details->{'end'};
        if ($sorted_SJ_cluster_ends_ref->[-1][0] > $sj_details->{'start'}) {
          $sorted_SJ_cluster_ends_ref->[-1][0] = $sj_details->{'start'};
        }
      }
    }
  }
}

sub write_sj_files {
  my ($chr, $n, $read_num_cutoff, $read_ratio_cutoff, $sjfs_dist_plus_add,
      $chr_seq_ref, $sorted_SJ_cluster_ref, $sorted_SJ_cluster_ends_ref,
      $SJ_updated_ref, $sj_has_been_counted_in_a_group_ref, $simplified_handle,
      $all_handle, $summary_data_ref) = @_;

  my @SJ_cluster_r_index_sort = sort {
    my $end_a = $sorted_SJ_cluster_ends_ref->[$a][1];
    my $end_b = $sorted_SJ_cluster_ends_ref->[$b][1];
    $end_a <=> $end_b;
  } (0 .. $#{$sorted_SJ_cluster_ref});

  for my $by_end_index (0 .. $#SJ_cluster_r_index_sort) {
    my $by_start_index = $SJ_cluster_r_index_sort[$by_end_index];
    my @SJs_current_cluster = @{$sorted_SJ_cluster_ref->[$by_start_index]};
    my $current_ends = $sorted_SJ_cluster_ends_ref->[$by_start_index];
    write_tsv_columns($simplified_handle,
                      ['SJ_cluster', $n, $by_start_index, $by_end_index, $chr,
                       $current_ends->[0], $current_ends->[1]]);
    for my $SJ2 (@SJs_current_cluster) {
      my $sj_details = $SJ_updated_ref->{$SJ2};
      $sj_details->{'by_start_index'} = $by_start_index;

      my $already_seen = exists $sj_has_been_counted_in_a_group_ref->{$SJ2};
      $sj_has_been_counted_in_a_group_ref->{$SJ2} = 1;

      my $SJ_start = $sj_details->{'start'};
      my $SJ_end = $sj_details->{'end'};
      my $perfect_count = $sj_details->{'perfect_read_count'};
      my $all_count = $sj_details->{'read_count'};
      my $tag = $sj_details->{'tag'};
      if (!$already_seen) {
        $summary_data_ref->{'perfect_splice_junction_read_count'} += (
          $perfect_count);
        $summary_data_ref->{'total_splice_junction_read_count'} += (
          $all_count);
      }
      my $has_count = $perfect_count >= $read_num_cutoff;
      my $has_ratio = $perfect_count >= ($all_count * $read_ratio_cutoff);
      if ($tag == 2 or ($tag == 1 and $has_count and $has_ratio)) {
        if (!$already_seen) {
          $summary_data_ref->{'num_high_confidence_splice_junctions'} ++;
        }
        $sj_details->{'is_high_confidence'} = 1;
        print $all_handle ">$SJ2 SJclst:$by_start_index: group:$n:\n";
        my $sj_seq = substr($chr_seq_ref->{$chr},
                            $SJ_start - $sjfs_dist_plus_add,
                            $sjfs_dist_plus_add);
        $sj_seq .= substr($chr_seq_ref->{$chr},
                          $SJ_end, $sjfs_dist_plus_add);
        print $all_handle $sj_seq . "\n";
      }
      # TODO empty column
      write_tsv_columns($simplified_handle,
                        [$n, $SJ2, '', $SJ_start, $SJ_end,
                         $sj_details->{'strand'}, $perfect_count, $all_count,
                         $sj_details->{'up_2nt'}, $sj_details->{'down_2nt'},
                         $tag, $sj_details->{'has_read'},
                         $sj_details->{'is_annotated'},
                         $sj_details->{'is_high_confidence'}, $by_start_index]);
    }
  }
}

sub summarize_sj_info_and_write_files {
  my ($out_dir, $SJFS_dist, $SJFS_dist_add, $read_num_cutoff, $read_ratio_cutoff,
      $sep_SJ, $chr_group_ref, $anno_SJ_ref, $SJ_group_ref, $SJ_all_info_ref,
      $group_all_info_ref, $chr_seq_ref, $summary_data_ref) = @_;

  my $two_sjfs_dist_plus_add = 2 * $SJFS_dist + $SJFS_dist_add;
  my $sjfs_dist_plus_add = $SJFS_dist + $SJFS_dist_add;
  my %splicing_signals_strand = ('GTAG' => '0', 'GCAG' => '0', 'ATAC' => '0',
                                 'CTAC' => '1', 'CTGC' => '1', 'GTAT' => '1');
  my $sj_all_path = "$out_dir/SJ_group_all.fa";
  open(my $all_handle, '>', $sj_all_path) or die "cannot write $sj_all_path: $!";
  while (my ($chr, $group_ID_ref) = each %{$chr_group_ref}) {
    my $chr_sj_ref = $anno_SJ_ref->{$chr};
    my $simplified_handle;
    if (length($chr) <= 5) {
      my $sj_simplified_path = "$out_dir/${chr}_SJ_simplified.list";
      open($simplified_handle, '>', $sj_simplified_path)
        or die "cannot open tmp $sj_simplified_path: $!";
    } else {
      my $sj_simplified_path = "$out_dir/other_SJ_simplified.list";
      open($simplified_handle, '>>', $sj_simplified_path)
        or die "cannot open tmp $sj_simplified_path: $!";
    }

    my @sorted_sjs_chr = sort {
      my $a_info = $chr_sj_ref->{$a};
      my $a_start = $a_info->{'start'};
      my $a_end = $a_info->{'end'};
      my $b_info = $chr_sj_ref->{$b};
      my $b_start = $b_info->{'start'};
      my $b_end = $b_info->{'end'};
      ($a_start <=> $b_start) or ($a_end <=> $b_end);
    } keys %{$chr_sj_ref};

    my %sj_has_been_counted_in_a_group = ();
    my $start_SJ_index = 0;
    for my $n (@{$group_ID_ref}) {
      my @sorted_sjs_group;
      my $group_start = $group_all_info_ref->[$n][0];
      my $group_end = $group_all_info_ref->[$n][1];
      for my $i ($start_SJ_index .. $#sorted_sjs_chr) {
        my $sj = $sorted_sjs_chr[$i];
        my $sj_info = $chr_sj_ref->{$sj};
        my $sj_start = $sj_info->{'start'};
        my $sj_end = $sj_info->{'end'};
        if ($sj_start > $group_end) {
          last;
        } elsif ($sj_end < $group_start) {
          $start_SJ_index = $i;
          next;
        }
        push @sorted_sjs_group, $sj;
      }

      my %SJ_updated;
      determine_sj_details($chr, $n, $sep_SJ, \%splicing_signals_strand,
                           $chr_seq_ref, $chr_sj_ref, \@sorted_sjs_group,
                           $SJ_group_ref, $SJ_all_info_ref, \%SJ_updated);
      if (scalar(keys %SJ_updated) == 0) {
        next;
      }

      my @sorted_SJ_cluster;
      my @sorted_SJ_cluster_ends;
      create_sj_clusters($two_sjfs_dist_plus_add, \%SJ_updated,
                         \@sorted_SJ_cluster, \@sorted_SJ_cluster_ends);
      write_sj_files($chr, $n, $read_num_cutoff, $read_ratio_cutoff,
                     $sjfs_dist_plus_add, $chr_seq_ref, \@sorted_SJ_cluster,
                     \@sorted_SJ_cluster_ends, \%SJ_updated,
                     \%sj_has_been_counted_in_a_group, $simplified_handle,
                     $all_handle, $summary_data_ref);
    }
    close $simplified_handle;
  }
  close $all_handle;
}

sub write_summary_file {
  my ($summary_path, $summary_data_ref) = @_;

  open(my $summary_handle, '>', $summary_path)
    or die "cannot write $summary_path: $!";
  print $summary_handle "$summary_data_ref->{'perl_command'}\n";
  print $summary_handle ('number of chromosomes only in input annotation:'
                         . " $summary_data_ref->{'num_chrs_only_in_anno'}\n");
  print $summary_handle ('number of chromosomes only in input FASTA:'
                         . " $summary_data_ref->{'num_chrs_only_in_fa'}\n");
  print $summary_handle ('number of chromosomes in both annotation and FASTA:'
                         . " $summary_data_ref->{'num_chrs_in_anno_and_fa'}\n");
  print $summary_handle ('number of isoforms in input annotation:'
                         . " $summary_data_ref->{'num_annotated_isoforms'}\n");
  print $summary_handle
    ('number of splice junctions in input annotation:'
     . " $summary_data_ref->{'num_annotated_splice_junctions'}\n");
  print $summary_handle
    ('number of high confidence splice junctions:'
     . " $summary_data_ref->{'num_high_confidence_splice_junctions'}\n");
  print $summary_handle
    ('total over all splice junctions of supporting reads:'
     . " $summary_data_ref->{'total_splice_junction_read_count'}\n");
  print $summary_handle
    ('total over all splice junctions of perfect reads:'
     . " $summary_data_ref->{'perfect_splice_junction_read_count'}\n");
  print $summary_handle ('number of read groups:'
                         . " $summary_data_ref->{'number_of_read_groups'}\n");
  print $summary_handle ('number of reads in output:'
                         . " $summary_data_ref->{'number_of_reads_output'}\n");
  print $summary_handle
    ('number of chrM alignments filtered:'
     . " $summary_data_ref->{'num_alignments_filtered_for_chrM'}\n");
  print $summary_handle
    ('number of secondary alignments filtered:'
     . " $summary_data_ref->{'num_alignments_filtered_for_secondary'}\n");
  print $summary_handle
    ('number of alignments filtered for mapping quality:'
     . " $summary_data_ref->{'num_alignments_filtered_for_mapping_quality'}\n");
  print $summary_handle
    ('number of alignments filtered for a long insertion:'
     . " $summary_data_ref->{'num_alignments_filtered_for_max_insertion'}\n");
  print $summary_handle
    ('number of alignments filtered for unrecognized coordinates:'
     . " $summary_data_ref->{'num_alignments_filtered_for_past_chr_end'}\n");
  print $summary_handle
    ('number of reads filtered for missing full sequence:'
     . " $summary_data_ref->{'num_reads_filtered_missing_full_sequence'}\n");
  print $summary_handle
    ('number of reads filtered for not being matched to a read group:'
     . " $summary_data_ref->{'num_reads_filtered_no_assigned_read_group'}\n");
  close $summary_handle;
}

sub main {
  my $args = parse_args();
  if (defined($args->{'help'})) {
    show_help_message();
    return;
  }
  if (!defined($args->{'list_samples'}) or !defined($args->{'fa'})) {
    print "The following parameter(s) are required:\n";
    if (!defined($args->{'list_samples'})) {
      print "\t--list_samples/-L";
    } if (!defined($args->{'fa'})) {
      print "\t--fa/-F";
    }
    print "\nPlease use the --help or -H option to get usage information.\n";
    return;
  }

  my $thread_details_ref = start_threads($args->{'num_thread'});

  my @warn_reason;
  my @die_reason;
  maybe_create_out_dir($args->{'out'}, \@warn_reason);

  # The summary file will be written as a final step
  my $summary_path = "$args->{'out'}/espresso_s_summary.txt";
  my %summary_data = ();
  initialize_summary_data(\%summary_data);
  $summary_data{'perl_command'} = ("$^X " . __FILE__
                                   . " $args->{'arguments_before_parsing'}");

  check_samtools_version(\@die_reason);
  check_storable_version(\@warn_reason);

  my %sam_size = ();
  my %file_ID = ();
  my %sam_sample = ();
  read_list_samples($args->{'list_samples'}, \%sam_size, \%file_ID, \%sam_sample,
                    \@die_reason);

  if (!$args->{'alignment_read_groups'} and !$args->{'anno'}) {
    push @die_reason, "--anno is required unless --alignment_read_groups\n";
  }

  if (@warn_reason >= 1) {
    print @warn_reason;
  }
  if (@die_reason >= 1) {
    die @die_reason;
  }

  # --gtf_read_groups was added to the code as an optional feature.
  # Now it's the default and --alignment_read_groups enables the old behavior.
  # $alignment_read_groups is only used for argument parsing.
  # $gtf_read_groups is used in the remainder of the code.
  $args->{'gtf_read_groups'} = !$args->{'alignment_read_groups'};
  delete $args->{'alignment_read_groups'};

  my $sam_distr_ref = assign_sam_files_to_threads(\%sam_size,
                                                  $args->{'num_thread'});
  write_samples_updated_file($args->{'list_samples'}, $args->{'out'}, \%file_ID,
                             \%sam_sample);
  remove_old_output_files($args->{'out'});

  my $sep_SJ = ':';
  print_with_timestamp("Loading reference");
  my %chr_seq;
  my %chr_seq_len;
  my %anno_SJ;
  my %anno_coords_by_chr_by_gene;
  my @chr_sort;
  load_reference_files(
    $args->{'fa'}, $args->{'anno'}, $args->{'SJ_bed'}, $sep_SJ,
    $args->{'gtf_read_groups'}, \%chr_seq, \%chr_seq_len, \%anno_SJ,
    \%anno_coords_by_chr_by_gene, \@chr_sort, \%summary_data);

  my @output_titles;
  my %output_titles_ID;
  get_output_titles(\@output_titles, \%output_titles_ID);

  start_threads_for_parallel_scan1(
    $sep_SJ, $args->{'out'}, $args->{'mapq_cutoff'}, $args->{'chrM'},
    $args->{'cont_del_max'}, $args->{'inserted_cont_cutoff'},
    $args->{'group_extra_range'}, $args->{'SJFS_dist'}, \%file_ID, \%chr_seq,
    \%chr_seq_len, \%sam_sample, $sam_distr_ref, $thread_details_ref);
  wait_for_threads_to_complete_work(
    scalar(@{$sam_distr_ref}), $thread_details_ref->{'thread_input_queues'},
    $thread_details_ref->{'thread_output_queues'},
    $thread_details_ref->{'worker_threads'}, \%summary_data);

  print_with_timestamp("Re-cluster all reads");
  my @sam_sort = sort sort_numeric_string keys %sam_size;
  my %group_sep;
  load_read_groups_from_sams($args->{'out'}, $args->{'keep_tmp'}, \@sam_sort,
                             \%file_ID, \%group_sep);

  my @group_all_info;
  my %chr_group;
  my %group_sep_in_all;
  my $num_read_groups = finalize_read_groups(
    $args->{'gtf_read_groups'}, \%group_sep, \@chr_sort,
    \%anno_coords_by_chr_by_gene, \@group_all_info, \%chr_group,
    \%group_sep_in_all);
  $summary_data{'number_of_read_groups'} = $num_read_groups;
  %anno_coords_by_chr_by_gene = ();
  %group_sep = ();

  my %SJ_all_info;
  my %SJ_group;
  load_sjs_from_sams($args->{'out'}, $args->{'gtf_read_groups'}, \@sam_sort,
                     \%file_ID, \%group_sep_in_all, \@group_all_info,
                     \%chr_group, \%SJ_all_info, \%SJ_group);

  print_with_timestamp(
    'Summarizing annotated splice junctions for each read group');
  summarize_sj_info_and_write_files(
    $args->{'out'}, $args->{'SJFS_dist'}, $args->{'SJFS_dist_add'},
    $args->{'read_num_cutoff'}, $args->{'read_ratio_cutoff'}, $sep_SJ,
    \%chr_group, \%anno_SJ, \%SJ_group, \%SJ_all_info, \@group_all_info,
    \%chr_seq, \%summary_data);
  %SJ_all_info = ();
  %SJ_group = ();
  %chr_seq = ();
  %anno_SJ = ();

  start_threads_for_parallel_scan2(
    $sep_SJ, $args->{'out'}, $args->{'gtf_read_groups'}, $args->{'keep_tmp'},
    $args->{'sort_buffer_size'}, \%file_ID, \%output_titles_ID,
    \%group_sep_in_all, \@group_all_info, \%chr_group, $sam_distr_ref,
    $thread_details_ref);
  @group_all_info = ();
  %chr_group = ();
  %group_sep_in_all = ();

  wait_for_threads_to_complete_work(
    scalar(@{$sam_distr_ref}), $thread_details_ref->{'thread_input_queues'},
    $thread_details_ref->{'thread_output_queues'},
    $thread_details_ref->{'worker_threads'}, \%summary_data);

  cleanup_threads_and_exit_if_error($thread_details_ref->{'thread_input_queues'},
                                    $thread_details_ref->{'worker_threads'});
  write_summary_file($summary_path, \%summary_data);

  print_with_timestamp('ESPRESSO_S finished its work.');
}

main();

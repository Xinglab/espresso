package Parasail;

use strict;
use warnings;

use File::Basename qw(dirname);
use Cwd qw(abs_path);

require Exporter;
require DynaLoader;
our @ISA = qw(Exporter DynaLoader);
our @EXPORT_OK = qw(cigar_decode_len
                    cigar_decode_op
                    matrix_create
                    matrix_set_value
                    profile_create_16
                    sw_trace_striped_profile_16
                    result_get_cigar
                    get_cigar_len
                    get_cigar_query_start
                    get_cigar_ref_start
                    get_cigar_int_i
                    cigar_free
                    result_get_score
                    result_get_end_query
                    result_get_end_ref
                    result_free
                    profile_free
                    matrix_free
                    get_score_matrix_p
                    result_get_cigar_info);

our $VERSION = '0.01';

my $file_dir = dirname(__FILE__);
$file_dir = abs_path($file_dir);

my @found_parasail = DynaLoader::dl_findfile('-lparasail');
if (scalar(@found_parasail) == 0) {
  # If the parasail library isn't found in the default search, then
  # check in the directory of this file and load the library.
  # If it was found in the default search then it *should* get
  # loaded during the regular bootstrap process.
  my @found_parasail = DynaLoader::dl_findfile("-L$file_dir", '-lparasail');
  if (scalar(@found_parasail) == 1) {
    DynaLoader::dl_load_file($found_parasail[0]);
  }
}

DynaLoader::bootstrap("Parasail", $VERSION);

sub get_score_matrix_p {
  my ($match, $mismatch) = @_;
  my $alphabet = 'ACGTN';
  my $ambiguous_base = 'N';
  my $score_matrix_p = matrix_create($alphabet, $match, $mismatch);
  for my $alpha_i (0 .. (length($alphabet) - 1)) {
    my $alpha = substr($alphabet, $alpha_i, 1);
    if ($alpha eq $ambiguous_base) {
      for my $alpha_j (0 .. (length($alphabet) - 1)) {
        matrix_set_value($score_matrix_p, $alpha_i, $alpha_j, 0);
        matrix_set_value($score_matrix_p, $alpha_j, $alpha_i, 0);
      }
    }
  }
  return $score_matrix_p;
}

sub result_get_cigar_info {
  my ($result_p, $query_seq, $ref_seq, $score_matrix_p) = @_;

  my $cigar_p = result_get_cigar($result_p,
                                 $query_seq,
                                 length($query_seq),
                                 $ref_seq,
                                 length($ref_seq),
                                 $score_matrix_p);
  my @cigar_parts = ();
  my $cigar_len = get_cigar_len($cigar_p);
  my $query_start = get_cigar_query_start($cigar_p);
  my $ref_start = get_cigar_ref_start($cigar_p);
  for my $cigar_i (0 .. ($cigar_len - 1)) {
    my $cigar_int = get_cigar_int_i($cigar_p, $cigar_i);
    my $cigar_num = cigar_decode_len($cigar_int);
    my $cigar_op = cigar_decode_op($cigar_int);
    push @cigar_parts, [$cigar_num, $cigar_op];
  }
  cigar_free($cigar_p);

  # It seems that the traceback for smith waterman can give an extra
  # operation at the start of the cigar even though that operation
  # didn't contribute to the alignment score.
  if ($cigar_len >= 1) {
    my $first_num = $cigar_parts[0][0];
    my $first_op = $cigar_parts[0][1];
    if ($first_op eq 'I') {
      shift(@cigar_parts);
      $query_start += $first_num;
    }
    if ($first_op eq 'D') {
      shift(@cigar_parts);
      $ref_start += $first_num
    }
    if ($first_op eq 'S') {
      shift(@cigar_parts);
      $query_start += $first_num;
    }
  }

  my %cigar_info;
  $cigar_info{'parts'} = \@cigar_parts;
  $cigar_info{'query_start'} = $query_start;
  $cigar_info{'ref_start'} = $ref_start;
  return %cigar_info;
}

1;

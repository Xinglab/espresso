use strict;
use warnings;

use Test::More tests => 17;

BEGIN { use_ok('Parasail') };

# Expect to align like:
#   ACGTAACCGGTTTGCA
#     GTA--CGCTTT
#   3=2D2=1X3=
my $ref_seq = "ACGTAACCGGTTTGCA";
my $query_seq = "GTACGCTTT";
my $match = 5;
my $mismatch = -4;
my $gap_open = 8;
my $gap_extension = 6;
my $score_matrix_p = Parasail::get_score_matrix_p($match, $mismatch);
my $profile_p = Parasail::profile_create_16($query_seq,
                                            length($query_seq),
                                            $score_matrix_p);
my $result_p = Parasail::sw_trace_striped_profile_16($profile_p,
                                                     $ref_seq,
                                                     length($ref_seq),
                                                     $gap_open,
                                                     $gap_extension);
my $score = Parasail::result_get_score($result_p);
is($score, 22);

my %cigar_info = Parasail::result_get_cigar_info($result_p, $query_seq,
                                                 $ref_seq, $score_matrix_p);
my @cigar_parts = @{$cigar_info{'parts'}};
my $query_start = $cigar_info{'query_start'};
my $ref_start = $cigar_info{'ref_start'};

is(scalar(@cigar_parts), 5);
is($query_start, 0);
is($ref_start, 2);
is($cigar_parts[0][0], 3);
is($cigar_parts[0][1], '=');
is($cigar_parts[1][0], 2);
is($cigar_parts[1][1], 'D');
is($cigar_parts[2][0], 2);
is($cigar_parts[2][1], '=');
is($cigar_parts[3][0], 1);
is($cigar_parts[3][1], 'X');
is($cigar_parts[4][0], 3);
is($cigar_parts[4][1], '=');

my $query_end = Parasail::result_get_end_query($result_p);
my $ref_end = Parasail::result_get_end_ref($result_p);
is($query_end, 8);
is($ref_end, 12);

Parasail::result_free($result_p);
Parasail::profile_free($profile_p);
Parasail::matrix_free($score_matrix_p);

#define PERL_NO_GET_CONTEXT
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include <stdint.h>
#include "parasail.h"

MODULE = Parasail		PACKAGE = Parasail

PROTOTYPES: DISABLE

TYPEMAP: <<HERE
uint32_t T_UV
int32_t T_IV
parasail_matrix_t* T_PTR
parasail_profile_t* T_PTR
parasail_result_t* T_PTR
parasail_cigar_t* T_PTR
HERE

uint32_t
cigar_decode_len(cigar_int)
    uint32_t cigar_int
  CODE:
    RETVAL = parasail_cigar_decode_len(cigar_int);
  OUTPUT:
    RETVAL

char
cigar_decode_op(cigar_int)
    uint32_t cigar_int
  CODE:
    RETVAL = parasail_cigar_decode_op(cigar_int);
  OUTPUT:
    RETVAL

parasail_matrix_t*
matrix_create(alphabet, match, mismatch)
    char* alphabet
    int32_t match
    int32_t mismatch
  CODE:
    RETVAL = parasail_matrix_create(alphabet, match, mismatch);
  OUTPUT:
    RETVAL

void
matrix_set_value(matrix, row, col, value)
    parasail_matrix_t* matrix
    int32_t row
    int32_t col
    int32_t value
  CODE:
    parasail_matrix_set_value(matrix, row, col, value);

parasail_profile_t*
profile_create_16(s1, s1Len, matrix)
    char* s1
    int32_t s1Len
    parasail_matrix_t* matrix
  CODE:
    RETVAL = parasail_profile_create_16(s1, s1Len, matrix);
  OUTPUT:
    RETVAL

parasail_result_t*
sw_trace_striped_profile_16(profile, s2, s2Len, open, gap)
    parasail_profile_t* profile
    char* s2
    int32_t s2Len
    int32_t open
    int32_t gap
  CODE:
    RETVAL = parasail_sw_trace_striped_profile_16(profile, s2, s2Len, open, gap);
  OUTPUT:
    RETVAL

parasail_cigar_t*
result_get_cigar(result, seqA, lena, seqB, lenb, matrix)
    parasail_result_t* result
    char* seqA
    int32_t lena
    char* seqB
    int32_t lenb
    parasail_matrix_t* matrix
  CODE:
    RETVAL = parasail_result_get_cigar(result, seqA, lena, seqB, lenb, matrix);
  OUTPUT:
    RETVAL

int32_t
get_cigar_len(cigar_p)
    parasail_cigar_t* cigar_p
  CODE:
    RETVAL = cigar_p->len;
  OUTPUT:
    RETVAL

int32_t
get_cigar_query_start(cigar_p)
    parasail_cigar_t* cigar_p
  CODE:
    RETVAL = cigar_p->beg_query;
  OUTPUT:
    RETVAL

int32_t
get_cigar_ref_start(cigar_p)
    parasail_cigar_t* cigar_p
  CODE:
    RETVAL = cigar_p->beg_ref;
  OUTPUT:
    RETVAL

uint32_t
get_cigar_int_i(cigar_p, cigar_i)
    parasail_cigar_t* cigar_p
    int32_t cigar_i
  CODE:
    RETVAL = cigar_p->seq[cigar_i];
  OUTPUT:
    RETVAL

void
cigar_free(cigar)
    parasail_cigar_t* cigar
  CODE:
    parasail_cigar_free(cigar);

int32_t
result_get_score(result)
    parasail_result_t* result
  CODE:
    RETVAL = parasail_result_get_score(result);
  OUTPUT:
    RETVAL

int32_t
result_get_end_query(result)
    parasail_result_t* result
  CODE:
    RETVAL = parasail_result_get_end_query(result);
  OUTPUT:
    RETVAL

int32_t
result_get_end_ref(result)
    parasail_result_t* result
  CODE:
    RETVAL = parasail_result_get_end_ref(result);
  OUTPUT:
    RETVAL

void
result_free(result)
    parasail_result_t* result
  CODE:
    parasail_result_free(result);

void
profile_free(profile)
    parasail_profile_t* profile
  CODE:
    parasail_profile_free(profile);

void
matrix_free(matrix)
    parasail_matrix_t* matrix
  CODE:
    parasail_matrix_free(matrix);

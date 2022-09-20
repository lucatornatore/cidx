
#include "cidx.h"


#if !defined(USE_LIBC_QSORT)
typedef struct { PID_t pid; int gen; } pidgen_t;
#define QS_NAME particle_id_gen
#define QS_DATA_ALIGN  ALIGN
#define QS_DATA_TYPE   particle_t
#define QS_DATA_SIZE   (sizeof( particle_t ))
#define KEY_DATA_TYPE  pidgen_t
#define KEY_COPY(P, A) { (P).pid = (A)->pid; (P).gen = (A)->gen; }
#define CMP_KEY(A, B)  (int)( (((A).pid > (B)->pid) - ((A).pid < (B)->pid)) + ((((A).pid > (B)->pid) - ((A).pid < (B)->pid))==0)*((A).gen - (B)->gen) )
#define CMP(A, B)      (int)( (((A)->pid > (B)->pid) - ((A)->pid < (B)->pid)) + ((((A)->pid > (B)->pid) - ((A)->pid < (B)->pid))==0)*((A)->gen - (B)->gen) )
#define QS_SWAP_TYPE 0
#include "qsort_template.h"


#define QS_NAME idtype_id
#define QS_DATA_ALIGN  ALIGN
#define QS_DATA_TYPE   pidtype_t
#define QS_DATA_SIZE   (sizeof( pidtype_t ))
#define KEY_DATA_TYPE  PID_t
#define KEY_COPY(P, A) { (P) = (A)->pid; }
#define CMP_KEY(A, B)  (((A) > (B)->pid) - ((A) < (B)->pid))
#define CMP(A, B)      (((A)->pid > (B)->pid) - ((A)->pid < (B)->pid))
#define QS_SWAP_TYPE 0
#include "qsort_template.h"


/* #define QS_NAME fofparticle_type */
/* #define QS_DATA_ALIGN  ALIGN */
/* #define QS_DATA_TYPE   particle_t */
/* #define QS_DATA_SIZE   (sizeof( particle_t )) */
/* #define KEY_DATA_TYPE  int */
/* #define KEY_COPY(P, A) { (P).type = (A)->type; (P).type = (A)->type; } */
/* #define CMP_KEY(A, B)  (int)( (((A).type > (B)->type) - ((A).type < (B)->type)) + ((((A).type > (B)->type) - ((A).type < (B)->type))==0)) */
/* #define CMP(A, B)      (int)( (((A)->type > (B)->type) - ((A)->type < (B)->type)) + ((((A)->type > (B)->type) - ((A)->type < (B)->type))==0)) */
/* #define QS_SWAP_TYPE 0 */
/* #include "qsort_template.h" */


#endif  // closes USE_LIBC_QSORT


#if !defined(USE_LIBC_BSEARCH)
#warning "using inline bsearch"

#include "mybsearch.h"
#endif

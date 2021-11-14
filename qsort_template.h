


// :: to be defined ::
// QS_DATA_TYPE
// QS_DATA_SIZE
// KEY_DATA_TYPE
// KEY_COPY ( P, A)
// CMP ( A, B )

#if !defined(QS_NAME)
#define QS_NAME
#endif
#define QS_ROUTINE_BASENAME inline_qsort_
#define QS_STACKNODE_BASENAME stack_node
#define QS_PASTE(x, y) x##y
#define QS_EVAL(x, y) QS_PASTE(x, y)
#define QS_ROUTINE_NAME QS_EVAL(QS_ROUTINE_BASENAME, QS_NAME)
#define QS_STACK_NODE QS_EVAL(QS_STACKNODE_BASENAME, QS_NAME)


#if !(defined QS_DATA_TYPE)
#error "you must define QS_DATA_TYPE before to include this header"
#endif
#if !(defined QS_DATA_SIZE)
#error "you must define QS_DATA_TYPE before to include this header"
#endif

#if !defined(QS_SWAP_TYPE)
#define QS_SWAP_TYPE 0
#endif

#if ( QS_SWAP_TYPE == 0 )

#warning "using direct-type unrolling when swapping"
#define QS_SWAP(A,B) do { QS_DATA_TYPE t = *(A); *(A)=*(B); *(B)=t; } while(0)

#elif ( QS_SWAP_TYPE == 1 )

/* Byte-wise swap two items of size SIZE. */
#warning "using sigle-byte unrolling when swapping"
#define QS_SWAP(A,B) do {int sz = (QS_DATA_SIZE);			\
    char *a = (char*)(A); char *b = (char*)(B);				\
    do { char _temp = *a;*a++ = *b;*b++ = _temp;}			\
    while (--sz);} while (0)

#elif ( QS_SWAP_TYPE == 4 )

#warning "using 4-bytes unrolling when swapping"
#define QS_DATA_SIZE_4   (QS_DATA_SIZE / 4)
#define QS_DATA_SIZE_4_r (QS_DATA_SIZE % 4)
#define QS_SWAP(A,B) do {char * restrict a = (char*)(__builtin_assume_aligned(A, QS_DATA_ALIGN)); \
    char * restrict b = (char*)(__builtin_assume_aligned(B, QS_DATA_ALIGN)); \
    _Pragma("ivdep")							\
      for( int sz = (QS_DATA_SIZE_4); sz--; ) {				\
	char _temp_a0 = *a, _temp_b0 = *b;				\
	char _temp_a1 = *(a+1), _temp_b1 = *(b+1);			\
	char _temp_a2 = *(a+2), _temp_b2 = *(b+2);			\
	char _temp_a3 = *(a+3), _temp_b3 = *(b+3);			\
	*a++ = _temp_b0;*b++ = _temp_a0;				\
	*a++ = _temp_b1;*b++ = _temp_a1;				\
	*a++ = _temp_b2;*b++ = _temp_a2;				\
	*a++ = _temp_b3;*b++ = _temp_a3;}				\
    for( int sz = QS_DATA_SIZE_4_r; sz--; ) {				\
      char _temp_a0 = *a, _temp_b0 = *b;				\
      *a++ = _temp_b0;*b++ = _temp_a0;}} while(0)

#elif ( QS_SWAP_TYPE == 8 )

#warning "using double unrolling when swapping"
#define QS_DATA_SIZE_8   (QS_DATA_SIZE / 8)
#define QS_DATA_SIZE_8_r (QS_DATA_SIZE % 8)
#define QS_SWAP(A,B) do {double * restrict ad = (double*)(__builtin_assume_aligned(A, QS_DATA_ALIGN)); \
    double * restrict bd = (double*)(__builtin_assume_aligned(B, QS_DATA_ALIGN)); \
    _Pragma("ivdep")							\
      for( int sz = (QS_DATA_SIZE_8); sz--; ) {				\
	double _temp_a = *ad, _temp_b = *bd;				\
	*ad++ = _temp_b; *bd++ = _temp_a;}				\
    char * restrict a = (char*)(--ad);					\
    char * restrict b = (char*)(--bd);					\
    _Pragma("ivdep")							\
      for( int sz = (QS_DATA_SIZE_8_r); sz--; ) {			\
	char _temp_a = *a, _temp_b = *b;				\
	*a++ = _temp_b;*b++ = _temp_a;}} while(0)

#endif

/* This should be replaced by a standard ANSI macro. */
#define QS_BYTES_PER_WORD 8

/* The next 4 #defines implement a very fast in-line stack abstraction. */
#define QS_STACK_SIZE     (QS_BYTES_PER_WORD * sizeof (double))
#define QS_PUSH(LOW,HIGH) do {top->lo = LOW;top++->hi = HIGH;} while (0)
#define QS_POP(LOW,HIGH)  do {LOW = (--top)->lo;HIGH = top->hi;} while (0)
#define QS_STACK_NOT_EMPTY (stack < top)

/* Discontinue quicksort algorithm when partition gets below this size.
   This particular magic number was chosen to work best on a Sun 4/260. */
#define QS_MAX_THRESH 8

/* Stack node declarations used to store unfulfilled partition obligations. */
typedef struct
{
  QS_DATA_TYPE *lo;
  QS_DATA_TYPE *hi;
} QS_STACK_NODE;

/* Order size using quicksort.  This implementation incorporates
   four optimizations discussed in Sedgewick:

   1. Non-recursive, using an explicit stack of pointer that store the
      next array partition to sort.  To save time, this maximum amount
      of space required to store an array of MAX_INT is allocated on the
      stack.  Assuming a 32-bit integer, this needs only 32 *
      sizeof (stack_node) == 136 bits.  Pretty cheap, actually.

   2. Chose the pivot element using a median-of-three decision tree.
      This reduces the probability of selecting a bad pivot value and
      eliminates certain extraneous comparisons.

   3. Only quicksorts TOTAL_ELEMS / MAX_THRESH partitions, leaving
      insertion sort to order the MAX_THRESH items within each partition.
      This is a big win, since insertion sort is faster for small, mostly
      sorted array segments.

   4. The larger of the two sub-partitions is always pushed onto the
      stack first, with the algorithm then concentrating on the
      smaller partition.  This *guarantees* no more than log (n)
      stack size is needed (actually O(1) in this case)! */

int
QS_ROUTINE_NAME (void* restrict vbase_ptr, num_t total_elems);
int
QS_ROUTINE_NAME (void* restrict vbase_ptr, num_t total_elems)
{
  /* Allocating SIZE bytes for a pivot buffer facilitates a better
     algorithm below since we can do comparisons directly on the pivot. */
  QS_DATA_TYPE  *base_ptr   = vbase_ptr;
  KEY_DATA_TYPE  pivot;

  if (total_elems > QS_MAX_THRESH)
    {
      QS_DATA_TYPE  *restrict lo = base_ptr;
      QS_DATA_TYPE  *restrict hi = lo + (total_elems - 1);
      QS_STACK_NODE  stack[QS_STACK_SIZE]; /* Largest size needed for 32-bit int!!! */
      QS_STACK_NODE *top = stack + 1;

      while (QS_STACK_NOT_EMPTY)
        {
	  QS_DATA_TYPE * restrict left_ptr;
	  QS_DATA_TYPE * restrict right_ptr;
	  QS_DATA_TYPE * restrict mid;
          {
            {
              /* Select median value from among LO, MID, and HI. Rearrange
                 LO and HI so the three values are sorted. This lowers the
                 probability of picking a pathological pivot value and
                 skips a comparison for both the LEFT_PTR and RIGHT_PTR. */

              mid = lo + ((hi - lo) >> 1);

              if (CMP (mid, lo) < 0)
                QS_SWAP (mid, lo);
              if (CMP (hi, mid) < 0)
                QS_SWAP (mid, hi);
              else
                goto jump_over;
              if (CMP (mid, lo) < 0)
                QS_SWAP (mid, lo);
            jump_over:
              KEY_COPY (pivot, mid);
            }

	    left_ptr = lo + 1;
	    //right_ptr = lo;
	    right_ptr = hi - 1;
	    
            /* Here's the famous ``collapse the walls'' section of quicksort.
               Gotta like those tight inner loops!  They are the main reason
               that this algorithm runs much faster than others. */
	    do
              {
                while (CMP_KEY (pivot, left_ptr) > 0)
                  left_ptr ++;
		
                while (CMP_KEY (pivot, right_ptr) < 0)
                  right_ptr --;
		
                if (left_ptr < right_ptr)
                  {
                    QS_SWAP (left_ptr, right_ptr);
		    int right_ptr_dec = (left_ptr != mid);
                    left_ptr += (right_ptr != mid);
                    right_ptr -= right_ptr_dec;
                  }
                else if (left_ptr == right_ptr)
                  {
                    left_ptr ++;
                    right_ptr --;
                    break;
                  }
              }
            while (left_ptr <= right_ptr);

          }

          /* Set up pointers for next iteration.  First determine whether
             left and right partitions are below the threshold size. If so,
             ignore one or both.  Otherwise, push the larger partition's
             bounds on the stack and continue sorting the smaller one. */

          if ((size_t)(right_ptr - lo) <= QS_MAX_THRESH)
            {
              if ((size_t)(hi - left_ptr) <= QS_MAX_THRESH) /* Ignore both small partitions. */
                QS_POP (lo, hi);
              else              /* Ignore small left partition. */
                lo = left_ptr;
            }
          else if ((size_t)(hi - left_ptr) <= QS_MAX_THRESH) /* Ignore small right partition. */
            hi = right_ptr;
          else if ((right_ptr - lo) > (hi - left_ptr)) /* Push larger left partition indices. */
            {
              QS_PUSH (lo, right_ptr);
              lo = left_ptr;
            }
          else                  /* Push larger right partition indices. */
            {
              QS_PUSH (left_ptr, hi);
              hi = right_ptr;
            }
        }
    }

  /* Once the BASE_PTR array is partially sorted by quicksort the rest
     is completely sorted using insertion sort, since this is efficient
     for partitions below MAX_THRESH size. BASE_PTR points to the beginning
     of the array to sort, and END_PTR points at the very last element in
     the array (*not* one beyond it!). */


#define QS_MIN(X,Y) ((X) < (Y) ? (X) : (Y))

  {
    QS_DATA_TYPE *end_ptr = base_ptr + (total_elems - 1);
    QS_DATA_TYPE *run_ptr;
    QS_DATA_TYPE *tmp_ptr = base_ptr;
    QS_DATA_TYPE *thresh  = QS_MIN (end_ptr, base_ptr + QS_MAX_THRESH);

    /* Find smallest element in first threshold and place it at the
       array's beginning.  This is the smallest array element,
       and the operation speeds up insertion sort's inner loop. */

    for (run_ptr = tmp_ptr + 1; run_ptr <= thresh; run_ptr++ )
      if (CMP (run_ptr, tmp_ptr) < 0)
        tmp_ptr = run_ptr;

    if (tmp_ptr != base_ptr)
      QS_SWAP (tmp_ptr, base_ptr);

    /* Insertion sort, running from left-hand-side up to `right-hand-side.'
       Pretty much straight out of the original GNU qsort routine. */

    for (run_ptr = base_ptr + 1; (tmp_ptr = ++run_ptr) <= end_ptr; )
      {
	tmp_ptr = run_ptr-1;
	while (CMP (run_ptr, tmp_ptr) < 0)
          tmp_ptr--;

        if ((++tmp_ptr) != run_ptr)
          {
            QS_DATA_TYPE *trav;

            for (trav = run_ptr + 1; --trav >= run_ptr;)
              {
                QS_DATA_TYPE c = *trav;
                QS_DATA_TYPE *hi, *lo;

                for (hi = lo = trav; (lo -= 1) >= tmp_ptr; hi = lo)
                  *hi = *lo;
                *hi = c;
              }
          }

      }
  }

  
  /* { */
  /*   // bose-nelson sorting network for MAX_THRESH elements */
  /*   QS_DATA_TYPE *end_ptr = base_ptr + (total_elems - 1); */
  /*   for(QS_DATA_TYPE * restrict run = base_ptr; run < end_ptr; run += QS_MAX_THRESH) */
  /*     { */
  /*      #if ( QS_MAX_THRESH == 2 ) */
	
  /* 	if(CMP(run, run + 1) > 0) */
  /* 	  QS_SWAP(run, run + 1); */

  /*      #elif ( QS_MAX_THRESH == 3 ) */
  /* 	if(CMP(run + 1, run + 2) > 0) */
  /* 	  QS_SWAP(run + 1, run + 2); */
  /* 	if(CMP(run, run + 2) > 0) */
  /* 	  QS_SWAP(run, run + 2); */
  /* 	if(CMP(run, run + 1) > 0) */
  /* 	  QS_SWAP(run, run + 1); */
	
  /*      #elif (QS_MAX_THRESH == 4 || QS_MAX_THRESH == 8) */
    
  /* 	if(CMP(run, run + 1) > 0) */
  /* 	  QS_SWAP(run, run + 1); */
  /* 	if(CMP(run + 2, run + 3) > 0) */
  /* 	  QS_SWAP(run + 2, run + 3); */
  /* 	if(CMP(run, run + 2) > 0) */
  /* 	  QS_SWAP(run, run + 2); */
  /* 	if(CMP(run + 1, run + 3) > 0) */
  /* 	  QS_SWAP(run + 1, run + 3); */
  /* 	if(CMP(run + 1, run + 2) > 0) */
  /* 	  QS_SWAP(run + 1, run + 2); */

  /*      #if ( QS_MAX_THRESH == 8 ) */
  /* 	if(CMP(run + 4, run + 5) > 0) */
  /* 	  QS_SWAP(run + 4, run + 5);	 */
  /* 	if(CMP(run + 6, run + 7) > 0) */
  /* 	  QS_SWAP(run + 6, run + 7); */
  /* 	if(CMP(run + 4, run + 6) > 0) */
  /* 	  QS_SWAP(run + 4, run + 6);	 */
  /* 	if(CMP(run + 5, run + 7) > 0) */
  /* 	  QS_SWAP(run + 5, run + 7); */
  /* 	if(CMP(run + 5, run + 6) > 0) */
  /* 	  QS_SWAP(run + 5, run + 6); */

  /* 	if(CMP(run, run + 4) > 0) */
  /* 	  QS_SWAP(run, run + 4); */
  /* 	if(CMP(run + 1, run + 5) > 0) */
  /* 	  QS_SWAP(run + 1, run + 5); */
  /* 	if(CMP(run + 1, run + 4) > 0) */
  /* 	  QS_SWAP(run + 1, run + 4);	 */
  /* 	if(CMP(run + 2, run + 6) > 0) */
  /* 	  QS_SWAP(run + 2, run + 6); */
  /* 	if(CMP(run + 3, run + 7) > 0) */
  /* 	  QS_SWAP(run + 3, run + 7); */
  /* 	if(CMP(run + 3, run + 6) > 0) */
  /* 	  QS_SWAP(run + 3, run + 6); */
  /* 	if(CMP(run + 2, run + 4) > 0) */
  /* 	  QS_SWAP(run + 2, run + 4); */
  /* 	if(CMP(run + 3, run + 5) > 0) */
  /* 	  QS_SWAP(run + 3, run + 5); */
  /* 	if(CMP(run + 3, run + 4) > 0) */
  /* 	  QS_SWAP(run + 3, run + 4); */

  /*      #endif */
  /*      #elif ( QS_MAX_THRESH == 5 ) */

  /* 	if(CMP(run, run + 1) > 0) */
  /* 	  QS_SWAP(run, run + 1); */
  /* 	if(CMP(run + 3, run + 4) > 0) */
  /* 	  QS_SWAP(run + 3, run + 4); */
  /* 	if(CMP(run + 2, run + 4) > 0) */
  /* 	  QS_SWAP(run + 2, run + 4); */
  /* 	if(CMP(run + 2, run + 3) > 0) */
  /* 	  QS_SWAP(run + 2, run + 3); */
  /* 	if(CMP(run, run + 3) > 0) */
  /* 	  QS_SWAP(run, run + 3); */
  /* 	if(CMP(run, run + 2) > 0) */
  /* 	  QS_SWAP(run, run + 2); */
  /* 	if(CMP(run + 1, run + 4) > 0) */
  /* 	  QS_SWAP(run + 1, run + 4); */
  /* 	if(CMP(run + 1, run + 3) > 0) */
  /* 	  QS_SWAP(run + 1, run + 3); */
  /* 	if(CMP(run + 1, run + 2) > 0) */
  /* 	  QS_SWAP(run + 1, run + 2); */

  /*      #elif ( QS_MAX_THRESH == 6 ) */

  /* 	if(CMP(run + 1, run + 2) > 0) */
  /* 	  QS_SWAP(run + 1, run + 2); */
  /* 	if(CMP(run, run + 2) > 0) */
  /* 	  QS_SWAP(run, run + 2); */
  /* 	if(CMP(run, run + 1) > 0) */
  /* 	  QS_SWAP(run, run + 1); */
  /* 	if(CMP(run + 4, run + 5) > 0) */
  /* 	  QS_SWAP(run + 4, run + 5); */
  /* 	if(CMP(run + 3, run + 5) > 0) */
  /* 	  QS_SWAP(run + 3, run + 5); */
  /* 	if(CMP(run + 3, run + 4) > 0) */
  /* 	  QS_SWAP(run + 3, run + 4); */
  /* 	if(CMP(run, run + 3) > 0) */
  /* 	  QS_SWAP(run, run + 3); */
  /* 	if(CMP(run + 1, run + 4) > 0) */
  /* 	  QS_SWAP(run + 1, run + 4); */
  /* 	if(CMP(run + 2, run + 5) > 0) */
  /* 	  QS_SWAP(run + 2, run + 5); */
  /* 	if(CMP(run + 2, run + 4) > 0) */
  /* 	  QS_SWAP(run + 2, run + 4); */
  /* 	if(CMP(run + 1, run + 3) > 0) */
  /* 	  QS_SWAP(run + 1, run + 3); */
  /* 	if(CMP(run + 2, run + 3) > 0) */
  /* 	  QS_SWAP(run + 2, run + 3); */

  /*      #elif ( QS_MAX_THRESH == 7 ) */

  /* 	if(CMP(run + 1, run + 2) > 0) */
  /* 	  QS_SWAP(run + 1, run + 2); */
  /* 	if(CMP(run, run + 2) > 0) */
  /* 	  QS_SWAP(run, run + 2); */
  /* 	if(CMP(run, run + 1) > 0) */
  /* 	  QS_SWAP(run, run + 1); */
  /* 	if(CMP(run + 3, run + 4) > 0) */
  /* 	  QS_SWAP(run + 3, run + 4); */
  /* 	if(CMP(run + 5, run + 6) > 0) */
  /* 	  QS_SWAP(run + 5, run + 6); */
  /* 	if(CMP(run + 3, run + 5) > 0) */
  /* 	  QS_SWAP(run + 3, run + 5); */
  /* 	if(CMP(run + 4, run + 6) > 0) */
  /* 	  QS_SWAP(run + 4, run + 6);        */
  /* 	if(CMP(run + 4, run + 5) > 0) */
  /* 	  QS_SWAP(run + 4, run + 5); */
  /* 	if(CMP(run, run + 4) > 0) */
  /* 	  QS_SWAP(run, run + 4); */
  /* 	if(CMP(run, run + 3) > 0) */
  /* 	  QS_SWAP(run, run + 3); */
  /* 	if(CMP(run + 1, run + 5) > 0) */
  /* 	  QS_SWAP(run + 1, run + 5); */
  /* 	if(CMP(run + 2, run + 6) > 0) */
  /* 	  QS_SWAP(run + 2, run + 6); */
  /* 	if(CMP(run + 2, run + 5) > 0) */
  /* 	  QS_SWAP(run + 2, run + 5); */
  /* 	if(CMP(run + 1, run + 3) > 0) */
  /* 	  QS_SWAP(run + 1, run + 3); */
  /* 	if(CMP(run + 2, run + 4) > 0) */
  /* 	  QS_SWAP(run + 2, run + 4); */
  /* 	if(CMP(run + 2, run + 3) > 0) */
  /* 	  QS_SWAP(run + 2, run + 3); */

  /* 	#endif */
  /*     } */
  /* } */

  return 1;
}

#undef QS_ROUTINE_BASENAME
#undef QS_PASTE
#undef QS_EVAL
#undef QS_ROUTINE_NAME
#undef QS_STACK_NODE
#undef QS_DATA_SIZE_8
#undef QS_DATA_SIZE_8_r
#undef QS_SWAP
#undef QS_BYTES_PER_WORD
#undef QS_STACK_SIZE
#undef QS_PUSH
#undef QS_POP
#undef QS_STACK_NOT_EMPTY
#undef QS_MAX_THRESH
#undef QS_MIN
#undef QS_NAME
#undef QS_DATA_SIZE
#undef QS_DATA_TYPE
#undef KEY_DATA_TYPE
#undef KEY_COPY
#undef CMP_KEY
#undef CMP
#undef SWAP_TYPE
#undef QS_SWAP_TYPE

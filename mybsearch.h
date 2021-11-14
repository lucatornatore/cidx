

inline num_t mybsearch_in_ids(const pidtype_t *data, const num_t N, const PID_t Key, int *err)
 {
   register num_t low = 0;
   register num_t high = N;
   register num_t mid;
   *err = 0;
   
   mid = (low + high) / 2;
   while(low <= high) {     

     __builtin_prefetch (&data[(low + mid - 1)/2].pid, 0, 1);
     __builtin_prefetch (&data[(mid + 1 + high)/2].pid, 0, 1);

     low  = ( data[mid].pid < Key ? mid + 1 : low  );
     high = ( data[mid].pid > Key ? mid - 1 : high );
     
     if(data[mid].pid == Key)
       return mid;

     mid = (low + high) / 2;
   }

   *err = -1;
   return 0;
 }


inline void* mybsearch_in_P(const particle_t *data, const num_t N, const PID_t Key )
 {
   if( data[0].pid > Key )
     return NULL;
   
   register num_t low = 0;
   register num_t high = N;
   register num_t mid;
   
   mid = (low + high) / 2;
  #if defined(MASKED_ID_DBG)
   if( Key == MASKED_ID_DBG )
     printf(">> (low, high, mid) = (%llu, %llu, %llu)\n", low, high, mid);
  #endif

   while(low <= high) {

     __builtin_prefetch (&data[(low + mid - 1)/2].pid, 0, 1);
     __builtin_prefetch (&data[(mid + 1 + high)/2].pid, 0, 1);

     low  = ( data[mid].pid < Key ? mid + 1 : low  );
     high = ( data[mid].pid > Key ? mid - 1 : high );
     
     if(data[mid].pid == Key)
       return (void*)(data+mid);

     mid = (low + high) / 2;
     
    #if defined(MASKED_ID_DBG)
     if( Key == MASKED_ID_DBG )
       printf(">> (low, high, mid) = (%llu, %llu, %llu)\n", low, high, mid);
    #endif

   }
   
   return NULL; 
 }

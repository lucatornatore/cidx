

inline ull_t mybsearch_in_ids(const pidtype_t *data, const ull_t N, const PID_t Key, int *err)
 {
   register ull_t low = 0;
   register ull_t high = N;
   register ull_t mid;
   *err = 0;
   
   mid = (low + high) / 2;
   while(low <= high) {     

     __builtin_prefetch (&data[(low + mid - 1)/2].pid, 0, 1);
     __builtin_prefetch (&data[(mid + 1 + high)/2].pid, 0, 1);
     
     if(data[mid].pid < Key)
       low = mid + 1; 
     else if(data[mid].pid > Key)
       high = mid - 1;
     else 
       return mid;

     mid = (low + high) / 2;
   }

   *err = -1;
   return 0;
 }


inline ull_t mybsearch_in_P(const particle_t *data, const ull_t N, const PID_t Key, int *err)
 {
   register ull_t low = 0;
   register ull_t high = N;
   register ull_t mid;
   *err = 0;
   
   mid = (low + high) / 2;
   while(low <= high) {     

     __builtin_prefetch (&data[(low + mid - 1)/2].pid, 0, 1);
     __builtin_prefetch (&data[(mid + 1 + high)/2].pid, 0, 1);
     
     if(data[mid].pid < Key)
       low = mid + 1; 
     else if(data[mid].pid > Key)
       high = mid - 1;
     else 
       return mid;

     mid = (low + high) / 2;
   }

   *err = -1;
   return 0;
 }

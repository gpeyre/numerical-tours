#include <stdio.h>
#include <stdlib.h>
#include "fib.h"

int main(void)
{
  
    struct fibheap *a;
    void *arr[10];
    int i;

    a = fh_makekeyheap();
	
    for (i=1 ; i < 8 ; i++)
     {
      arr[i]= fh_insertkey(a,0,(void *)i);
      printf("adding: 0 %d \n",i);
     }
     
     printf(" \n");

     fh_replacekey(a, arr[1],-2);
     fh_replacekey(a, arr[4],-3);
     fh_replacekey(a, arr[7],-5);
	  
     printf("Wert(minkey) %d\n",fh_minkey(a));
     printf("Knoten: %d\n\n", (int)fh_extractmin(a)); 
        
     fh_replacekey(a, arr[3],-2);
     fh_replacekey(a, arr[6],-1);
	
     printf("Wert(minkey) %d\n",fh_minkey(a));
     printf("Knoten: %d\n\n", (int)fh_extractmin(a));
	
     fh_replacekey(a, arr[1],-9);
     fh_replacekey(a, arr[5],-3);

     printf("Wert(minkey) %d\n",fh_minkey(a));
     printf("Knoten: %d\n\n", (int)fh_extractmin(a));

     fh_replacekey(a, arr[2],-4);
     fh_replacekey(a, arr[5],-5);
     fh_replacekey(a, arr[6],-3);
        	
     printf("Wert(minkey) %d\n",fh_minkey(a));
     printf("Knoten: %d\n\n", (int)fh_extractmin(a));

     fh_replacekey(a, arr[2],-6);
     fh_replacekey(a, arr[6],-6);

     printf("Wert(minkey) %d\n",fh_minkey(a));
     printf("Knoten: %d\n\n", (int)fh_extractmin(a));

     printf("Wert(minkey) %d\n",fh_minkey(a));
     printf("Knoten: %d\n\n", (int)fh_extractmin(a));

     printf("Wert(minkey) %d\n",fh_minkey(a));
     printf("Knoten: %d\n\n", (int)fh_extractmin(a));
       	
     fh_deleteheap(a);

	return 0;
}



#define MAXELT          6000000

                   

#include "quicksort.c" 


int main()                           // overhead!
{
    int i=-1,j;
    double n;
    char t[10];
    void quicksort(double *,int *, int);
    
    double *  list = (double *) calloc(sizeof(double), MAXELT);
    int    *  dIndex = (int *) calloc(sizeof(int),MAXELT);
    
    for (i = 0; i < MAXELT;i++) {
       list[i] = rand(); 
        
    }    
  

   
        
    quicksort(list,dIndex,i); //Le paso el numero de elementos a ordenar
     
    
 
    

        
    printf("\nThe list obtained is ");
    for (j=0;j<i;j++)
        printf("\n %lf",list[j]);

    printf("\n\nProgram over.");
    
    return 0;       // successful termination.
}



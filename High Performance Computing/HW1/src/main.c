#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Function prototypes */
void fill_array(double*, size_t, const double, const double);
double simpson(double*, size_t, const double);


/* -------- Main Function -------- */
int main(int argc, char **argv){

    if (argc != 2){ // We need 1 argument
        fprintf(stderr, "\tYou have entered %d arguments. This executable requires exactly 1 argument.\n", argc-1);
        return EXIT_FAILURE;
    }

    // const int n_max = atoi(argv[1]);

    /* atoi is too limited...it will yield garbage 
       if the argument passed to the exectubale is too large. 
       The strto* family of functions is more reliable..

       I found the following on StackOverflow:
       "
          atoi is deprecated and thread-unsafe on many platforms;
          better ditch it in favour of strto(l|ll|ul|ull)
       "
       There's also a good explanation here:
       https://stackoverflow.com/questions/3420629/what-is-the-difference-between-sscanf-or-atoi-to-convert-a-string-to-an-integer
    */

    char *char_ptr;
    const long long n_max = strtoll(argv[1], &char_ptr, 10);

    /* The address of char_ptr must be passed to strtoll.
       What it does is grab a string of characters immediately
       following the entered integer, if the user makes the
       mistake of entering characters along with the integer.
        For instance, if I enter
            ./main 10garbage
        then char_ptr will store the string "garbage"
       You can test this by trying
            printf("String part is \"%s\"\n", char_ptr);
    */
    
    // printf("\n\tn_max = %lld\n", n_max);

    if (n_max <= 0 || n_max % 2 != 0){
        fprintf(stderr,"\tThe integer you enter must be positive and even! Let's try again ...\n");
        return EXIT_FAILURE;
    }


    // Allocate space for an array of doubles of length n_max + 1 
    double* f_array  = malloc((n_max+1) * sizeof *f_array);

    // make sure the memory has been successfully  allocated
    if (f_array == NULL) {
        printf("\tMemory not allocated. Perhaps your n_max value is too large...\n");
        exit(0);
    }

    // Integral parameters
    const double x_min = 0.;
    const double x_max = 1./sqrt(2);
    const double dx    = (x_max - x_min) / n_max;

    fill_array(f_array, n_max, x_min, dx);

    double num_sol   = simpson(f_array, n_max, dx);
    double exact_sol = .6426990816987241548;    // the exact solution = (2+Pi)/8
    double error     = fabs(num_sol - exact_sol);

    // printf("\n\t The numerical solution is = %.20g\n", num_sol);
    // printf("\n\t The exact solution is = %.20g\n", exact_sol);
    // printf("\n\t The absolute value of the error is = %.20f\n", error);

    printf("%.20f, %.20f, %.20f\n", dx, error, num_sol);

    // Free the memory
    free(f_array);

    return EXIT_SUCCESS;
}


/* Function implementations */
void fill_array(double* mem_ptr, size_t size, const double min_val, const double h){

    double x_grid [size+1];

    for (size_t i = 0; i <= size; ++i){
            x_grid[i]  = min_val + i * h;
            mem_ptr[i] = sqrt(1 - x_grid[i] * x_grid[i]);
    }
}

double simpson (double* mem_ptr, size_t size, const double h){
    
    /* Initially add just the endpoints */
    double S_n = h/3 * (mem_ptr[0] + mem_ptr[size]);  

    /* Now add the odd elements */
    for (size_t i = 1; i < size; i += 2)
        S_n += h/3 * 4. * mem_ptr[i];

    /* And lastly the even elements */
    for (size_t i = 2; i < size; i += 2)
        S_n += h/3 * 2. * mem_ptr[i];

    return S_n;  
}

#include <stdio.h>
#include <util.h>

/*
 * example driver for de Bruijn sequences.
 *
 * this code generates all n-mers over a protein
 * or dna alphabet. useful for creating fasta test sequences.
 */

/* k = 22 */
char proteinalphabet[] = "arndcqeghilkmfpstwyvbz";

/* k = 4 */
char dnaalphabet[] = "acgt";

int main(int argc, char *argv[])
{
  int i;
  int n, k;
  char *output;
  int outputsize;
  char *alphabet;

  if (argc != 4)
    {
      fprintf(stderr, "usage: %s n k prot-or-dna\n", argv[0]);
      fprintf(stderr, "where prot-or-dna = 0 for protein, 1 for dna\n");
      fprintf(stderr, "example: %s 3 22 0\n",argv[0]);
      fprintf(stderr, "example: %s 11 4 1\n",argv[0]);
      exit(1);
    }
  
  n = atoi(argv[1]);
  k = atoi(argv[2]);
  
  if (atoi(argv[3]) == 0)
    alphabet = proteinalphabet;
  else
    alphabet = dnaalphabet;

  /* output array needs:
   * k^n bytes - to store the de Bruijn sequence
   * n-1 bytes - to unwrap (see below)
   * 1   byte  - for the terminating NUL
   */

  outputsize = iexp(k,n) + (n-1) + 1;
  output = (char *) calloc( outputsize, sizeof(char));

  /* compute the (n,k) de Bruijn sequence */  
  debruijn(n,k,output,alphabet);
  
  /* but, we don't want a true cyclical de Bruijn sequence, we want
   * all words in a straight line. so we copy the first n-1 letters
   * to the end.
   */

  for(i=0;i<(n-1);i++)
    output[outputsize-((n-1)+1)+i] = output[i];

  /* don't forget to NUL-terminate it */

  output[outputsize-1] = '\0';

  fprintf(stderr,"n (word size) = %d\n",n);
  fprintf(stderr,"k (alphabet size) = %d\n",k);
  fprintf(stderr,"output size = k^n + (n-1) + 1 = %d + %d + 1 = %d\n",iexp(k,n),n-1,outputsize);
  fprintf(stderr,"naive size would have been %d\n",n * iexp(k,n));

  puts(output);

  free(output);
}

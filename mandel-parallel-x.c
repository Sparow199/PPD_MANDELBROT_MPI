#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>

#include "mpi.h"

/* Valeur par defaut des parametres */
#define N_ITER  255		/* nombre d'iterations */

#define X_MIN   -1.78		/* ensemble de Mandelbrot */
#define X_MAX   0.78
#define Y_MIN   -0.961
#define Y_MAX   0.961

#define X_SIZE  2048		/* dimension image */
#define Y_SIZE  1536
#define FILENAME "mandel.ppm"	/* image resultat */

typedef struct {
    int x_size, y_size;		/* dimensions */
    char *pixels;		/* matrice linearisee de pixels */
} picture_t;

static void
usage()
{
    fprintf(stderr, "usage : ./mandel [options]\n\n");
    fprintf(stderr, "Options \t Signification \t\t Val. defaut\n\n");
    fprintf(stderr, "-n \t\t Nbre iter. \t\t %d\n", N_ITER);
    fprintf(stderr, "-b \t\t Bornes \t\t %f %f %f %f\n",
	    X_MIN, X_MAX, Y_MIN, Y_MAX);
    fprintf(stderr, "-d \t\t Dimensions \t\t %d %d\n", X_SIZE, Y_SIZE);
    fprintf(stderr, "-f \t\t Fichier \t\t %s\n", FILENAME);

    exit(EXIT_FAILURE);
}

static void
parse_argv (int argc, char *argv[],
	    int *n_iter,
	    double *x_min, double *x_max, double *y_min, double *y_max,
	    int *x_size, int *y_size,
	    char **path)
{
    const char *opt = "b:d:n:f:";
    int c;

    /* Valeurs par defaut */
    *n_iter = N_ITER;
    *x_min  = X_MIN;
    *x_max  = X_MAX;
    *y_min  = Y_MIN;
    *y_max  = Y_MAX;
    *x_size = X_SIZE;
    *y_size = Y_SIZE;
    *path   = FILENAME;

    /* Analyse arguments */
    while ((c = getopt(argc, argv, opt)) != EOF) {
	switch (c) {
	    case 'b': 		/* domaine */
		sscanf(optarg, "%lf", x_min);
		sscanf(argv[optind++], "%lf", x_max);
		sscanf(argv[optind++], "%lf", y_min);
		sscanf(argv[optind++], "%lf", y_max);
		break;
	    case 'd':		/* largeur hauteur */
		sscanf(optarg, "%d", x_size);
		sscanf(argv[optind++], "%d", y_size);
		break;
	    case 'n':		/* nombre d'iterations */
		*n_iter = atoi(optarg);
		break;
	    case 'f':		/* fichier de sortie */
		*path = optarg;
		break;
	    default:
		usage();
	}
    }
}

static void
init_picture (picture_t *pict, int x_size, int y_size)
{
    pict->y_size = y_size;
    pict->x_size = x_size;
    pict->pixels = malloc(y_size * x_size); /* allocation espace memoire */
}

/* Enregistrement de l'image au format ASCII .ppm */
static void
save_picture (const picture_t *pict, const char *pathname)
{
    unsigned i;
    FILE *f = fopen(pathname, "w");

    fprintf(f, "P6\n%d %d\n255\n", pict->x_size, pict->y_size);
    for (i = 0 ; i < pict->x_size * pict->y_size; i++) {
	char c = pict->pixels[i];
	fprintf(f, "%c%c%c", c, c, c); /* monochrome blanc */
    }

    fclose (f);
}

static void
compute (picture_t *pict,
	 int nb_iter,
	 double x_min, double x_max, double y_min, double y_max)
{
    int pos = 0;
    int iy, ix, i;
    double pasx = (x_max - x_min) / pict->x_size, /* discretisation */
	   pasy = (y_max - y_min) / pict->y_size;

    /* Calcul en chaque point de l'image */
    for (iy = 0 ; iy < pict->y_size ; iy++) {
	for (ix = 0 ; ix < pict->x_size; ix++) {
	    double a = x_min + ix * pasx,
		b = y_max - iy * pasy,
		x = 0, y = 0;
	    for (i = 0 ; i < nb_iter ; i++) {
		double tmp = x;
		x = x * x - y * y + a;
		y = 2 * tmp * y + b;
		if (x * x + y * y > 4) /* divergence ! */
		    break;
	    }

	    pict->pixels[pos++] = (double) i / nb_iter * 255;
	}
    }
}

int
main (int argc, char *argv[])
{

    int n_iter,			/* degre de nettete  */
	x_size, y_size;		/* & dimensions de l'image */
    double x_min, x_max, y_min, y_max; /* bornes de la representation */
    char *pathname;		/* fichier destination */
    picture_t pict;   /* structure image finale*/
    picture_t seg_pict; /*structure image intermediaire sur chaque processus*/


/**************Initialisation de MPI****************************/

  int nb_process; /*Le nombre de processus*/
  int rank_process; /*Le rang du processus courant*/
  double segment_size; /*Taille du segment verticale de chaque processus*/


  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &nb_process);
  MPI_Comm_rank (MPI_COMM_WORLD, &rank_process);



/**************************************************************/


    parse_argv(argc, argv,
	       &n_iter,
	       &x_min, &x_max, &y_min, &y_max,
	       &x_size, &y_size, &pathname);

    segment_size = (x_max-x_min) / nb_process;

    init_picture ( &seg_pict, x_size / nb_process, y_size);

    if ( rank_process == 0 ){
        init_picture( &pict, x_size, y_size );
    }


    compute (&seg_pict, n_iter,x_min + segment_size * ( nb_process - rank_process - 1 ), x_min + segment_size * ( nb_process - rank_process ), y_min , y_max );

    MPI_Gather( seg_pict.pixels, seg_pict.x_size * seg_pict.y_size , MPI_CHAR,
               pict.pixels, seg_pict.x_size * seg_pict.y_size, MPI_CHAR,
               0, MPI_COMM_WORLD);

    if ( rank_process == 0 ) {
        save_picture (&pict, pathname);
    }

    MPI_Finalize();

    exit(EXIT_SUCCESS);
}

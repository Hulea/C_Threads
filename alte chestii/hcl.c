#include <pthread.h>
#include <stdio.h>
#include <semaphore.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>

int hydrogen_counter = 0;
int chlorine_counter = 0;
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
sem_t hydrogen, chlorine;

typedef struct
{
    int max_val;

} ATOM_STRUCT;

int i1 = 0;
int FLAG = 0;

pthread_mutex_t mutexxx = PTHREAD_MUTEX_INITIALIZER;
void hcl_bond(int maxi)
{
    pthread_mutex_lock(&mutexxx);

    if (i1 >= maxi)
    {
        FLAG = 1;
    }
    else
    {
        i1++;
        printf("**ch %d h %d Molecule 2xH,Cl no. %d created**\n\n",
               chlorine_counter, hydrogen_counter, i1);
    }

    pthread_mutex_unlock(&mutexxx);
}

void *hydrogen_func(void *sent_value)
{
    ATOM_STRUCT *aux = (ATOM_STRUCT *)sent_value;

    pthread_mutex_lock(&mutex1);

    if (FLAG == 1)
        pthread_exit(NULL);

    if (chlorine_counter >= 2 && hydrogen_counter >= 1)
    {
        hydrogen_counter++;
        printf("h+ : %d\n", hydrogen_counter);
        hydrogen_counter -= 2;
        printf("h-2 : %d\n", hydrogen_counter);

        sem_post(&hydrogen);
        sem_post(&hydrogen);

        chlorine_counter -= 2;
        printf("c-2 : %d\n", hydrogen_counter);

        sem_post(&chlorine);
        sem_post(&chlorine);

        pthread_mutex_unlock(&mutex1);

        hcl_bond(aux->max_val);

        pthread_exit(NULL);
    }

    else
    {
        hydrogen_counter++;
        printf("h+ : %d\n", hydrogen_counter);

        pthread_mutex_unlock(&mutex1);

        sem_wait(&hydrogen);
    }
}

void *chlorine_func(void *sent_value)
{
    ATOM_STRUCT *aux = (ATOM_STRUCT *)sent_value;

    pthread_mutex_lock(&mutex1);

    if (FLAG == 1)
        pthread_exit(NULL);

    if (hydrogen_counter >= 2 && chlorine_counter >= 1)

    {
        chlorine_counter++;
        printf("ch+ : %d\n", chlorine_counter);
        chlorine_counter -= 2;
        printf("ch-2 : %d\n", chlorine_counter);

        sem_post(&chlorine);
        sem_post(&chlorine);

        hydrogen_counter -= 2;
        printf("h-2 : %d\n", hydrogen_counter);

        sem_post(&hydrogen);
        sem_post(&hydrogen);

        pthread_mutex_unlock(&mutex1);
        hcl_bond(aux->max_val);
        pthread_exit(NULL);
    }

    else
    {
        chlorine_counter++;
        printf("ch+ : %d\n", chlorine_counter);
        pthread_mutex_unlock(&mutex1);

        sem_wait(&chlorine);
    }
}

int main(int argc, char **argv)
{

    int number_of_fecl2_molecules = atoi(argv[1]);

    if (sem_init(&hydrogen, 0, 0) == -1)
    {
        perror("error initilalizing h semaphore\n");
    }

    if (sem_init(&chlorine, 0, 0) == -1)
    {
        perror("error initilalizing ch semaphore\n");
    }

    srand(time(NULL));
    int h_nr = 2 * number_of_fecl2_molecules;
    int cl_nr = 2 * number_of_fecl2_molecules;

    printf("h nr : %d\n", h_nr);
    printf("cl nr : %d\n", cl_nr);

    ATOM_STRUCT atom_time_h[h_nr];
    ATOM_STRUCT atom_time_cl[cl_nr];

    for (int i = 0; i < h_nr; i++)
    {
        atom_time_h[i].max_val = number_of_fecl2_molecules;
    }
    for (int i = 0; i < cl_nr; i++)
    {
        atom_time_cl[i].max_val = number_of_fecl2_molecules;
    }

    pthread_t h_th[h_nr], cl_th[cl_nr];

    for (int i = 0; i < h_nr; i++)
    {
        pthread_create(&h_th[i], NULL, hydrogen_func, (void *)&atom_time_h[i]);
    }
    for (int i = 0; i < cl_nr; i++)
    {
        pthread_create(&cl_th[i], NULL, chlorine_func, (void *)&atom_time_cl[i]);
    }

    for (int i = 0; i < h_nr; i++)
    {
        pthread_join(h_th[i], NULL);
    }
    for (int i = 0; i < cl_nr; i++)
    {
        pthread_join(cl_th[i], NULL);
    }

    sem_destroy(&hydrogen);
    sem_destroy(&chlorine);
}
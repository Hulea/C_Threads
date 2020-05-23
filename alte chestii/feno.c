#include <pthread.h>
#include <stdio.h>
#include <semaphore.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>

int iron_counter = 0;
int nitrogen_counter = 0;
int oxygen_counter = 0;
int hydrochloric_acid_counter = 0;
int iron_hydroxide_counter = 0;
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
sem_t iron, nitrogen, oxygen, first_molecule, second_molecule, final_molecule;

typedef struct
{
    int max_val;

} ATOM_STRUCT;

int i1 = 0;
int FLAG = 0;

pthread_mutex_t feno_mutex = PTHREAD_MUTEX_INITIALIZER;
void feno_bond(int maxi)
{
    pthread_mutex_lock(&feno_mutex);

    if (i1 >= maxi)
    {
        FLAG = 1;
    }
    else
    {
        i1++;
        printf("**fe %d n %d o %d Molecule Fe(NO3)2 no. %d created**\n\n",
               iron_counter, nitrogen_counter, oxygen_counter, i1);
    }

    pthread_mutex_unlock(&feno_mutex);
}

void *iron_func(void *sent_value)
{

    ATOM_STRUCT *aux = (ATOM_STRUCT *)sent_value;

    pthread_mutex_lock(&mutex2);

    if (FLAG == 1)
        pthread_exit(NULL);
    iron_counter++;
    printf("fe+: %d\n", iron_counter);
    if (nitrogen_counter >= 2 && oxygen_counter >= 6 && iron_counter >= 1)
    {

        iron_counter--;
        printf("fe-: %d\n", iron_counter);
        sem_post(&iron);

        nitrogen_counter -= 2;
        printf("n-2: %d\n", nitrogen_counter);
        sem_post(&nitrogen);
        sem_post(&nitrogen);

        oxygen_counter -= 6;
        printf("o-6: %d\n", oxygen_counter);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);

        pthread_mutex_unlock(&mutex2);
        feno_bond(aux->max_val);
    }

    else
    {
        pthread_mutex_unlock(&mutex2);
        sem_wait(&iron);
    }
}

void *nitrogen_func(void *sent_value)
{

    ATOM_STRUCT *aux = (ATOM_STRUCT *)sent_value;

    pthread_mutex_lock(&mutex2);

    if (FLAG == 1)
        pthread_exit(NULL);

    nitrogen_counter++;
    printf("n+: %d\n", nitrogen_counter);
    if (iron_counter >= 1 && oxygen_counter >= 6 && nitrogen_counter >= 2)
    {

        nitrogen_counter -= 2;
        printf("n-2: %d\n", nitrogen_counter);
        sem_post(&nitrogen);
        sem_post(&nitrogen);

        iron_counter--;
        printf("fe-: %d\n", iron_counter);
        sem_post(&iron);

        oxygen_counter -= 6;
        printf("o-6: %d\n", oxygen_counter);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);

        pthread_mutex_unlock(&mutex2);
        feno_bond(aux->max_val);
    }
    else
    {

        pthread_mutex_unlock(&mutex2);
        sem_wait(&nitrogen);
    }
}

void *oxygen_func(void *sent_value)
{
    ATOM_STRUCT *aux = (ATOM_STRUCT *)sent_value;

    pthread_mutex_lock(&mutex2);

    if (FLAG == 1)
        pthread_exit(NULL);
    oxygen_counter++;
    printf("o+: %d\n", oxygen_counter);
    if (iron_counter >= 1 && nitrogen_counter >= 2 && oxygen_counter >= 6)
    {

        oxygen_counter -= 6;
        printf("o-6: %d\n", oxygen_counter);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);

        iron_counter--;
        printf("fe-: %d\n", iron_counter);
        sem_post(&iron);

        nitrogen_counter -= 2;
        printf("n-2: %d\n", nitrogen_counter);
        sem_post(&nitrogen);
        sem_post(&nitrogen);

        pthread_mutex_unlock(&mutex2);
        feno_bond(aux->max_val);
    }
    else
    {

        pthread_mutex_unlock(&mutex2);
        sem_wait(&oxygen);
    }
}

int main(int argc, char **argv)
{

    int number_of_fecl2_molecules = atoi(argv[1]);

    if (sem_init(&iron, 0, 0) == -1)
    {
        perror("error initilalizing fe semaphore\n");
    }

    if (sem_init(&nitrogen, 0, 0) == -1)
    {
        perror("error initilalizing semaphore\n");
    }

    if (sem_init(&oxygen, 0, 0) == -1)
    {
        perror("error initilalizing o semaphore\n");
    }

    if (sem_init(&first_molecule, 0, 0) == -1)
    {
        perror("error initilalizing 3hcl semaphore\n");
    }

    if (sem_init(&second_molecule, 0, 0) == -1)
    {
        perror("error initilalizing fe(oh)3 semaphore\n");
    }

    if (sem_init(&final_molecule, 0, 0) == -1)
    {
        perror("error initilalizing + semaphore\n");
    }

    srand(time(NULL));

    int fe_nr = 1 * number_of_fecl2_molecules ;
    int n_nr = 2 * number_of_fecl2_molecules ;
    int o_nr = 6 * number_of_fecl2_molecules ;

    printf("fe nr : %d\n", fe_nr);
    printf("n nr : %d\n", n_nr);
    printf("o nr : %d\n\n", o_nr);

    ATOM_STRUCT atom_time_fe[fe_nr];
    ATOM_STRUCT atom_time_n[n_nr];
    ATOM_STRUCT atom_time_o[o_nr];

    for (int i = 0; i < fe_nr; i++)
    {
        atom_time_fe[i].max_val = number_of_fecl2_molecules;
    }
    for (int i = 0; i < n_nr; i++)
    {
        atom_time_n[i].max_val = number_of_fecl2_molecules;
    }
    for (int i = 0; i < o_nr; i++)
    {
        atom_time_o[i].max_val = number_of_fecl2_molecules;
    }
    /////////////////////

    pthread_t fe_th[fe_nr], n_th[n_nr], o_th[o_nr];

    for (int i = 0; i < fe_nr; i++)
    {
        pthread_create(&fe_th[i], NULL, iron_func, (void *)&atom_time_fe[i]);
    }
    for (int i = 0; i < n_nr; i++)
    {
        pthread_create(&n_th[i], NULL, nitrogen_func, (void *)&atom_time_n[i]);
    }
    for (int i = 0; i < o_nr; i++)
    {
        pthread_create(&o_th[i], NULL, oxygen_func, (void *)&atom_time_o[i]);
    }

    for (int i = 0; i < fe_nr; i++)
    {
        pthread_join(fe_th[i], NULL);
    }
    for (int i = 0; i < n_nr; i++)
    {
        pthread_join(n_th[i], NULL);
    }
    for (int i = 0; i < o_nr; i++)
    {
        pthread_join(o_th[i], NULL);
    }

    sem_destroy(&iron);
    sem_destroy(&nitrogen);
    sem_destroy(&oxygen);
}
//Hulea Andrei-Florin
//Grupa 30225
//Specificatie individuala : 2HCl + Fe(NO3)2 → FeCl2 + 2HNO3

#include <pthread.h>
#include <stdio.h>
#include <semaphore.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <string.h>

int fd;
int fs;
struct timespec start;
#define bill 1000000000L;
double nanoseconds()
{
    struct timespec stop;

    double s;

    if (clock_gettime(CLOCK_REALTIME, &stop) == -1)
    {
        printf("Error generating time\n");
        exit(0);
    }

    s = (stop.tv_sec - start.tv_sec) + (double)(stop.tv_nsec - start.tv_nsec) / (double)bill;
    return s;
}

struct Queue
{
    int head, tail, size;
    int capacity;
    int *elements;
};

struct Queue *create(int capacity)
{
    struct Queue *queue = (struct Queue *)malloc(sizeof(struct Queue));
    queue->capacity = capacity;
    queue->head = queue->size = 0;
    queue->tail = capacity - 1; // This is important, see the enqueue
    queue->elements = (int *)malloc(queue->capacity * sizeof(int));
    return queue;
}

int full_q(struct Queue *queue)
{
    return (queue->size == queue->capacity);
}

int empty_q(struct Queue *queue)
{
    return (queue->size == 0);
}

int returnTail(struct Queue *queue)
{
    if (empty_q(queue))
        return INT_MIN;
    return queue->elements[queue->tail];
}

void enqueue(struct Queue *queue, int item)
{
    if (full_q(queue))
        return;
    queue->tail = (queue->tail + 1) % queue->capacity;
    queue->elements[queue->tail] = item;
    queue->size = queue->size + 1;
}

int dequeue(struct Queue *queue)
{
    if (empty_q(queue))
        return INT_MIN;
    int item = queue->elements[queue->head];
    queue->head = (queue->head + 1) % queue->capacity;
    queue->size = queue->size - 1;
    return item;
}

struct Queue *h_q;
struct Queue *cl_q;
struct Queue *fe_q;
struct Queue *n_q;
struct Queue *o_q;

struct Queue *hcl_q;
struct Queue *feno_q;

char output[1000000];

int hydrogen_counter = 0;
int chlorine_counter = 0;
int iron_counter = 0;
int nitrogen_counter = 0;
int oxygen_counter = 0;
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
sem_t hydrogen, chlorine, iron, nitrogen, oxygen;

typedef struct
{
    int max_val;
    int molecule_id;
    char *molecule_type;
    double creation_time;
    double end_time;
    long time_in_nanoseconds;
} MOLECULE_STRUCT;

typedef struct
{
    int max_val;
    int atom_id;
    char *atom_type;
    double creation_time;
    double end_time;
    long time_in_nanoseconds;
} ATOM_STRUCT;

typedef struct
{
    int nr_reactie;
    double reaction_time;
    char *reaction_type;
    long time_in_nanoseconds;
} REACTION_STRUCT;

int big_equation_counter = 0;
int hcl_counter = 0;
int feno_counter = 0;
pthread_mutex_t final_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t final_th_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;
void *final_func(void *sent_value)
{
    long aux = (long)sent_value;

    pthread_mutex_lock(&final_th_mutex);

    if (big_equation_counter >= (int)aux)
        printf("Stopping\n");
    else
    {

        big_equation_counter++;
        REACTION_STRUCT reac;
        reac.nr_reactie = big_equation_counter;
        reac.reaction_type = "2HCl + Fe(NO3)2 -> FeCl2 + 2HNO3";
        reac.reaction_time = nanoseconds();

        fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777);
        write(fs, &reac, sizeof(reac));

        fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
        char toWrite[150];
        sprintf(toWrite, "\nREACTIE - TipReactie: %s, NrReactie: %d, Moment de timp: %lf\n",
                reac.reaction_type, reac.nr_reactie, reac.reaction_time);
        write(fd, toWrite, strlen(toWrite));
        toWrite[0] = '\0';

        int hcl1 = dequeue(hcl_q);
        int feno1 = dequeue(feno_q);

        sprintf(toWrite, "->TipMolecula: %s, NrMolecula: %d\n->TipMolecula: %s, NrMolecula: %d\n",
                "HCl", hcl1, "Fe(NO3)2", feno1);
        write(fd, toWrite, strlen(toWrite));

        toWrite[0] = '\0';
        sprintf(toWrite, "TERMINARE_MOLECULA - TipMolecula %s, NrMolecula: %d, Moment de timp: %lf\n",
                "2HCl", hcl1, nanoseconds());
        write(fd, toWrite, strlen(toWrite));
        toWrite[0] = '\0';
        sprintf(toWrite, "TERMINARE_MOLECULA - TipMolecula %s, NrMolecula: %d, Moment de timp: %lf\n",
                "Fe(NO3)2", feno1, nanoseconds());
        write(fd, toWrite, strlen(toWrite));
    }
    pthread_mutex_unlock(&final_th_mutex);
}

void final_bond(int maxi)
{
    pthread_mutex_lock(&final_mutex);

    pthread_t final_th;
    pthread_create(&final_th, NULL, final_func, (void *)(long)maxi);
    pthread_join(final_th, NULL);

    pthread_mutex_unlock(&final_mutex);
}

int enough_hcl = 0, enough_feno = 0;

pthread_mutex_t hcl_mutex = PTHREAD_MUTEX_INITIALIZER;
int no_of_hcl = 0;
pthread_mutex_t hcl_th_mutex = PTHREAD_MUTEX_INITIALIZER;

void *hcl_func(void *sent_value)
{
    long aux = (long)sent_value;

    pthread_mutex_lock(&hcl_th_mutex);
    hcl_counter++;
    no_of_hcl++;

    MOLECULE_STRUCT molec;
    molec.creation_time = nanoseconds();
    molec.max_val = (int)aux;
    molec.molecule_id = no_of_hcl;
    molec.molecule_type = "2HCl";

    fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777);
    write(fs, &molec, sizeof(molec));

    fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
    char toWrite[150];
    sprintf(toWrite, "\nFORMARE - TipMolecula: %s, NrMolecula: %d, Moment de timp: %lf\n",
            molec.molecule_type, molec.molecule_id, molec.creation_time);

    write(fd, toWrite, strlen(toWrite));

    enqueue(hcl_q, molec.molecule_id);

    int h1 = dequeue(h_q);
    int h2 = dequeue(h_q);
    int cl1 = dequeue(cl_q);
    int cl2 = dequeue(cl_q);

    toWrite[0] = '\0';
    sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "hydrogen", h1, "hydrogen", h2);
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "chlorine", cl1, "chlorine", cl2);
    write(fd, toWrite, strlen(toWrite));

    if (hcl_counter > 0 && feno_counter > 0)
    {
        hcl_counter--;
        feno_counter--;
        final_bond((int)aux);
    }

    toWrite[0] = '\0';
    sprintf(toWrite, "\nTERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "hydrogen", h1, nanoseconds());

    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "hydrogen", h2, nanoseconds());

    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "chlorine", cl1, nanoseconds());

    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "chlorine", cl2, nanoseconds());

    write(fd, toWrite, strlen(toWrite));

    pthread_mutex_unlock(&hcl_th_mutex);
}

void hcl_bond(ATOM_STRUCT *maxi)
{
    pthread_mutex_lock(&hcl_mutex);
    {
        pthread_t hcl_th;
        pthread_create(&hcl_th, NULL, hcl_func, (void *)(long)maxi->max_val);
        pthread_join(hcl_th, NULL);
    }
    pthread_mutex_unlock(&hcl_mutex);
}

pthread_mutex_t feno_mutex = PTHREAD_MUTEX_INITIALIZER;
int no_of_feno = 0;
pthread_mutex_t feno_th_mutex = PTHREAD_MUTEX_INITIALIZER;

void *feno_func(void *sent_value)
{
    long aux = (long)sent_value;

    pthread_mutex_lock(&feno_th_mutex);

    feno_counter++;
    no_of_feno++;

    MOLECULE_STRUCT molec;
    molec.creation_time = nanoseconds();
    molec.max_val = (int)aux;
    molec.molecule_id = no_of_feno;
    molec.molecule_type = "Fe(NO3)2";

    fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777);
    write(fs, &molec, sizeof(molec));

    fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
    char toWrite[150];
    sprintf(toWrite, "\nFORMARE - TipMolecula: %s, NrMolecula: %d, Moment de timp: %lf\n",
            molec.molecule_type, molec.molecule_id, molec.creation_time);

    write(fd, toWrite, strlen(toWrite));

    enqueue(feno_q, molec.molecule_id);

    int fe1 = dequeue(fe_q);
    int n1 = dequeue(n_q);
    int n2 = dequeue(n_q);
    int o1 = dequeue(o_q);
    int o2 = dequeue(o_q);
    int o3 = dequeue(o_q);
    int o4 = dequeue(o_q);
    int o5 = dequeue(o_q);
    int o6 = dequeue(o_q);

    toWrite[0] = '\0';
    sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n", "iron", fe1);
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "nitrogen", n1, "nitrogen", n2);
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "oxygen", o1, "oxygen", o2);
    write(fd, toWrite, strlen(toWrite));
    toWrite[0] = '\0';
    sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "oxygen", o3, "oxygen", o4);
    write(fd, toWrite, strlen(toWrite));
    toWrite[0] = '\0';
    sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "oxygen", o5, "oxygen", o6);
    write(fd, toWrite, strlen(toWrite));

    if (hcl_counter > 0 && feno_counter > 0)
    {
        feno_counter--;
        hcl_counter--;
        final_bond(aux);
    }

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "iron", fe1, nanoseconds());
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "nitrogen", n1, nanoseconds());
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "nitrogen", n2, nanoseconds());
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "oxygen", o1, nanoseconds());
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "oxygen", o2, nanoseconds());
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "oxygen", o3, nanoseconds());
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "oxygen", o4, nanoseconds());
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "oxygen", o5, nanoseconds());
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
            "oxygen", o6, nanoseconds());
    write(fd, toWrite, strlen(toWrite));

    pthread_mutex_unlock(&feno_th_mutex);
}

void feno_bond(ATOM_STRUCT *maxi)
{
    pthread_mutex_lock(&feno_mutex);
    {
        pthread_t feno_th;
        pthread_create(&feno_th, NULL, feno_func, (void *)(long)maxi->max_val);
        pthread_join(feno_th, NULL);
    }
    pthread_mutex_unlock(&feno_mutex);
}

void *hydrogen_func(void *sent_value)
{
    ATOM_STRUCT *aux = (ATOM_STRUCT *)sent_value;

    pthread_mutex_lock(&mutex1);

    if (chlorine_counter >= 2 && hydrogen_counter >= 1)
    {

        hydrogen_counter++;
        enqueue(h_q, aux->atom_id);

        hydrogen_counter -= 2;

        sem_post(&hydrogen);
        sem_post(&hydrogen);

        chlorine_counter -= 2;

        sem_post(&chlorine);
        sem_post(&chlorine);

        pthread_mutex_unlock(&mutex1);

        hcl_bond(aux);
    }

    else
    {
        hydrogen_counter++;

        enqueue(h_q, aux->atom_id);

        pthread_mutex_unlock(&mutex1);

        sem_wait(&hydrogen);
    }
}

void *chlorine_func(void *sent_value)
{
    ATOM_STRUCT *aux = (ATOM_STRUCT *)sent_value;

    pthread_mutex_lock(&mutex1);

    if (hydrogen_counter >= 2 && chlorine_counter >= 1)

    {
        chlorine_counter++;
        enqueue(cl_q, aux->atom_id);

        chlorine_counter -= 2;

        sem_post(&chlorine);
        sem_post(&chlorine);

        hydrogen_counter -= 2;

        sem_post(&hydrogen);
        sem_post(&hydrogen);

        pthread_mutex_unlock(&mutex1);

        hcl_bond(aux);
    }

    else
    {
        chlorine_counter++;

        enqueue(cl_q, aux->atom_id);

        pthread_mutex_unlock(&mutex1);

        sem_wait(&chlorine);
    }
}

void *iron_func(void *sent_value)
{

    ATOM_STRUCT *aux = (ATOM_STRUCT *)sent_value;

    pthread_mutex_lock(&mutex2);

    if (nitrogen_counter >= 2 && oxygen_counter >= 6)
    {
        iron_counter++;
        enqueue(fe_q, aux->atom_id);

        iron_counter--;
        nitrogen_counter -= 2;
        sem_post(&nitrogen);
        sem_post(&nitrogen);

        oxygen_counter -= 6;
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);

        pthread_mutex_unlock(&mutex2);

        feno_bond(aux);
    }

    else
    {
        iron_counter++;

        enqueue(fe_q, aux->atom_id);

        pthread_mutex_unlock(&mutex2);
        sem_wait(&iron);
    }
}

void *nitrogen_func(void *sent_value)
{

    ATOM_STRUCT *aux = (ATOM_STRUCT *)sent_value;

    pthread_mutex_lock(&mutex2);

    if (iron_counter >= 1 && oxygen_counter >= 6 && nitrogen_counter >= 1)
    {
        nitrogen_counter++;
        enqueue(n_q, aux->atom_id);

        nitrogen_counter -= 2;
        sem_post(&nitrogen);
        sem_post(&nitrogen);

        iron_counter--;
        sem_post(&iron);

        oxygen_counter -= 6;
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);

        pthread_mutex_unlock(&mutex2);

        feno_bond(aux);
    }
    else
    {
        nitrogen_counter++;

        enqueue(n_q, aux->atom_id);

        pthread_mutex_unlock(&mutex2);
        sem_wait(&nitrogen);
    }
}

void *oxygen_func(void *sent_value)
{
    ATOM_STRUCT *aux = (ATOM_STRUCT *)sent_value;

    pthread_mutex_lock(&mutex2);

    if (iron_counter >= 1 && nitrogen_counter >= 2 && oxygen_counter >= 5)
    {

        oxygen_counter++;

        enqueue(o_q, aux->atom_id);

        oxygen_counter -= 6;
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);
        sem_post(&oxygen);

        iron_counter--;

        sem_post(&iron);

        nitrogen_counter -= 2;

        sem_post(&nitrogen);
        sem_post(&nitrogen);

        pthread_mutex_unlock(&mutex2);

        feno_bond(aux);
    }
    else
    {
        oxygen_counter++;
        //printf("o+: %d\n", oxygen_counter);

        enqueue(o_q, aux->atom_id);

        pthread_mutex_unlock(&mutex2);
        sem_wait(&oxygen);
    }
}

int main(int argc, char **argv)
{
         int number_of_fecl2_molecules;
    if (argc == 2)
    {
        number_of_fecl2_molecules= atoi(argv[1]);
        system(argv[2]);
    }
    else if (argc == 3)
    {
        number_of_fecl2_molecules = atoi(argv[2]);
        system(argv[3]);
    }
    else
    {
        printf("Numar gresit de argumente\n");
        exit(0);
    }

    if (clock_gettime(CLOCK_REALTIME, &start) == -1)
    {
        perror("clock gettime");
        return EXIT_FAILURE;
    }

    h_q = create(10000);
    cl_q = create(10000);
    fe_q = create(10000);
    n_q = create(10000);
    o_q = create(10000);

    hcl_q = create(10000);
    feno_q = create(10000);

    if (sem_init(&hydrogen, 0, 0) == -1)
    {
        perror("error initilalizing h semaphore\n");
    }

    if (sem_init(&chlorine, 0, 0) == -1)
    {
        perror("error initilalizing ch semaphore\n");
    }

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

    int h_nr = 2 * number_of_fecl2_molecules;
    int cl_nr = 2 * number_of_fecl2_molecules;
    int fe_nr = 1 * number_of_fecl2_molecules;
    int n_nr = 2 * number_of_fecl2_molecules;
    int o_nr = 6 * number_of_fecl2_molecules;

    ATOM_STRUCT atom_time_h[h_nr];
    ATOM_STRUCT atom_time_cl[cl_nr];
    ATOM_STRUCT atom_time_fe[fe_nr];
    ATOM_STRUCT atom_time_n[n_nr];
    ATOM_STRUCT atom_time_o[o_nr];

    pthread_t h_th[h_nr], cl_th[cl_nr], fe_th[fe_nr], n_th[n_nr], o_th[o_nr];

    int ok = 0;

    for (int i = 0; i < h_nr; i++)
    {

        atom_time_h[i].max_val = number_of_fecl2_molecules;
        atom_time_h[i].atom_id = i;
        atom_time_h[i].atom_type = "hydrogen";
        atom_time_h[i].creation_time = nanoseconds();

        pthread_create(&h_th[i], NULL, hydrogen_func, (void *)&atom_time_h[i]);

        if (ok == 0)
            fd = open("log.txt", O_CREAT | O_TRUNC | O_RDWR, 0777);
        else
            fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
        ok = 1;
        char aux[150];
        sprintf(aux, "CREARE - TipAtom: %s, IdAtom: %d,  Moment de timp: %lf \n",
                atom_time_h[i].atom_type, atom_time_h[i].atom_id, atom_time_h[i].creation_time);

        write(fd, aux, strlen(aux));
    }

    fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777);
    write(fs, &atom_time_h, sizeof(atom_time_h));

    for (int i = 0; i < cl_nr; i++)
    {

        atom_time_cl[i].max_val = number_of_fecl2_molecules;
        atom_time_cl[i].atom_id = i;
        atom_time_cl[i].atom_type = "chlorine";
        atom_time_cl[i].creation_time = nanoseconds();

        pthread_create(&cl_th[i], NULL, chlorine_func, (void *)&atom_time_cl[i]);

        fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
        char aux[150];
        sprintf(aux, "CREARE - TipAtom: %s, IdAtom: %d, Moment de timp: %lf \n",
                atom_time_cl[i].atom_type, atom_time_cl[i].atom_id, atom_time_cl[i].creation_time);

        write(fd, aux, strlen(aux));
    }

    fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777);
    write(fs, &atom_time_cl, sizeof(atom_time_cl));

    for (int i = 0; i < fe_nr; i++)
    {

        atom_time_fe[i].max_val = number_of_fecl2_molecules;
        atom_time_fe[i].atom_id = i;
        atom_time_fe[i].atom_type = "iron";
        atom_time_fe[i].creation_time = nanoseconds();
        pthread_create(&fe_th[i], NULL, iron_func, (void *)&atom_time_fe[i]);

        fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
        char aux[150];
        sprintf(aux, "CREARE - TipAtom: %s, IdAtom: %d, Moment de timp: %lf \n",
                atom_time_fe[i].atom_type, atom_time_fe[i].atom_id, atom_time_fe[i].creation_time);

        write(fd, aux, strlen(aux));
    }

    fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777);
    write(fs, &atom_time_fe, sizeof(atom_time_fe));

    for (int i = 0; i < n_nr; i++)
    {

        atom_time_n[i].max_val = number_of_fecl2_molecules;
        atom_time_n[i].atom_id = i;
        atom_time_n[i].atom_type = "nitrogen";
        atom_time_n[i].creation_time = nanoseconds();
        pthread_create(&n_th[i], NULL, nitrogen_func, (void *)&atom_time_n[i]);

        fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
        char aux[150];
        sprintf(aux, "CREARE - TipAtom: %s, IdAtom: %d, Moment de timp: %lf \n",
                atom_time_n[i].atom_type, atom_time_n[i].atom_id, atom_time_n[i].creation_time);

        write(fd, aux, strlen(aux));
    }

    fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777);
    write(fs, &atom_time_n, sizeof(atom_time_n));

    for (int i = 0; i < o_nr; i++)
    {

        atom_time_o[i].max_val = number_of_fecl2_molecules;
        atom_time_o[i].atom_id = i;
        atom_time_o[i].atom_type = "oxygen";
        atom_time_o[i].creation_time = nanoseconds();
        pthread_create(&o_th[i], NULL, oxygen_func, (void *)&atom_time_o[i]);

        fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
        char aux[150];
        sprintf(aux, "CREARE - TipAtom: %s, IdAtom: %d, Moment de timp: %lf \n",
                atom_time_o[i].atom_type, atom_time_o[i].atom_id, atom_time_o[i].creation_time);

        write(fd, aux, strlen(aux));
    }

    fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777);
    write(fs, &atom_time_o, sizeof(atom_time_o));

    for (int i = 0; i < h_nr; i++)
        pthread_join(h_th[i], NULL);

    for (int i = 0; i < cl_nr; i++)
        pthread_join(cl_th[i], NULL);

    for (int i = 0; i < fe_nr; i++)
        pthread_join(fe_th[i], NULL);

    for (int i = 0; i < n_nr; i++)
        pthread_join(n_th[i], NULL);

    for (int i = 0; i < o_nr; i++)
        pthread_join(o_th[i], NULL);

    if (big_equation_counter == number_of_fecl2_molecules)
    {
        char aux[150];
        sprintf(aux, "\nCONDITIE_DE_FINAL - %lf\n", nanoseconds());
        fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
        write(fd, aux, strlen(aux));
    }

    sem_destroy(&hydrogen);
    sem_destroy(&chlorine);
    sem_destroy(&iron);
    sem_destroy(&nitrogen);
    sem_destroy(&oxygen);

    close(fd);
}
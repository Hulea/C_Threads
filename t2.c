/*Hulea Andrei-Florin
Grupa 30225
Specificatie individuala : 2HCl + Fe(NO3)2 → FeCl2 + 2HNO3

Mod compilare: gcc t2.c -o t2 -lpthread
Mod rulare:  ./t2 N sau ./t2 -o N
N reprezinta numarul de molecule de FeCl2 ce trebuie generate.

In urma rularii vor fi create 2 fisiere: log.dat si log.txt
In log.dat sunt stocate structurile cu informatiile necesare pentru atomi,molecule si reactii.
In log.txt sunt stocate informatiile sub forma de string.
*/

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

//Genereaza timpul in microsecunde de la inceperea programului pana la apelarea functiei
double microseconds()
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

//Structura pentru coada si functiile pentru aceasta
struct Queue
{
    int head, tail, size;
    int qsize;
    int *elements;
};

struct Queue *create(int qsize)
{
    struct Queue *queue = (struct Queue *)malloc(sizeof(struct Queue));
    queue->qsize = qsize;
    queue->head = queue->size = 0;
    queue->tail = qsize - 1;
    queue->elements = (int *)malloc(queue->qsize * sizeof(int));
    return queue;
}

int full_q(struct Queue *queue)
{
    return (queue->size == queue->qsize);
}

int empty_q(struct Queue *queue)
{
    return (queue->size == 0);
}

void enqueue(struct Queue *queue, int item)
{
    if (full_q(queue))
        return;
    queue->tail = (queue->tail + 1) % queue->qsize;
    queue->elements[queue->tail] = item;
    queue->size = queue->size + 1;
}

int dequeue(struct Queue *queue)
{
    if (empty_q(queue))
        return INT_MIN;
    int item = queue->elements[queue->head];
    queue->head = (queue->head + 1) % queue->qsize;
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

int hydrogen_counter = 0;
int chlorine_counter = 0;
int iron_counter = 0;
int nitrogen_counter = 0;
int oxygen_counter = 0;
pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
sem_t hydrogen, chlorine, iron, nitrogen, oxygen;

//Structuri pentru atom,molecula si reactie
typedef struct
{
    int max_val;
    int molecule_id;
    char *molecule_type;
    double creation_time;
    double end_time;
    long time_in_microseconds;
} MOLECULE_STRUCT;

typedef struct
{
    int max_val;
    int atom_id;
    char *atom_type;
    double creation_time;
    double end_time;
    long time_in_microseconds;
} ATOM_STRUCT;

typedef struct
{
    int nr_reactie;
    double reaction_time;
    char *reaction_type;
    long time_in_microseconds;
} REACTION_STRUCT;

int big_equation_counter = 0;
int hcl_counter = 0;
int feno_counter = 0;
pthread_mutex_t final_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t final_th_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;

//Creeaza reactia finala si afiseaza moleculele din compozitia acesteia
//In caz ca se formeaza mai multe reactii decat trebuie este afisat mesajul 'Too many reactions!!!'
//Numarul de reactii este echivalent cu N-ul pentru conditia de terminare,
//adica pentru a genera N molecule ce FeCl2 avem nevoie de N*(2Hcl + Fe(NO3)2)
void *reaction_func(void *sent_value)
{
    long aux = (long)sent_value;

    pthread_mutex_lock(&final_th_mutex);

    if (big_equation_counter >= (int)aux)
        printf("Too many reactions!!!\n");
    else
    {

        big_equation_counter++;
        REACTION_STRUCT reac;
        reac.nr_reactie = big_equation_counter;
        reac.reaction_type = "2HCl + Fe(NO3)2 -> FeCl2 + 2HNO3";
        reac.reaction_time = microseconds();

        if ((fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
        {
            printf("Error opening log.dat\n");
            exit(0);
        }
        write(fs, &reac, sizeof(reac));

        if ((fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
        {
            printf("Error opening log.txt\n");
            exit(0);
        }
        char toWrite[150];
        if (sprintf(toWrite, "\nREACTIE - TipReactie: %s, NrReactie: %d, Moment de timp: %lf\n",
                    reac.reaction_type, reac.nr_reactie, reac.reaction_time) < 0)
        {
            printf("Error at sprintf\n");
            exit(0);
        }
        write(fd, toWrite, strlen(toWrite));
        toWrite[0] = '\0';

        int hcl1 = dequeue(hcl_q);
        int feno1 = dequeue(feno_q);

        if (sprintf(toWrite, "->TipMolecula: %s, NrMolecula: %d\n->TipMolecula: %s, NrMolecula: %d\n",
                    "HCl", hcl1, "Fe(NO3)2", feno1) < 0)
        {
            printf("Error at sprintf\n");
            exit(0);
        }
        write(fd, toWrite, strlen(toWrite));

        MOLECULE_STRUCT molec;
        molec.molecule_type = "2HCl";
        molec.molecule_id = hcl1;
        molec.end_time = microseconds();

        toWrite[0] = '\0';
        if (sprintf(toWrite, "TERMINARE_MOLECULA - TipMolecula %s, NrMolecula: %d, Moment de timp: %lf\n",
                    molec.molecule_type, molec.molecule_id, molec.end_time) < 0)
        {
            printf("Error at sprintf\n");
            exit(0);
        }
        write(fd, toWrite, strlen(toWrite));
        write(fs, &molec, sizeof(molec));

        molec.molecule_type = "Fe(NO3)2";
        molec.molecule_id = feno1;
        molec.end_time = microseconds();
        toWrite[0] = '\0';
        if (sprintf(toWrite, "TERMINARE_MOLECULA - TipMolecula %s, NrMolecula: %d, Moment de timp: %lf\n",
                    molec.molecule_type, molec.molecule_id, molec.end_time) < 0)
        {
            printf("Error at sprintf\n");
            exit(0);
        }
        write(fd, toWrite, strlen(toWrite));
        write(fs, &molec, sizeof(molec));
    }
    pthread_mutex_unlock(&final_th_mutex);
}

//Creeaza threadul pentru reactie
void reaction(int maxi)
{
    pthread_mutex_lock(&final_mutex);

    pthread_t final_th;
    if (pthread_create(&final_th, NULL, reaction_func, (void *)(long)maxi))
    {
        printf("Error creating reaction thread\n");
        exit(0);
    }
    if (pthread_join(final_th, NULL))
    {
        printf("Error joining reaction thread\n");
        exit(0);
    }

    pthread_mutex_unlock(&final_mutex);
}

int enough_hcl = 0, enough_feno = 0;
pthread_mutex_t hcl_mutex = PTHREAD_MUTEX_INITIALIZER;
int no_of_hcl = 0;
pthread_mutex_t hcl_th_mutex = PTHREAD_MUTEX_INITIALIZER;

//Creeaza o molecula de 2HCl si afiseaza atomii din compozitia acesteia
void *hcl_func(void *sent_value)
{
    long aux = (long)sent_value;

    pthread_mutex_lock(&hcl_th_mutex);
    hcl_counter++;
    no_of_hcl++;

    MOLECULE_STRUCT molec;
    molec.creation_time = microseconds();
    molec.max_val = (int)aux;
    molec.molecule_id = no_of_hcl;
    molec.molecule_type = "2HCl";

    if ((fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
    {
        printf("Error opening log.dat\n");
        exit(0);
    }
    write(fs, &molec, sizeof(molec));

    if ((fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
    {
        printf("Error opening log.txt\n");
        exit(0);
    }
    char toWrite[150];
    if (sprintf(toWrite, "\nFORMARE - TipMolecula: %s, NrMolecula: %d, Moment de timp: %lf\n",
                molec.molecule_type, molec.molecule_id, molec.creation_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }

    write(fd, toWrite, strlen(toWrite));

    enqueue(hcl_q, molec.molecule_id);

    int h1 = dequeue(h_q);
    int h2 = dequeue(h_q);
    int cl1 = dequeue(cl_q);
    int cl2 = dequeue(cl_q);

    toWrite[0] = '\0';
    if (sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "hydrogen", h1, "hydrogen", h2) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    if (sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "chlorine", cl1, "chlorine", cl2) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));

    ATOM_STRUCT atom;
    atom.atom_id = h1;
    atom.atom_type = "hydrogen";
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "\nTERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    atom.atom_id = h2;
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    atom.atom_id = cl1;
    atom.atom_type = "chlorine";
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    atom.atom_id = cl2;
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    if (hcl_counter > 0 && feno_counter > 0)
    {
        hcl_counter--;
        feno_counter--;
        reaction((int)aux);
    }

    pthread_mutex_unlock(&hcl_th_mutex);
}

//Creeaza threadul pentru o molecula de 2HCl
void hcl_bond(ATOM_STRUCT *maxi)
{
    pthread_mutex_lock(&hcl_mutex);
    {
        pthread_t hcl_th;
        if (pthread_create(&hcl_th, NULL, hcl_func, (void *)(long)maxi->max_val))
        {
            printf("Error creating hcl thread\n");
            exit(0);
        }
        if (pthread_join(hcl_th, NULL))
        {
            printf("Error joining hcl thread\n");
            exit(0);
        }
    }
    pthread_mutex_unlock(&hcl_mutex);
}

pthread_mutex_t feno_mutex = PTHREAD_MUTEX_INITIALIZER;
int no_of_feno = 0;
pthread_mutex_t feno_th_mutex = PTHREAD_MUTEX_INITIALIZER;

//Creeaza o molecula de Fe(NO3)2 si afiseaza atomii din compozitia acesteia
void *feno_func(void *sent_value)
{
    long aux = (long)sent_value;

    pthread_mutex_lock(&feno_th_mutex);

    feno_counter++;
    no_of_feno++;

    MOLECULE_STRUCT molec;
    molec.creation_time = microseconds();
    molec.max_val = (int)aux;
    molec.molecule_id = no_of_feno;
    molec.molecule_type = "Fe(NO3)2";

    if ((fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
    {
        printf("Error opening log.dat\n");
        exit(0);
    }
    write(fs, &molec, sizeof(molec));

    if ((fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
    {
        printf("Error opening log.txt\n");
        exit(0);
    }
    char toWrite[150];
    if (sprintf(toWrite, "\nFORMARE - TipMolecula: %s, NrMolecula: %d, Moment de timp: %lf\n",
                molec.molecule_type, molec.molecule_id, molec.creation_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }

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
    if (sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n", "iron", fe1) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    if (sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "nitrogen", n1, "nitrogen", n2) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));

    toWrite[0] = '\0';
    if (sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "oxygen", o1, "oxygen", o2) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    toWrite[0] = '\0';
    if (sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "oxygen", o3, "oxygen", o4) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    toWrite[0] = '\0';
    if (sprintf(toWrite, "->TipAtom: %s, IdAtom: %d\n->TipAtom: %s, IdAtom: %d\n", "oxygen", o5, "oxygen", o6) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));

    ATOM_STRUCT atom;
    atom.atom_type = "iron";
    atom.atom_id = fe1;
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "\nTERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    atom.atom_type = "nitrogen";
    atom.atom_id = n1;
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    atom.atom_id = n2;
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    atom.atom_type = "oxygen";
    atom.atom_id = o1;
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    atom.atom_id = o2;
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    atom.atom_id = o3;
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    atom.atom_id = o4;
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    atom.atom_id = o5;
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    atom.atom_id = o6;
    atom.end_time = microseconds();

    toWrite[0] = '\0';
    if (sprintf(toWrite, "TERMINARE_ATOM - TipAtom: %s, IdAtom: %d, Moment de timp: %lf\n",
                atom.atom_type, atom.atom_id, atom.end_time) < 0)
    {
        printf("Error at sprintf\n");
        exit(0);
    }
    write(fd, toWrite, strlen(toWrite));
    write(fs, &atom, sizeof(atom));

    if (hcl_counter > 0 && feno_counter > 0)
    {
        feno_counter--;
        hcl_counter--;
        reaction(aux);
    }

    pthread_mutex_unlock(&feno_th_mutex);
}

//Creeaza threadul pentru o molecula de Fe(NO3)2
void feno_bond(ATOM_STRUCT *maxi)
{
    pthread_mutex_lock(&feno_mutex);
    {
        pthread_t feno_th;
        if (pthread_create(&feno_th, NULL, feno_func, (void *)(long)maxi->max_val))
        {
            printf("Error creating feno thread\n");
            exit(0);
        }
        if (pthread_join(feno_th, NULL))
        {
            printf("Error joining feno thread\n");
            exit(0);
        }
    }
    pthread_mutex_unlock(&feno_mutex);
}

//Creeaza un atom de hidrogen si o pune in asteptare pana cand poate fi folosita la o molecula
//Daca conditiile pentru formarea unei molecule sunt indeplinite se da sem_post
//pentru un numar suficient de mare de atomi si se apeleaza functia de hcl_bond
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

//Creeaza un atom de clor si o pune in asteptare pana cand poate fi folosita la o molecula
//Daca conditiile pentru formarea unei molecule sunt indeplinite se da sem_post
//pentru un numar suficient de mare de atomi si se apeleaza functia de hcl_bond
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

//Creeaza un atom de fier si o pune in asteptare pana cand poate fi folosita la o molecula
//Daca conditiile pentru formarea unei molecule sunt indeplinite se da sem_post
//pentru un numar suficient de mare de atomi si se apeleaza functia de feno_bond
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

//Creeaza un atom de azot si o pune in asteptare pana cand poate fi folosita la o molecula
//Daca conditiile pentru formarea unei molecule sunt indeplinite se da sem_post
//pentru un numar suficient de mare de atomi si se apeleaza functia de feno_bond
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

//Creeaza un atom de oxigen si o pune in asteptare pana cand poate fi folosita la o molecula
//Daca conditiile pentru formarea unei molecule sunt indeplinite se da sem_post
//pentru un numar suficient de mare de atomi si se apeleaza functia de feno_bond
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
        number_of_fecl2_molecules = atoi(argv[1]);
    }
    else if (argc == 3)
    {
        number_of_fecl2_molecules = atoi(argv[2]);
    }
    else
    {
        printf("Numar gresit de argumente\n");
        exit(0);
    }

    //Genereaza timpul de start
    //Timpul pentru o operatie este diferenta
    //dintre timpul cand a fost apelata functia microseconds
    //si acest timp de start, care porneste din momentul
    //de cand am introdus executia programului in terminal
    system(argv[argc]);
    if (clock_gettime(CLOCK_REALTIME, &start) == -1)
    {
        perror("clock gettime");
        return EXIT_FAILURE;
    }

    char *empty = " ";
    if ((fs = open("log.dat", O_CREAT | O_TRUNC | O_RDWR, 0777)) < 0)
    {
        printf("Error opening log.dat\n");
        exit(0);
    }
    write(fs, "", sizeof(empty));

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
        exit(0);
    }

    if (sem_init(&chlorine, 0, 0) == -1)
    {
        perror("error initilalizing ch semaphore\n");
        exit(0);
    }

    if (sem_init(&iron, 0, 0) == -1)
    {
        perror("error initilalizing fe semaphore\n");
        exit(0);
    }

    if (sem_init(&nitrogen, 0, 0) == -1)
    {
        perror("error initilalizing semaphore\n");
        exit(0);
    }

    if (sem_init(&oxygen, 0, 0) == -1)
    {
        perror("error initilalizing o semaphore\n");
        exit(0);
    }

    int h_nr = 2 * number_of_fecl2_molecules;
    int cl_nr = 2 * number_of_fecl2_molecules;
    int fe_nr = 1 * number_of_fecl2_molecules;
    int n_nr = 2 * number_of_fecl2_molecules;
    int o_nr = 6 * number_of_fecl2_molecules;

    ATOM_STRUCT hydrogen_atom[h_nr];
    ATOM_STRUCT chlorine_atom[cl_nr];
    ATOM_STRUCT iron_atom[fe_nr];
    ATOM_STRUCT nitrogen_atom[n_nr];
    ATOM_STRUCT oxygen_atom[o_nr];

    pthread_t h_th[h_nr], cl_th[cl_nr], fe_th[fe_nr], n_th[n_nr], o_th[o_nr];

    int ok = 0;

    for (int i = 0; i < h_nr; i++)
    {

        hydrogen_atom[i].max_val = number_of_fecl2_molecules;
        hydrogen_atom[i].atom_id = i;
        hydrogen_atom[i].atom_type = "hydrogen";
        hydrogen_atom[i].creation_time = microseconds();

        if (pthread_create(&h_th[i], NULL, hydrogen_func, (void *)&hydrogen_atom[i]))
        {
            printf("Error creating h thread\n");
            exit(0);
        }

        if (ok == 0)
        {
            if ((fd = open("log.txt", O_CREAT | O_TRUNC | O_RDWR, 0777)) < 0)
            {
                printf("Error opening log.txt\n");
                exit(0);
            }
        }
        else
        {
            if ((fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
            {
                printf("Error opening log.txt\n");
                exit(0);
            }
        }
        ok = 1;
        char aux[150];
        if (sprintf(aux, "CREARE - TipAtom: %s, IdAtom: %d,  Moment de timp: %lf \n",
                    hydrogen_atom[i].atom_type, hydrogen_atom[i].atom_id, hydrogen_atom[i].creation_time) < 0)
        {
            printf("Error during sprintf\n");
            exit(0);
        }

        write(fd, aux, strlen(aux));
    }

    if ((fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
    {
        printf("Error opening log.dat\n");
        exit(0);
    }
    write(fs, &hydrogen_atom, sizeof(hydrogen_atom));

    for (int i = 0; i < cl_nr; i++)
    {

        chlorine_atom[i].max_val = number_of_fecl2_molecules;
        chlorine_atom[i].atom_id = i;
        chlorine_atom[i].atom_type = "chlorine";
        chlorine_atom[i].creation_time = microseconds();

        if (pthread_create(&cl_th[i], NULL, chlorine_func, (void *)&chlorine_atom[i]))
        {
            printf("Error creating cl thread\n");
            exit(0);
        }

        fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
        char aux[150];
        if (sprintf(aux, "CREARE - TipAtom: %s, IdAtom: %d, Moment de timp: %lf \n",
                    chlorine_atom[i].atom_type, chlorine_atom[i].atom_id, chlorine_atom[i].creation_time) < 0)
        {
            printf("Error during sprintf\n");
            exit(0);
        }

        write(fd, aux, strlen(aux));
    }

    if ((fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
    {
        printf("Error opening log.dat\n");
        exit(0);
    }
    write(fs, &chlorine_atom, sizeof(chlorine_atom));

    for (int i = 0; i < fe_nr; i++)
    {

        iron_atom[i].max_val = number_of_fecl2_molecules;
        iron_atom[i].atom_id = i;
        iron_atom[i].atom_type = "iron";
        iron_atom[i].creation_time = microseconds();
        if (pthread_create(&fe_th[i], NULL, iron_func, (void *)&iron_atom[i]))
        {
            printf("Error creating cl thread\n");
            exit(0);
        }

        fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
        char aux[150];
        if (sprintf(aux, "CREARE - TipAtom: %s, IdAtom: %d, Moment de timp: %lf \n",
                    iron_atom[i].atom_type, iron_atom[i].atom_id, iron_atom[i].creation_time) < 0)
        {
            printf("Error during sprintf\n");
            exit(0);
        }
        write(fd, aux, strlen(aux));
    }

    if ((fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
    {
        printf("Error opening log.dat\n");
        exit(0);
    }
    write(fs, &iron_atom, sizeof(iron_atom));

    for (int i = 0; i < n_nr; i++)
    {

        nitrogen_atom[i].max_val = number_of_fecl2_molecules;
        nitrogen_atom[i].atom_id = i;
        nitrogen_atom[i].atom_type = "nitrogen";
        nitrogen_atom[i].creation_time = microseconds();
        if (pthread_create(&n_th[i], NULL, nitrogen_func, (void *)&nitrogen_atom[i]))
        {
            printf("Error creating n thread\n");
            exit(0);
        }

        fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
        char aux[150];
        if (sprintf(aux, "CREARE - TipAtom: %s, IdAtom: %d, Moment de timp: %lf \n",
                    nitrogen_atom[i].atom_type, nitrogen_atom[i].atom_id, nitrogen_atom[i].creation_time) < 0)
        {
            printf("Error during sprintf\n");
            exit(0);
        }

        write(fd, aux, strlen(aux));
    }

    if ((fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
    {
        printf("Error opening log.txt\n");
        exit(0);
    }
    write(fs, &nitrogen_atom, sizeof(nitrogen_atom));

    for (int i = 0; i < o_nr; i++)
    {

        oxygen_atom[i].max_val = number_of_fecl2_molecules;
        oxygen_atom[i].atom_id = i;
        oxygen_atom[i].atom_type = "oxygen";
        oxygen_atom[i].creation_time = microseconds();
        if (pthread_create(&o_th[i], NULL, oxygen_func, (void *)&oxygen_atom[i]))
        {
            printf("Error creating o thread\n");
            exit(0);
        }

        fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777);
        char aux[150];
        if (sprintf(aux, "CREARE - TipAtom: %s, IdAtom: %d, Moment de timp: %lf \n",
                    oxygen_atom[i].atom_type, oxygen_atom[i].atom_id, oxygen_atom[i].creation_time) < 0)
        {
            printf("Error during sprintf\n");
            exit(0);
        }

        write(fd, aux, strlen(aux));
    }

    if ((fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
    {
        printf("Error opening log.dat\n");
        exit(0);
    }
    write(fs, &oxygen_atom, sizeof(oxygen_atom));

    for (int i = 0; i < h_nr; i++)
        if (pthread_join(h_th[i], NULL))
        {
            printf("Error joining thread\n");
            exit(0);
        }

    for (int i = 0; i < cl_nr; i++)
        if (pthread_join(cl_th[i], NULL))
        {
            printf("Error joining thread\n");
            exit(0);
        }

    for (int i = 0; i < fe_nr; i++)
        if (pthread_join(fe_th[i], NULL))
        {
            printf("Error joining thread\n");
            exit(0);
        }

    for (int i = 0; i < n_nr; i++)
        if (pthread_join(n_th[i], NULL))
        {
            printf("Error joining thread\n");
            exit(0);
        }

    for (int i = 0; i < o_nr; i++)
        if (pthread_join(o_th[i], NULL))
        {
            printf("Error joining thread\n");
            exit(0);
        }

    if (big_equation_counter == number_of_fecl2_molecules)
    {
        char aux[150];
        if (sprintf(aux, "\nCONDITIE_DE_FINAL - %lf\n", microseconds()) < 0)
        {
            printf("Error during sprintf\n");
            exit(0);
        }
        if ((fd = open("log.txt", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
        {
            printf("Error opening log.txt\n");
            exit(0);
        }
        if ((write(fd, aux, strlen(aux))) < 0)
        {
            printf("Error writing in log.txt\n");
            exit(0);
        }

        if ((fs = open("log.dat", O_CREAT | O_APPEND | O_RDWR, 0777)) < 0)
        {
            printf("Error opening log.txt\n");
            exit(0);
        }
        if ((write(fs, aux, strlen(aux))) < 0)
        {
            printf("Error writing in log.dat\n");
            exit(0);
        }
    }

    sem_destroy(&hydrogen);
    sem_destroy(&chlorine);
    sem_destroy(&iron);
    sem_destroy(&nitrogen);
    sem_destroy(&oxygen);

    close(fd);
}
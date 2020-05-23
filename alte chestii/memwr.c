#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <string.h>

void write_with_mmap(char *input){

    const char *path = "log.dat";

    int fd = open(path, O_RDWR | O_CREAT, 0644);
    
    if (fd == -1)
    {
        perror("Error opening file\n");
        exit(0);
    }

    int size = strlen(input); 
    
    if (lseek(fd, size - 1, SEEK_SET) == -1)
    {
        perror("Error at lseek\n");
        exit(0);
    }
    
    write(fd, "", 1);
  
    char *log_map = mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
    if (log_map == MAP_FAILED)
    {
        perror("Error at mmap\n");
        exit(0);
    }
    
    for (int i = 0; i < size; i++)
        log_map[i] = input[i];

    close(fd);
    
}

int main(int argc, const char *argv[])
{   
    char primu[1000000]= "guasddas\nsadad\nasdasd\n";
    strcat(primu,"gioni");
    write_with_mmap(primu);
    
    return 0;
}
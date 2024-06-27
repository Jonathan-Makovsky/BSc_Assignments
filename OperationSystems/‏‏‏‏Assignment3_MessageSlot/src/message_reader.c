#include <fcntl.h>     
#include <unistd.h>     
#include <sys/ioctl.h>  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "message_slot.h"

int main(int argc, char *argv[]){
    if (argc != 3){
        perror("Error");
        exit(1);
    }
    
    unsigned int id = atoi(argv[2]);
    char message[128];
    int size = 0;

    int file_directory = open(argv[1], O_RDWR);
    if (file_directory < 0){
        perror("Error while doing open");
        exit(1);
    }
    if (ioctl(file_directory, MSG_SLOT_CHANNEL, id) != 0){
		perror("Error while doing ioctl");
		close(file_directory);
        exit(1);
    }
    size = read(file_directory, message, 128);
    if (size <= 0){
        perror("Error while doing read");
		close(file_directory);
        exit(1);
    }
    close(file_directory);
    if (write(1, message, size) <= 0){
        perror("Error while writing the message to the screen");
        exit(1);
    }
    exit(0);
}

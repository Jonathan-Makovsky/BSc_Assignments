#include <fcntl.h>     
#include <unistd.h>     
#include <sys/ioctl.h>  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include "message_slot.h"

int main(int argc, char *argv[]){
    if (argc != 4){
        perror("Error");
        exit(1);
    }
    
    unsigned int id = atoi(argv[2]);
    char *message = argv[3];
	
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
    
    if (write(file_directory, message, strlen(message)) <= 0){
         perror("Error while doing write");
		close(file_directory);
        exit(1);
    }
    
    close(file_directory);
    exit(0);
}

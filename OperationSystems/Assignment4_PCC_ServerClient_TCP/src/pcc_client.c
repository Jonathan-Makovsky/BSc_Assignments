#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <arpa/inet.h>
#include <stdint.h>

// Macro - size of Megabyte ****************
#define MEGA 1048576

int main(int argc, char *argv[])
{
  // General properties ******************
  int sockfd = -1;

  // For the writing proccess ********************
  int str_bytes_sent = 0;
  int curr_write = 0;
  int curr_buff_write = 0;
  int bytes_left = 0;
  char temp_buff[MEGA];

  // Retreiving C - number of printable characters
  uint32_t C = 0;
  int C_bytes_read = 0;

  // Socket atributes ****************
  struct sockaddr_in serv_addr; // where we Want to get to
  
  if (argc != 4)
  {
    perror("Number of arguments is incorrect");
    exit(1);
  }
  
	char *path = argv[3];
	FILE *fp = fopen(path, "r");
	if (fp == NULL)
	{
		error("error - file open didn't work!");
		exit(1);
	}
	
	// Retreiving N *************************
	int nsent;
	fseek(fp, 0L, SEEK_END);
	uint32_t N = ftell(fp); // N size
	rewind(fp); // fp reset to the beginning of the file
	

  if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
  {
    perror("Socket process didn't work!");
    exit(1);
  }
  
  memset(&serv_addr, 0, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_port = htons((unsigned short)atoi(argv[2]));  
  
  if (inet_pton(AF_INET, argv[1], &(serv_addr.sin_addr.s_addr)) != 1){
  perror("Error while converting ip address");
  exit(1);
  }
  
  // connect socket to the target address
  if (connect(sockfd,
              (struct sockaddr *)&serv_addr,
              sizeof(serv_addr)) < 0)
  {
    perror("connect didn't work well!\n");
    exit(1);
  }
	
  // Write N to the server ********************
  N = htonl(N);
  nsent = write(sockfd, &N, sizeof(uint32_t));
  if (nsent <= 0)
  {
    perror("Write didn't work well!");
    exit(1);
  }

  // Writing the bytes stream to the server ***********************
  while (0 < (curr_write = fread(temp_buff, 1, MEGA, fp)))
  {
    curr_buff_write = 0;
    bytes_left = curr_write;
    while (1)
    {
      str_bytes_sent = write(sockfd, temp_buff + curr_buff_write, bytes_left);
      if (str_bytes_sent <= 0)
      {
        perror("Write didn't work well!");
        exit(1);
      }
      curr_buff_write += str_bytes_sent;
      bytes_left -= str_bytes_sent;

      if (curr_buff_write == curr_write) // We wrote the entire buffer
      {
        break;
      }
    }
  }

  // Retreiving C from the server **************
  C_bytes_read = read(sockfd, &C, sizeof(uint32_t));
  if (C_bytes_read <= 0)
  {
    perror("Error while read!");
    exit(1);
  }

  C = ntohl(C);
  
  close(sockfd); 
  
  fprintf(stdout, "# of printable characters: %u\n", C);
  
  exit(0);
}

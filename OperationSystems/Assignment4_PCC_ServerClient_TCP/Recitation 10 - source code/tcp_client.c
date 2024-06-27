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
#include <errno.h>

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
  unsigned long C = 0;
  int C_bytes_read = 0;

  // Socket atributes ****************
  struct sockaddr_in serv_addr; // where we Want to get to
  struct sockaddr_in my_addr;   // where we actually connected through
  struct sockaddr_in peer_addr; // where we actually connected to
  socklen_t addrsize = sizeof(struct sockaddr_in);
  
  if (argc != 4)
  {
    perror("Number of arguments is incorrect");
    exit(1);
  }

  char *path = argv[3];
  FILE *fp = fopen(path, "r");
  if (fp == NULL)
  {
    perror("error - file open didn't work!");
    exit(1);
  }

  // Retreiving N *************************
  int nsent;
  fseek(fp, 0L, SEEK_END);
  unsigned long N = ftell(fp); // N size
  rewind(fp); // fp reset to the beginning of the file

  if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) < 0)
  {
    perror("Socket process didn't work!");
    exit(1);
  }

  // print socket details ***************
  getsockname(sockfd,
              (struct sockaddr *)&my_addr,
              &addrsize);
  printf("Client: socket created %s:%d\n",
         inet_ntoa((my_addr.sin_addr)),
         ntohs(my_addr.sin_port));

  memset(&serv_addr, 0, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_port = htons(10000);                  // Note: htons for endiannes
  serv_addr.sin_addr.s_addr = inet_addr("127.0.0.1"); // hardcoded...

  printf("Client: connecting...\n");
  // Note: what about the client port number?
  // connect socket to the target address
  if (connect(sockfd,
              (struct sockaddr *)&serv_addr,
              sizeof(serv_addr)) < 0)
  {
    printf("\n Error : Connect Failed. %s \n", strerror(errno));
    return 1;
  }

  // print socket details again
  getsockname(sockfd, (struct sockaddr *)&my_addr, &addrsize);
  getpeername(sockfd, (struct sockaddr *)&peer_addr, &addrsize);
  printf("Client: Connected. \n"
         "\t\tSource IP: %s Source Port: %d\n"
         "\t\tTarget IP: %s Target Port: %d\n",
         inet_ntoa((my_addr.sin_addr)), ntohs(my_addr.sin_port),
         inet_ntoa((peer_addr.sin_addr)), ntohs(peer_addr.sin_port));

  // Write N to the server ********************
  N = htonl(N);
  nsent = write(sockfd, &N, sizeof(unsigned long));
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
  C_bytes_read = read(sockfd, &C, unsigned long);
  if (C_bytes_read <= 0)
  {
    perror("Error while read!");
    exit(1);
  }

  C = htonl(C);
  close(sockfd); 
  return 0;
}

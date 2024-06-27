#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <assert.h>

// Macro - size of Megabyte ****************
#define MEGA 1048576

// Global *****************
unsigned long pcc_total[95];
int flag = 0;
int is_ok = 1;

// Function printing the stats *******************
void print_printable()
{
  int i = 0;
  for (i = 0; i < 95; i++)
  {
    printf("char '%c' : %u times\n", i + 32, pcc_total[i]);
  }
}

// Function to handel sigint ************************
static void sig_handler(int sig)
{
  if (flag != 1)
  {
    print_prinatble();
    exit(0);
  }
  is_ok = 0;
}

// Function that sets the properties of handling sigints
int sig_prop()
{
  struct sigaction sa;
  sa.sa_handler = sig_handler;
  sa.sa_flags = SA_RESTART;
  if (sigaction(SIGINT, &sa, NULL) == -1)
  {
    return -1;
  }
  return 0;
}

// main ********************************************************
int main(int argc, char *argv[])
{
  // General properties ******************
  int listenfd = -1;
  int connfd = -1;

  // Socket atributes ****************
  struct sockaddr_in serv_addr;
  struct sockaddr_in my_addr;
  struct sockaddr_in peer_addr;
  socklen_t addrsize = sizeof(struct sockaddr_in);

  // For the Reading of the bytes stream proccess ********************
  char recv_buff[MEGA];
  int curr = 0;
  int i = 0;
  unsigned long temp[95]; // For atomic reading of the bytes stream
  int bytes_read = 0;
  int bytes_left = 0;
  int total_bytes_read = 0;
  int str_bytes_read = 0;
  int curr_read = 0;

  // C - number of printable characters properties ****************
  unsigned long C = 0;
  int Csent = 0;

  // Retreiving N *************************
  unsigned long N;

  if (argc != 2)
  {
    perror("Number of arguments is incorrect");
    exit(1);
  }

  sig_prop(); // Sigint properties are now activated

  listenfd = socket(AF_INET, SOCK_STREAM, 0);
  if (listenfd == -1)
  {
    perror("socket initialize error");
    exit(1);
  }
  memset(&serv_addr, 0, addrsize);

  serv_addr.sin_family = AF_INET;
  // INADDR_ANY = any local machine address
  serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
  serv_addr.sin_port = htons(10000);

  if (0 != bind(listenfd,
                (struct sockaddr *)&serv_addr,
                addrsize))
  {
    perror("bind initialize error");
    exit(1);
  }
  if (0 != listen(listenfd, 10))
  {
    perror("listen initialize error");
    exit(1);
  }
  while (is_ok)
  {
    C = 0;
    unsigned long temp[95];
    // Accept a connection.
    connfd = accept(listenfd,
                    (struct sockaddr *)&peer_addr,
                    &addrsize);
    if (connfd < 0)
    {
      perror("accept initialize error");
      exit(1);
    }
    flag = 1; // A client is being processed

    bytes_read = read(connfd, &N, sizeof(unsigned long)); // Reading N

    if (bytes_read < 0)
    {
      if (errno == EINTR || errno == ETIMEDOUT || errno == ECONNRESET || errno == EPIPE)
      {
        flag = 0;
        continue;
      }
    }
    N = htonl(N);
    bytes_left = N;
    total_bytes_read = 0;
    str_bytes_read = 0;

    while (1)
    {
      curr_read = MEGA;
      if (bytes_left < MEGA)
      {
        curr_read = bytes_left;
      }
      str_bytes_read = read(connfd, recv_buff, curr_read);
      if (bytes_read < 0)
      {
        if (errno == EINTR || errno == ETIMEDOUT || errno == ECONNRESET || errno == EPIPE)
        {
          flag = 0;
          continue;
        }
      }
      bytes_left -= str_bytes_read;
      total_bytes_read += str_bytes_read;

      for (i = 0; i < str_bytes_read; i++)
      {
        curr = recv_buff[i];
        if (32 <= curr && curr <= 126)
        {
          C++;
          temp[curr - 32]++;
        }
      }

      if (total_bytes_read == N) // We read the entire buffer
      {
        break;
      }
    }

    C = htonl(C);
    
    // Copying the printable character temp array to the original one **************
    for (i = 0; i < 95; i++)
    {
      pcc_total[i] = temp[i];
    }

    // Sending C to the client ************
    Csent = write(connfd, &C, unsigned long);
    if (Csent <= 0)
    {
      if (errno == EINTR || errno == ETIMEDOUT || errno == ECONNRESET || errno == EPIPE)
      {
        flag = 0;
        continue;
      }
    }

    // close socket
    close(connfd);
    flag = 0;
  }
  print_prinatble();
  exit(0);
}

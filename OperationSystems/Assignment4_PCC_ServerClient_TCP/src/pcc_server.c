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
#include <stdint.h>
#include <sys/wait.h>

// Macro - size of Megabyte ****************
#define MEGA 1048576

// Global *****************
uint32_t pcc_total[95];
int flag = 0;
int is_ok = 1;
int flag_read = 1; //Flag for when bytes stream read failes and we dont want to exit

void print_printable();

// Function printing the stats *******************
void print_printable()
{
  int i = 0;
  for (i = 0; i < 95; i++)
  {
    fprintf(stdout, "char '%c' : %u times\n", i + 32, pcc_total[i]);
  }
}

// Function to handel sigint ************************
static void sig_handler(int sig)
{
  if (flag != 1)
  {
    print_printable();
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
  //struct sockaddr_in my_addr;
  struct sockaddr_in peer_addr;
  socklen_t addrsize = sizeof(struct sockaddr_in);

  // For the Reading of the bytes stream proccess ********************
  char recv_buff[MEGA];
  int curr = 0;
  int i = 0;
  int bytes_read = 0;
  int bytes_left = 0;
  int total_bytes_read = 0;
  int str_bytes_read = 0;
  int curr_read = 0;
  uint32_t temp[95];

  // C - number of printable characters properties ****************
  uint32_t C = 0;
  int Csent = 0;

  // Retreiving N *************************
  uint32_t N;

  if (argc != 2)
  {
    perror("Number of arguments is incorrect");
    exit(1);
  }

  if(sig_prop()!=0){ // Sigint properties are now activated
  	perror("sigint error");
    exit(1);
  }

  listenfd = socket(AF_INET, SOCK_STREAM, 0);
  if (listenfd == -1)
  {
    perror("socket initialize error");
    exit(1);
  }
  
  if (setsockopt(listenfd, SOL_SOCKET, SO_REUSEADDR, &(int){1}, sizeof(int)) < 0){
    	perror("setsockopt error");
    	exit(1);
  }
    
  memset(&serv_addr, 0, addrsize);

  serv_addr.sin_family = AF_INET;
  // INADDR_ANY = any local machine address
  serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
  serv_addr.sin_port = htons((unsigned short)atoi(argv[1]));

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
    memset(temp, 0, 95*4);
    flag_read = 1;
    // Accept a connection.
    connfd = accept(listenfd,(struct sockaddr *)&peer_addr, &addrsize);
    if (connfd < 0)
    {
      perror("accept initialize error");
      exit(1);
    }
    
    flag = 1; // A client is being processed

    bytes_read = read(connfd, &N, 4); // Reading N
    if (bytes_read <= 0)
    {
      if (bytes_read == 0 || errno == ETIMEDOUT || errno == ECONNRESET || errno == EPIPE)
      {
      	close(connfd);
        flag = 0;
        continue;
      }
      perror("read N error");
      exit(1);
      
    }
    
    N = ntohl(N);
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
      if (str_bytes_read <= 0)
      {
        if (str_bytes_read == 0 || errno == ETIMEDOUT || errno == ECONNRESET || errno == EPIPE)
        {
          flag = 0;
          flag_read = 0;
          break;
        }
        perror("read bytes_stream error");
      	exit(1);
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

      if (total_bytes_read == N) // We have read the entire buffer
      {
        break;
      }
    }
	if (flag_read == 1){
    	// Sending C to the client ************
    	C = htonl(C);
    	Csent = write(connfd, &C, sizeof(uint32_t));
    	if (Csent <= 0)
    	{
      		if (errno == ETIMEDOUT || errno == ECONNRESET || errno == EPIPE)
      		{
      			close(connfd);
        		flag = 0;
        		continue;
      		}
      		else{
      			perror("read bytes_stream error");
      			exit(1);
      		}
    	}
    	
    	// Copying the printable character temp array to the original one **************
    	for (i = 0; i < 95; i++)
		{
    		pcc_total[i] += temp[i];
		}
	}
    // close socket
    close(connfd);
    flag = 0;
    
  }
  print_printable();
  exit(0);
}

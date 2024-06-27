#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>

// arglist - a list of char* arguments (words) provided by the user
// it contains count+1 items, where the last item (arglist[count]) and *only* the last is NULL
// RETURNS - 1 if should continue, 0 otherwise
int process_arglist(int count, char **arglist);

// prepare and finalize calls for initialization and destruction of anything required
int prepare(void);
int finalize(void);

// waitpid decleration
int waitpid(pid_t, int *, int);


int ignore_sig() {
	struct sigaction sa;
	sa.sa_handler = SIG_IGN;
	sa.sa_flags = SA_RESTART;
	if (sigaction(SIGINT, &sa, NULL) == -1) {
        return -1;
    }

    return 0;
}

int acc_sig() {
	struct sigaction sa;
	sa.sa_handler = SIG_DFL;
	sa.sa_flags = SA_RESTART;
	if (sigaction(SIGINT, &sa, NULL) == -1) {
        return -1;
    }

    return 0;
}


// https://stackoverflow.com/questions/7171722/how-can-i-handle-sigchld/7171836#7171836
// THe Eran Trick
static void child_handler(int sig) {
    pid_t wait_pid;
    int status;
	
    /* EEEEXTEERMINAAATE! */
    while((wait_pid = waitpid(-1, &status, WNOHANG)) > 0);
	
    if(wait_pid == -1 && ECHILD != errno && EINTR != errno){
        perror("Zombie error");
        exit(1);
    }
}


int prepare() {
    // The Eran trick
    
    struct sigaction sa;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART;
    sa.sa_handler = child_handler;

    if (sigaction(SIGCHLD, &sa, NULL) == -1){
        perror("error");
        return 1;
    }
    if (ignore_sig() == -1){
        perror("error");
        return 1;
    }

    return 0;
}
int finalize(){
	return 0;
}

//process_arglist implementation
int process_arglist(int count, char **arglist)
{
    int i = 0, num_of_elements = 0, status = 0;
    int stat = 0, stat1 = 0, stat2 = 0, pid = 0, pid1 = 0, pid2 = 0;
    int fork_res = 0;
    char *op = "";

    while (arglist[i] != NULL)
    {
        if (strcmp(arglist[i], ">") == 0)
        {
            op = ">";
            break;
        }
        if (strcmp(arglist[i], "&") == 0)
        {
            op = "&";
            break;
        }
        if (strcmp(arglist[i], "|") == 0)
        {
            op = "|";
            break;
        }
        i++;
        num_of_elements++;
    }
    num_of_elements++;
    
    //////////////  Case where op == "" ////////////////////////
    if (strcmp(op, "") == 0)
    {
    	
        fork_res = fork();
        if (fork_res == -1) 
        {
            perror("fork error");
            return 0;
        }
        if (fork_res == 0)
        {
        	if (acc_sig() == -1){
                perror("error");
                exit(1);
            }
            status = execvp(arglist[0], arglist);
            if (status == -1) 
            {
                perror("Child error");
                exit(1);
            }
        }
        if (fork_res > 0)
        {
            pid = waitpid(fork_res, &stat, 0);
            if (pid == -1 && ECHILD != errno && EINTR != errno){
            	perror("error");
            	return 0;
            }
        }
    }
    //////////////////////  op == "&" ///////////////////
    if (strcmp(op, "&") == 0)
    {
    	
        fork_res = fork();
        if (fork_res == -1)
        {
            perror("fork error");
            return 0;
        }
        if (fork_res == 0)
        {
        	arglist[num_of_elements - 1] = NULL;
            status = execvp(arglist[0], arglist);
            if (status == -1)
            {
                perror("Child error");
                exit(1);
            }
        }
    }
    //////////////  Case where op == "|" ////////////////////////
    if (strcmp(op, "|") == 0)
    {
    	
        int l_counter = num_of_elements - 1;
            
        /// Creating the pipe
        int pipefd[2];
        pid_t cpid;
        pid_t cpid2;
        int k = pipe(pipefd);

        if (-1 == k)
        {
            perror("pipe");
            exit(-1);
        }
        cpid = fork();
 
        if (-1 == cpid)
        {
            perror("fork");
            close(pipefd[0]);
            close(pipefd[1]);
            return 0;
        }
        else if (0 == cpid)
        {
        	if (acc_sig() == -1){
                perror("error");
                exit(1);
            }
            // Child write from pipe
            // Close unused read end  
            close(pipefd[0]);  
            if (-1 == dup2(pipefd[1], 1)){
            	perror("Child error");
            	exit(1);
            }
            
            arglist[l_counter] = NULL;
            status = execvp(arglist[0], arglist);
            if (status == -1)
            {
                perror("Child error");
            	exit(1);
            }
        }
        else if(cpid > 0)
        {
            cpid2 = fork();
            
            if (-1 == cpid2)
            {
                perror("error");
                close(pipefd[0]);
                close(pipefd[1]);
                return 0;
            }
            else if (0 == cpid2)
            {
            	if (acc_sig() == -1){
                	perror("error");
                	exit(1);
            	}
                // Child reads from pipe
                // Close unused write end
                close(pipefd[1]);
                if (-1 == dup2(pipefd[0], 0)){
                	perror("Child error");
                	return 0;
                }
                
                status = execvp(arglist[l_counter + 1],  arglist + l_counter + 1);
                if (status == -1)
                { 
                	perror("Child error");
                	exit(1);
            	}
             }
             else if(cpid2 > 0)
             { // The father process, after we finished with both of the child processes
             
                close(pipefd[0]);
                close(pipefd[1]);
                pid1 = waitpid(cpid, &stat1, 0);
                if (pid1 == -1 && ECHILD != errno && EINTR != errno){
            		perror("error");
            		return 0;
            	}
            	pid2 = waitpid(cpid2, &stat2, 0);
            	if (pid2 == -1 && ECHILD != errno && EINTR != errno){
            		perror("error");
            		return 0;
            	}
                    
            }
		}    
    }
    //////////////  Case where op == ">" ////////////////////////
    if (strcmp(op, ">") == 0)
    {	
        int fd = open(arglist[count - 1], O_WRONLY | O_CREAT | O_TRUNC, S_IRWXU);
        if (fd < 0)
        { 
        	perror("error");
        	return 0;
        }
        fork_res = fork();
        if (fork_res == -1)
        { 
        	perror("fork error");
            return 0;
        }
        if (fork_res == 0)
        {	
        	if (acc_sig() == -1){
                perror("error");
                exit(1);
            }
        	arglist[num_of_elements - 1] = NULL;
            if (-1 == dup2(fd, 1)){ 
            	perror("Child error");
            	exit(1);
            }
            status = execvp(arglist[0], arglist);
            if (status == -1)
            {
                perror("Child error");
            	exit(1);
            }
        }
        if (fork_res > 0)
        {
            pid = waitpid(fork_res, &stat, 0);
            if (pid == -1 && ECHILD != errno && EINTR != errno){
            	perror("error");
            	return 0;
            }
            close(fd);
        }
    }
    return 1;
}



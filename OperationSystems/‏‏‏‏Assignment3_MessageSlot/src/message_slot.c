#undef __KERNEL__
#define __KERNEL__
#undef MODULE
#define MODULE

#include <linux/kernel.h>   /* We're doing kernel work */
#include <linux/module.h>   /* Specifically, a module */
#include <linux/fs.h>       /* for register_chrdev */
#include <linux/uaccess.h>  /* for get_user and put_user */
#include <linux/string.h>   /* for memset. NOTE - not string.h!*/
#include <linux/slab.h>

MODULE_LICENSE("GPL");

//Our custom definitions of IOCTL operations
#include "message_slot.h"
Device device_arr[256];



//================== DEVICE FUNCTIONS ===========================
static int device_open(struct inode* inode, struct file* file) {
    int minor_num = iminor(inode);
    PrivateData* curr_data = kmalloc(sizeof(PrivateData), GFP_KERNEL);
    if (curr_data == NULL){
        return -EIO;
    }
 
    curr_data->minor_num = minor_num;
    file->private_data = (void*)curr_data;
    return 0;
}


//---------------------------------------------------------------
static int device_release(struct inode* inode, struct file* file) {
    PrivateData* curr_data = (PrivateData*)(file->private_data);
    kfree(curr_data);
    return 0;
}

//---------------------------------------------------------------
// a process which has already opened
// the device file attempts to read from it
static ssize_t device_read(struct file* file, char __user* buffer, size_t length, loff_t* offset){
    PrivateData* data = (PrivateData*)file->private_data;
    Channel* channel;
    int msg_size;
    int i = 0;

    if (data == NULL){
        return -EINVAL;
    }
    channel = data->channel;
    if(channel == NULL) {
        return -EINVAL;
    }
    msg_size = channel->msg_size;
    if (msg_size == 0){
        return -EWOULDBLOCK;
    }
    if (length < msg_size){
        return -ENOSPC;
    }
    while (i < msg_size){ 
        if ((put_user((channel->msg)[i], buffer+i)) != 0){
            return -EINVAL;
        }
        i++;
    }
    return msg_size;
}

//---------------------------------------------------------------
// a processs which has already opened
// the device file attempts to write to it
static ssize_t device_write(struct file* file, const char __user* buffer, size_t length, loff_t* offset){
    PrivateData* data = (PrivateData*)(file->private_data);
    Channel* channel;
    int i = 0;
    char new_msg[128];
    
    if (data == NULL){
        return -EINVAL;
    }
    channel = data->channel;
    if(channel == NULL) {
        return -EINVAL;
    }
    if (length <= 0 || 128 < length ){
        return -EMSGSIZE;
    }

    while (i < length){
        if((get_user(new_msg[i], buffer+i)) != 0){
            return -EINVAL;
        }
        i++;
    }
    for (i = 0; i < length; i++){
        (channel->msg)[i] = new_msg[i];
    }
    channel->msg_size = length;
    return length;
}

//----------------------------------------------------------------
static long device_ioctl(struct file* file, unsigned int ioctl_command_id, unsigned long ioctl_param){
    PrivateData* data;
    int minor_num;
    Channel* head;
    Channel* new_channel;
    
   
    if(ioctl_command_id != MSG_SLOT_CHANNEL || ioctl_param == 0){
        return -1;
    }
    data = (PrivateData*)(file->private_data);
    if (data == NULL){
        return -EINVAL;
    }
    minor_num = data->minor_num;
    head = device_arr[minor_num].head;

    if (head == NULL){
        Channel* new_channel = kmalloc(sizeof(Channel), GFP_KERNEL);
        if (new_channel == NULL){
            return -EIO;
        }
        new_channel->channel_id = ioctl_param;
        device_arr[minor_num].head = new_channel;
        data->channel = new_channel;
        data->channel_id = ioctl_param;

        return 0;
    }

    while (head->next != NULL){
        int id = head->channel_id;
        if (id == ioctl_param){
            data->channel = head;
            data->channel_id = ioctl_param;
            return 0;
        }
        head = head->next;
    }

    if (head ->channel_id == ioctl_param){
        data->channel = head;
        data->channel_id = ioctl_param;
        return 0;
    }

    new_channel = kmalloc(sizeof(Channel), GFP_KERNEL);
    if (new_channel == NULL){
        return -EIO;
    }
    new_channel->channel_id = ioctl_param;
    head->next = new_channel;
    data->channel = new_channel;
    data->channel_id = ioctl_param;
    
    return 0;
}
//==================== DEVICE SETUP =============================

// This structure will hold the functions to be called
// when a process does something to the device we created
struct file_operations Fops = {
        .owner	  = THIS_MODULE,
        .read           = device_read,
        .write          = device_write,
        .open           = device_open,
        .unlocked_ioctl = device_ioctl,
        .release        = device_release,
};

//---------------------------------------------------------------
// Initialize the module - Register the character device
static int __init simple_init(void) {
    int rc = -1;

    // Register driver capabilities. Obtain major num
    rc = register_chrdev( MAJOR_NUM, DEVICE_RANGE_NAME, &Fops );

    // Negative values signify an error
    if( rc < 0 ) {
        printk(KERN_ALERT "%s registraion failed for  %d\n", DEVICE_RANGE_NAME, MAJOR_NUM );
        return rc;
    }
    /* Making sure that the devices are initialize*/
 

    return SUCCESS;
}

//---------------------------------------------------------------
static void __exit simple_cleanup(void) {
    // Unregister the device
    // Should always succeed
	int i = 0;
    for (i = 0; i < 256; i++){
        Device curr_device = device_arr[i];
        Channel* curr_channel = curr_device.head;
        while (curr_channel != NULL){
            Channel* next_channel = curr_channel->next;
            kfree(curr_channel);
            curr_channel = next_channel;
        }

    }
    unregister_chrdev(MAJOR_NUM, DEVICE_RANGE_NAME);
}


//---------------------------------------------------------------
module_init(simple_init);
module_exit(simple_cleanup);

//========================= END OF FILE =========================

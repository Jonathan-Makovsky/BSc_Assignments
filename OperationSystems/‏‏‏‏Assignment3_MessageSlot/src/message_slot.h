#ifndef MESSAGE_SLOT_H
#define MESSAGE_SLOT_H

#include <linux/ioctl.h>

#define MAJOR_NUM (235)
#define MSG_SLOT_CHANNEL _IOW(MAJOR_NUM, 0, unsigned int)

#define DEVICE_RANGE_NAME "message_slot"
#define BUF_LEN 128
#define SUCCESS 0

typedef struct Channel{
    int channel_id;
    char msg[128];
    int msg_size;
    struct Channel* next;
} Channel;

typedef struct Device{
   struct Channel* head;
   int minor_num;
} Device;  



typedef struct PrivateData{
    int minor_num;
    int channel_id;
    struct Channel* channel;
} PrivateData;



#endif

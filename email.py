#!/usr/bin/env python3
from envelope import Envelope
def send_report():
 
    Envelope()\
        .from_("hpcreports@cmd.su")\
        .subject("[hpc_tasks] report")\
        .to("schelkunov@cmd.su")\
        .message("report")\
        .smtp(host="email.cmd.su", port="25", user="hpcreports", password="tnXmk5k7KLy8", security="None", timeout=30, attempts=3, delay=3)\
        .send()
 
 
if __name__ == '__main__':
    send_report()
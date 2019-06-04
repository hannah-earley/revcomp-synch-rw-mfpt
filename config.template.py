#!/usr/bin/env python3

# Email configuration for batch status updates
class email:
    argv = ['/usr/sbin/sendmail', '-t', '-oi']
    recipient = 'Me <user@localhost>'
    sender = 'RWX Batch <user@localhost>'
    subject = 'Re: RWX Batch Update'
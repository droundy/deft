# General cheat sheet


Terminal commands:
-----------------

    2>&1 
takes the stderr of your program (one of the types of output) and puts it in stdout, along with what would also be in stdout.  The idea is that the "2" represents the stderr and the "&1" represents the fact that "1" is a file descriptor and "1" itself is the stdout.

    valgrind *program* | less
Will run program in valgrind and pipe results in less so can read them.  Valgrind is named after the main entrance to Valhalla, which is ruled by Odin and is where warriors go when they die in combat.  Runs the program on a virtual machine after it has converted it into a language called Intermediate Representation.  Once in IR, it can use "tools" on it, that do stuff that's nice for us when debugging.  Most common is "memcheck".  This keeps track of allocation of memory - it let's you know when memory is unallocated or "undefined" and also whether a memory address points to an allocated, non-freed memory block. Output is in the form:
   
    ==2088== Conditional jump or move depends on uninitialised value(s)
    ==2088==    at 0x56274A0: __printf_fp (printf_fp.c:404)
    ==2088==    by 0x562396A: vfprintf (vfprintf.c:1622)
    ==2088==    by 0x562C927: fprintf (fprintf.c:33)
    ==2088==    by 0x402D10: main (protein.cpp:307)
The problem in your code is given in the last line (the protein.cpp in example above)
 
    git diff
This views the changes you made in your working tree relative to the index (staging area for the next commit).  There are different forms or options with git diff that allow you to compare different things.  Worth looking up. Git has nothing to do with Odin *or* warriors.

## ssh keys

Private keys are used to decrypt messages sent over ssh.  Every host
has a private key.  Your account on quipu is a host that has a private
key.  This private key is never given to any other host.  The private
key on your host is never stored in decrypted form on your host's hard
disk.  It's encrypted form is stored on your hard disk and in order to
decrypt it you need to use your pass phrase.  You pass phrase will be
used by an ssh agent, which is a program on your operating system,
to decrypt your private key, which will only be stored in memory and
will be deleted once the ssh agent is done using it.  This way there
are never any decrypted private keys laying around on hard disk.

Associated with your private key is a public key.  This public key can
be shown to anyone.  When your host ssh-es into another host, you send
a public key, and that other host checks this public key against its
list of authorized keys to see if it will allow you to set up an ssh
channel with it.

There is a directory in your `/home/youraccount/` directory that is
called `.ssh`.  This directory contains four files:

* `authorized_keys` -
the public keys of hosts that your ssh agent will allow to ssh into
your host. (these other hosts are asking to ssh into your computer,
and if there public key is listed in this file, a connection will be
established)
* `id_rsa` - this is the encrypted form of the private key
that is stored on your hard disk.  The ssh agent demands a pass phrase
that it uses to decrypt this private key when it is establishing an
ssh connection with a remote host
* `id_rsa.pub` - this is your host's
public key.  This key is asked for by remote hosts when you try to set
up an ssh connection with a remote host.
* `known_hosts` - these are the
public keys of hosts that you have sshed into.  After you have used
your pass phrase once to decrypt your private key, you do not need to
re-enter your password in order to ssh into any of the hosts on this
list.

The command line program `ssh-add` will add a private key to your agent.

If I were an evil person >:) I would steal a laptop, then spend a long
time guessing at the pass phrase of the host of that laptop.  This
would allow me to dycrypt their private key, which would allow me to
ssh into any computer that recognizes my new laptop's public key.

## This cheat sheet

The `SConstruct` file in `deft` has a section which creates this
webpage.  It uses the html file index.html, which uses doc/index.html,
which uses doc/general-cheat-sheet.html, which is created from
doc/general-cheat-sheet.md.  Editing doc/general-cheat-sheet.md and
then running a scons from the deft directory will remake the websites'
html that can in turn can be read by a browser.

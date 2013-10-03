# Configuring your editor

You can use any text editor to work on deft, but there are some common
mistakes that you can avoid by appropriately configuring your text
editor.  This page will first outline those mistakes, and then explain
how to configure your editor to avoid them.

 * We prefer to use space characters rather than tab characters in
   deft.  This makes the indentation conisitent across all editors,
   and is generally nicer.  Your editor should *never* put tabs into a
   file.  If it does so, it is very hard for you to notice that this
   has happened.

 * If your editor supports automatic indentation, set it to indent by
   two spaces.

 * We would prefer never to end a line with spaces.  The single
   exception is that when a blank line occurs within source code, it
   is reasonable (but not necessary) for it to be indented just as
   much as the surrounding code.  Your editor should ideally show
   trailing whitespace so that you can avoid adding it.

The preceding rules are valuable because every change to white space
that is made could conflict with real changes made by another
developer to the same file.  This is a pain.

 * The final problem that can occur (and can be alleviated by an
   editor) is that if you have a file open while doing a git pull, you
   can end up accidentally saving the file, overwriting the version
   that you just pulled.  At a minimum, your editor should strongly
   warn you if a file is modified (e.g. by git) that you are editing.

Configuring emacs
=================

You configure emacs by editing the .emacs file in your home
directory.  The most relevant lines you will want to add are:

    (setq-default indent-tabs-mode nil)
    (setq-default show-trailing-whitespace t)

By default, emacs will warn you if you try to edit a buffer for a file
whose contents have changed on disk.  You could go one step further
and have emacs update the buffer automatically (if there were no
unsaved changes) by adding the following to your .emacs file:

    (global-auto-revert-mode)

Configuring kate
================

Kate is a reasonably nice text editor to use in working with deft (or
other programming projects.  It's pretty familiar to non-unixy folks,
but is less familiar to old-timers.  The following changes in kate's
configuration should make kate behave reasonably well.

 * Under "General" configuration, set kate to warn about files
   modified by foreign processes.
 * Under "Editing" set kate to insert spaces instead of tabs.
 * Under "Editing" set tab width to 2.
 * Under "Editing" set kate to highlight trailing spaces.
